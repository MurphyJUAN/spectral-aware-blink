# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Demo GIF: per-pass blink of raw stacks vs residual SNR maps (B2, C) around WISE 0855 in D5+D6.
Usage: python blink_gif.py <cache_dir> <ra0> <dec0> <grid_n> <out_gif>
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import astropy.units as u

from core import load_cache, preprocess, reproject_frames, tan_grid
from bench import run_methods, pass_snr_maps

cache, ra0, dec0, n, outgif = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
fr = load_cache(cache)
for e in fr:
    e["pre"] = preprocess(e)
w, ra, dec = tan_grid(ra0, dec0, n)
Y, V, L = reproject_frames(fr, ra, dec)
N = len(fr)
det = np.array([e["det"] for e in fr])
pas = np.array([e["pass"] for e in fr])
mjd = np.array([e["mjd"] for e in fr])
res = run_methods(Y.reshape(N, -1), V.reshape(N, -1), L.reshape(N, -1), det, pas, mjd, eta={"B2_detstack": 0.03, "C_spectral": 0.03},
                  methods=["B2_detstack", "C_spectral"], alpha=100.0)
ds = [5, 6]
Z = {m: pass_snr_maps(*res[m], pas, (n, n), dets=ds, det=det) for m in res}
c = SkyCoord(133.81890825 * u.deg, -7.247048 * u.deg, pm_ra_cosdec=-8123.7 * u.mas / u.yr, pm_dec=673.2 * u.mas / u.yr,
             distance=Distance(parallax=439.0 * u.mas), obstime=Time("J2000"))
DIV = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#e8e7e3", "#f2a27f", "#b8461b"])
frames = []
m56 = np.isin(det, ds)
vmax = None
for k in np.unique(pas):
    mk = m56 & (pas == k)
    wts = np.where(np.isfinite(Y[mk]) & np.isfinite(V[mk]), 1 / V[mk], 0)
    raw = (np.nan_to_num(Y[mk]) * wts).sum(0) / np.maximum(wts.sum(0), 1e-30)
    if vmax is None:
        vmax = np.nanpercentile(raw, 99.5)
    tk = float(np.mean(mjd[pas == k]))
    p = c.apply_space_motion(Time(tk, format="mjd"))
    x, y = w.wcs_world2pix(p.ra.deg, p.dec.deg, 0)
    lam = np.nanmedian(L[mk][:, n // 2, n // 2])
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.9), facecolor="#fcfcfb")
    ax[0].imshow(np.arcsinh(raw / (0.05 * vmax)), origin="lower", cmap="Greys_r", vmin=0, vmax=np.arcsinh(20))
    ax[0].set_title(f"naive pass stack D5+D6\n(⟨λ⟩ at centre {lam:.2f} µm)", fontsize=9)
    ax[1].imshow(Z["B2_detstack"][int(k)], origin="lower", cmap=DIV, vmin=-5, vmax=5)
    ax[1].set_title("B2 detector-stack residual (SNR)", fontsize=9)
    ax[2].imshow(Z["C_spectral"][int(k)], origin="lower", cmap=DIV, vmin=-5, vmax=5)
    ax[2].set_title("C spectral-LOPO residual (SNR)", fontsize=9)
    for a in ax:
        a.add_patch(plt.Circle((x, y), 3, fill=False, color="#eb6834", lw=1.5))
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"Pass {k+1}  ({Time(tk, format='mjd').iso[:10]})  —  circle: WISE 0855 predicted position", fontsize=10)
    fig.tight_layout()
    fn = f"/tmp/blink_{k}.png"
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    frames.append(Image.open(fn).convert("P", palette=Image.ADAPTIVE))
frames[0].save(outgif, save_all=True, append_images=frames[1:], duration=900, loop=0)
print("wrote", outgif)
