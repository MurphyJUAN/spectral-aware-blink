# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Export a static data pack for the web viewer.
Usage: python export_web.py <cache_dir> <ra0> <dec0> <grid_n> <e0_dir> <out_json> <name>
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np

from core import load_cache, preprocess, reproject_frames, tan_grid
from bench import run_methods, pass_snr_maps

cache, ra0, dec0, n, e0dir, out, name = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
fr = load_cache(cache)
for e in fr:
    e["pre"] = preprocess(e)
w, ra, dec = tan_grid(ra0, dec0, n)
Y, V, L = reproject_frames(fr, ra, dec)
N = len(fr)
det = np.array([e["det"] for e in fr])
pas = np.array([e["pass"] for e in fr])
mjd = np.array([e["mjd"] for e in fr])
res = run_methods(Y.reshape(N, -1), V.reshape(N, -1), L.reshape(N, -1), det, pas, mjd, eta={"B2_detstack": 0.03, "C_spectral": 0.03, "B1_samecolor": 0.03},
                  methods=["B1_samecolor", "B2_detstack", "C_spectral"], alpha=100.0)
GROUPS = {"D123": [1, 2, 3], "D456": [4, 5, 6]}
img = {}
q = lambda a, s: np.round(np.nan_to_num(a) * s).astype(int).tolist()  # compact ints
for g, ds in GROUPS.items():
    for m in res:
        Z = pass_snr_maps(*res[m], pas, (n, n), dets=ds, det=det)
        img[f"{m}|{g}"] = [q(np.clip(Z[k], -20, 20), 10) for k in sorted(Z)]
    raws = []
    for k in np.unique(pas):
        mk = np.isin(det, ds) & (pas == k)
        wt = np.where(np.isfinite(Y[mk]) & np.isfinite(V[mk]), 1 / V[mk], 0)
        raw = (np.nan_to_num(Y[mk]) * wt).sum(0) / np.maximum(wt.sum(0), 1e-30)
        raws.append(q(raw, 1000))
    img[f"raw|{g}"] = raws
z = np.load(os.path.join(e0dir, "e0_phot.npz"))
sx, sy = w.wcs_world2pix(z["ra"], z["dec"], 0)
stars = []
for j in range(len(z["ra"])):
    ok = np.isfinite(z["F"][:, j]) & np.isfinite(z["LAM"][:, j])
    if ok.sum() < 20:
        continue
    stars.append(dict(id=j, x=round(float(sx[j]), 2), y=round(float(sy[j]), 2), ra=round(float(z["ra"][j]), 5), dec=round(float(z["dec"][j]), 5),
                      snr=round(float(z["snr"][j]), 1),
                      f=[round(float(v), 4) if np.isfinite(v) else None for v in z["F"][:, j]],
                      e=[round(float(v), 4) if np.isfinite(v) else None for v in z["E"][:, j]],
                      l=[round(float(v), 4) if np.isfinite(v) else None for v in z["LAM"][:, j]]))
# predicted W0855 positions per pass
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import astropy.units as u
c = SkyCoord(133.81890825 * u.deg, -7.247048 * u.deg, pm_ra_cosdec=-8123.7 * u.mas / u.yr, pm_dec=673.2 * u.mas / u.yr,
             distance=Distance(parallax=439.0 * u.mas), obstime=Time("J2000"))
marks = []
for k in np.unique(pas):
    p = c.apply_space_motion(Time(float(np.mean(mjd[pas == k])), format="mjd"))
    x, y = w.wcs_world2pix(p.ra.deg, p.dec.deg, 0)
    marks.append([round(float(x), 2), round(float(y), 2)])
pack = dict(name=name, ra0=ra0, dec0=dec0, n=n, pix_arcsec=6.15, mjd=[round(float(t), 4) for t in z["mjd"]], det=z["det"].tolist(), pas=z["pas"].tolist(),
            pass_dates=[Time(float(np.mean(mjd[pas == k])), format="mjd").iso[:10] for k in np.unique(pas)],
            img=img, stars=stars, mover={"name": "WISE 0855-0714", "xy": marks},
            scale={"Z": 10, "raw": 1000})
json.dump(pack, open(out, "w"), separators=(",", ":"))
print("wrote", out, os.path.getsize(out) // 1024, "KB", len(stars), "stars")
