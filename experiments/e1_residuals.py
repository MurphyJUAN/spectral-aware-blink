# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""E1: residual cleanliness of stationary sky (chi2_nu by pixel class) for B0/B1/B2/C, eta calibration,
false-positive rate in per-pass SNR maps.
Usage: python e1_residuals.py <cache_dir> <ra0> <dec0> <grid_n> <outdir>
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np

from core import load_cache, preprocess, reproject_frames, tan_grid
from bench import run_methods, pass_snr_maps, chi2_stats, METHODS

cache, ra0, dec0, n, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
os.makedirs(out, exist_ok=True)
t0 = time.time()
fr = load_cache(cache)
for e in fr:
    e["pre"] = preprocess(e)
w, ra, dec = tan_grid(ra0, dec0, n)
if os.environ.get("REGISTER", "0") == "1":
    from phot import detect_sources, register
    _src, _ = detect_sources(fr, ra0, dec0, n, nsig=30)
    _off = register(fr, _src["ra"], _src["dec"])
    print("registered on", len(_src["ra"]), "stars; shift rms dx=%.3f dy=%.3f px" % (np.std(_off[:, 0]), np.std(_off[:, 1])), flush=True)
Y, V, L = reproject_frames(fr, ra, dec)
det = np.array([e["det"] for e in fr])
pas = np.array([e["pass"] for e in fr])
mjd = np.array([e["mjd"] for e in fr])
N = len(fr)
Yf, Vf, Lf = Y.reshape(N, -1), V.reshape(N, -1), L.reshape(N, -1)
PH = reproject_frames.last_phase.reshape(N, -1, 2)
print("loaded", N, "frames", time.time() - t0)

# pixel classes from the per-detector median SNR stack
snr_stack = np.zeros(n * n)
for d in np.unique(det):
    m = det == d
    med = np.nanmedian(Yf[m], 0)
    noise = np.sqrt(np.nanmedian(Vf[m], 0) / m.sum())
    snr_stack = np.maximum(snr_stack, np.nan_to_num(med / noise))
cls = np.where(snr_stack > 50, "bright", np.where(snr_stack > 10, "faint", "bkg"))

res = run_methods(Yf, Vf, Lf, det, pas, mjd, PH=PH, shape=(n, n))
summary = {}
for meth in METHODS:
    R, S = res[meth]
    c2, use = chi2_stats(R, S)
    summary[meth] = {c: dict(chi2_med=float(np.nanmedian(c2[cls == c])), use=float(np.nanmean(use[cls == c])), npix=int((cls == c).sum())) for c in ["bkg", "faint", "bright"]}

# eta: fractional floor so that median chi2 in bright pixels ~1 (per method)
eta = {}
for meth in METHODS:
    R, S = res[meth]
    pred = Yf - R
    best = None
    for e_ in [0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.065, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]:
        z2 = R ** 2 / (S ** 2 + (e_ * pred) ** 2)
        c = np.nanmedian(np.nanmean(z2, 0)[cls == "bright"])
        if c <= 1.5:
            best = e_
            break
    eta[meth] = best if best is not None else 0.5
print("eta", eta)
res2 = run_methods(Yf, Vf, Lf, det, pas, mjd, eta=eta, PH=PH, shape=(n, n))
fp = {}
for meth in METHODS:
    R, S = res2[meth]
    c2, use = chi2_stats(R, S)
    summary[meth]["eta"] = eta[meth]
    for c in ["bkg", "faint", "bright"]:
        summary[meth][c]["chi2_med_eta"] = float(np.nanmedian(c2[cls == c]))
    Z = pass_snr_maps(R, S, pas, (n, n))
    zall = np.concatenate([z[3:-3, 3:-3].ravel() for z in Z.values()])
    zall = zall[np.isfinite(zall)]
    fp[meth] = dict(n_gt5=int((zall > 5).sum()), n_lt_m5=int((zall < -5).sum()), std=float(np.std(zall)),
                    mad=float(1.4826 * np.median(np.abs(zall))), npix=int(zall.size))
    np.save(os.path.join(out, f"Z_{meth}.npy"), np.array([Z[k] for k in sorted(Z)]))
summary["fp"] = fp
json.dump(summary, open(os.path.join(out, "e1.json"), "w"), indent=1)
np.savez(os.path.join(out, "grid.npz"), Y=Y, L=L, det=det, pas=pas, mjd=mjd, cls=cls.reshape(n, n), snr=snr_stack.reshape(n, n))
for meth in METHODS:
    s = summary[meth]
    print(f"{meth:14s} eta={s['eta']:.3f} | " + " | ".join(f"{c}: chi2={s[c]['chi2_med']:.2f}->{s[c]['chi2_med_eta']:.2f} use={s[c]['use']:.2f}" for c in ["bkg", "faint", "bright"]) + f" | Z>5: {fp[meth]['n_gt5']} Z<-5: {fp[meth]['n_lt_m5']} std={fp[meth]['std']:.2f}")
print("time", time.time() - t0)
