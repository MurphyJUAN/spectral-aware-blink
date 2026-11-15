# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""E0: is per-frame flux variation of stationary stars explained by wavelength (LVF) rather than time?
Also: within-pass time<->wavelength aliasing at fixed sky position.
Usage: python e0_problem.py <cache_dir> <ra0> <dec0> <grid_n> <outdir>
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np
from scipy.stats import spearmanr

from core import load_cache, fit_predict
from phot import detect_sources, forced_phot

cache, ra0, dec0, n, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
os.makedirs(out, exist_ok=True)
fr = load_cache(cache)
print("frames", len(fr))
src, (w, ra, dec, Y, V, L, snr) = detect_sources(fr, ra0, dec0, n, nsig=10)
print("sources", len(src["ra"]))
F, E, LAM, XX, YY = forced_phot(fr, src["ra"], src["dec"])
det = np.array([e["det"] for e in fr])
mjd = np.array([e["mjd"] for e in fr])
pas = np.array([e["pass"] for e in fr])
np.savez(os.path.join(out, "e0_phot.npz"), F=F, E=E, LAM=LAM, det=det, mjd=mjd, pas=pas, ra=src["ra"], dec=src["dec"], snr=src["snr"])

# --- time<->lambda aliasing within a pass (at the grid centre) ---
alias = []
lamc = L[:, n // 2, n // 2]
for d in range(1, 7):
    for p in np.unique(pas):
        m = (det == d) & (pas == p)
        if m.sum() >= 5:
            rho = spearmanr(mjd[m], lamc[m]).statistic
            alias.append(dict(det=d, pas=int(p), n=int(m.sum()), rho=float(rho),
                              lam_min=float(lamc[m].min()), lam_max=float(lamc[m].max()),
                              span_days=float(mjd[m].max() - mjd[m].min())))
print("alias", json.dumps(alias, indent=0))

# --- H1: variance explained by lambda-spline vs per-pass (time) means ---
res = []
for d in range(1, 7):
    m = det == d
    if m.sum() < 10:
        continue
    f, e, l, p = F[m], E[m], LAM[m], pas[m]
    for j in range(F.shape[1]):
        ok = np.isfinite(f[:, j]) & np.isfinite(l[:, j])
        if ok.sum() < 12:
            continue
        y, s, x, pp = f[ok, j], e[ok, j], l[ok, j], p[ok]
        snr_med = np.median(y / s)
        if snr_med < 5:
            continue
        tot = np.var(y)
        # lambda model: spline (6 knots) fitted in-sample
        pred, _ = fit_predict(y[:, None], (s ** 2)[:, None], x[:, None], np.ones(len(y), bool), np.ones(len(y), bool), model="spline", nk=6, alpha=0.1, iters=1)
        r2l = 1 - np.var(y - pred[:, 0]) / tot
        # time model: per-pass mean
        tp = np.zeros_like(y)
        for q in np.unique(pp):
            tp[pp == q] = np.mean(y[pp == q])
        r2t = 1 - np.var(y - tp) / tot
        # LOPO spline chi2 vs LOPO const chi2
        chis = {}
        for mod in ["const", "spline"]:
            r = []
            for q in np.unique(pp):
                te = pp == q
                if (~te).sum() < 6:
                    continue
                pr, pv = fit_predict(y[:, None], (s ** 2)[:, None], x[:, None], ~te, te, model=mod, nk=6, alpha=0.1, iters=2)
                r.append((y[te] - pr[:, 0]) / np.sqrt(s[te] ** 2 + pv[:, 0]))
            r = np.concatenate(r) if r else np.array([np.nan])
            chis[mod] = float(np.mean(r ** 2))
        res.append(dict(det=d, src=j, n=int(ok.sum()), snr=float(snr_med), r2_lam=float(r2l), r2_time=float(r2t),
                        frac_rms=float(np.std(y) / np.mean(y)), chi2_const=chis["const"], chi2_spline=chis["spline"],
                        lam_span=float(x.max() - x.min())))
json.dump(dict(alias=alias, h1=res), open(os.path.join(out, "e0.json"), "w"), indent=1)
import pandas as pd
df = pd.DataFrame(res)
if len(df):
    print(df.groupby("det")[["r2_lam", "r2_time", "frac_rms", "chi2_const", "chi2_spline", "snr"]].median().round(3))
    print("N star-det", len(df))
