# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""E6: is the pass-to-pass 'variability' of detector stacks (B2) predicted by the SED x wavelength sampling?
Observed: ln(F_k / F_all) of per-pass inverse-variance stacks of stationary stars.
Predicted: ln(mean_{i in k} s(lam_i) / mean_all s(lam_i)) with s = spline SED fitted to ALL OTHER passes (LOPO),
and the 1st-order law beta * (<ln lam>_k - <ln lam>_all).
Usage: python e6_b2law.py <e0_phot.npz> <outdir>
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np

from core import fit_predict

z = np.load(sys.argv[1])
out = sys.argv[2]
os.makedirs(out, exist_ok=True)
F, E, LAM, det, pas = z["F"], z["E"], z["LAM"], z["det"], z["pas"]
rows = []
for d in range(1, 7):
    m = det == d
    for j in range(F.shape[1]):
        f, e, l, p = F[m, j], E[m, j], LAM[m, j], pas[m]
        ok = np.isfinite(f) & np.isfinite(l) & (f > 0)
        if ok.sum() < 15:
            continue
        y, s, x, pp = f[ok], e[ok], l[ok], p[ok]
        snr = np.median(y / s)
        if snr < 15:
            continue
        w = 1 / s ** 2
        Fall = np.sum(w * y) / np.sum(w)
        # full-data spline (for the 'all' reference) and slope beta
        pr_all, _ = fit_predict(y[:, None], (s ** 2)[:, None], x[:, None], np.ones(len(y), bool), np.ones(len(y), bool), nk=6, alpha=float(os.environ.get("SALPHA", "100")))
        beta = np.polyfit(np.log(x), np.log(np.clip(pr_all[:, 0], 1e-9, None)), 1)[0]
        for k in np.unique(pp):
            te = pp == k
            if te.sum() < 3 or (~te).sum() < 8:
                continue
            Fk = np.sum(w[te] * y[te]) / np.sum(w[te])
            obs = np.log(Fk / Fall)
            # LOPO prediction: SED learned without pass k, evaluated at pass-k wavelengths vs all wavelengths
            pr, _ = fit_predict(y[:, None], (s ** 2)[:, None], x[:, None], ~te, np.ones(len(y), bool), nk=6, alpha=float(os.environ.get("SALPHA", "100")))
            sk = pr[te, 0]
            pred = np.log(np.sum(w[te] * sk) / np.sum(w[te]) / (np.sum(w * pr[:, 0]) / np.sum(w)))
            lin = beta * (np.sum(w[te] * np.log(x[te])) / np.sum(w[te]) - np.sum(w * np.log(x)) / np.sum(w))
            err = np.sqrt(1 / np.sum(w[te])) / Fk
            rows.append(dict(det=d, src=int(j), pas=int(k), n=int(te.sum()), snr=float(snr), beta=float(beta), obs=float(obs), pred=float(pred), lin=float(lin), err=float(err)))
json.dump(rows, open(os.path.join(out, "e6.json"), "w"), indent=0)
import pandas as pd
df = pd.DataFrame(rows)
df["resid"] = df.obs - df.pred
print("N", len(df), "stars", df.groupby(["det", "src"]).ngroups)
print("rms obs (B2 false variability)      :", round(float(np.sqrt(np.mean(df.obs ** 2))), 4))
print("rms obs - pred (after spectral corr):", round(float(np.sqrt(np.mean(df.resid ** 2))), 4))
print("rms obs - lin (1st-order beta law)  :", round(float(np.sqrt(np.mean((df.obs - df.lin) ** 2))), 4))
print("median stat err                     :", round(float(df.err.median()), 4))
print("corr(obs,pred)", round(float(np.corrcoef(df.obs, df.pred)[0, 1]), 3), " slope", round(float(np.polyfit(df.pred, df.obs, 1)[0]), 3))
print(df.groupby("det").apply(lambda g: pd.Series(dict(rms_obs=np.sqrt(np.mean(g.obs ** 2)), rms_res=np.sqrt(np.mean(g.resid ** 2)), med_err=g.err.median(), n=len(g)))).round(4))
