# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""E7: source-level variability injection-recovery on REAL stationary-star photometry.
For each (star, detector, pass k): multiply pass-k fluxes by (1+A) and test whether the pass-k brightening is
detected at equal false-alarm probability (thresholds from the un-injected statistic distribution).
Statistics (all signed, in ln flux):
  B2: ln(stack_k / stack_others)                           (detector stack, lambda-blind)
  B1: mean over pass-k frames of ln(F_i / F_partner) with |dlam/lam|<tau partners in other passes
  C : ln(stack_k / stack of LOPO SED prediction at pass-k wavelengths)  (spectral model of other passes)
Each statistic is normalised by its empirical (MAD) scatter -> thresholds at FAP.
Usage: python e7_varinj.py <e0_phot.npz> <outdir> [tau]
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np

from core import fit_predict

z = np.load(sys.argv[1])
out = sys.argv[2]
tau = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
os.makedirs(out, exist_ok=True)
F, E, LAM, det, pas = z["F"], z["E"], z["LAM"], z["det"], z["pas"]
KEEP = float(os.environ.get("KEEP", "1"))
if KEEP < 1:  # thin exposures (simulate sparse coverage / early mission), same subset for all methods
    rk = np.random.default_rng(int(os.environ.get("SEED", "0")))
    keep = rk.uniform(size=len(det)) < KEEP
    F, E, LAM, det, pas = F[keep], E[keep], LAM[keep], det[keep], pas[keep]
    print("kept", keep.sum(), "of", len(keep), "frames")
AMPS = [0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3]


def stats(y, s, x, pp, k, A):
    y = y.copy()
    te = pp == k
    y[te] *= 1 + A
    w = 1 / s ** 2
    # B2
    b2 = np.log((np.sum(w[te] * y[te]) / np.sum(w[te])) / (np.sum(w[~te] * y[~te]) / np.sum(w[~te])))
    # B1
    r = []
    for i in np.nonzero(te)[0]:
        mm = (~te) & (np.abs(np.log(x / x[i])) < tau)
        if mm.any():
            ref = np.sum(w[mm] * y[mm]) / np.sum(w[mm])
            if ref > 0 and y[i] > 0:
                r.append(np.log(y[i] / ref))
    b1 = float(np.mean(r)) if r else np.nan
    nb1 = len(r)
    # C
    pr, _ = fit_predict(y[:, None], (s ** 2)[:, None], x[:, None], ~te, te, nk=6, alpha=float(os.environ.get("SALPHA", "100")))
    pk = pr[:, 0]
    c = np.log((np.sum(w[te] * y[te]) / np.sum(w[te])) / (np.sum(w[te] * pk) / np.sum(w[te])))
    return b2, b1, c, nb1


rows = []
for d in range(1, 7):
    m = det == d
    for j in range(F.shape[1]):
        f, e, l, p = F[m, j], E[m, j], LAM[m, j], pas[m]
        ok = np.isfinite(f) & np.isfinite(l) & (f > 0)
        if ok.sum() < 15:
            continue
        y, s, x, pp = f[ok], e[ok], l[ok], p[ok]
        snr = float(np.median(y / s))
        if snr < 10:
            continue
        for k in np.unique(pp):
            if (pp == k).sum() < 3 or (pp != k).sum() < 8:
                continue
            for A in AMPS:
                b2, b1, c, nb1 = stats(y, s, x, pp, k, A)
                rows.append(dict(det=d, src=int(j), pas=int(k), snr=snr, A=A, B2=b2, B1=b1, C=c, nb1=nb1, nk=int((pp == k).sum())))
import pandas as pd
df = pd.DataFrame(rows)
df.to_csv(os.path.join(out, "e7.csv"), index=False)
null = df[df.A == 0]
res = {}
pd.set_option("display.width", 200)
for fap in [0.01, 0.05]:
    thr = {m: float(np.nanpercentile(null[m].abs(), 100 * (1 - fap))) for m in ["B2", "B1", "C"]}
    t = pd.DataFrame({m: (df[m] > thr[m]).groupby(df.A).mean() for m in ["B2", "B1", "C"]})
    # B1 is undefined when a pass has no same-colour partner -> counts as non-detection (data loss)
    print(f"FAP={fap}: |stat| thresholds (ln flux):", {k: round(v, 4) for k, v in thr.items()})
    print("completeness vs injected brightening A of one pass:")
    print(t.round(3))
    res[str(fap)] = dict(thr=thr, compl=t.to_dict())
for lo, hi in [(10, 30), (30, 100), (100, 1e9)]:
    s_ = null[(null.snr >= lo) & (null.snr < hi)]
    print(f"null scatter (1.4826*MAD) snr {lo}-{hi}: N={len(s_)}", {m: round(float(1.4826 * np.nanmedian(np.abs(s_[m]))), 4) for m in ["B2", "B1", "C"]},
          "B1 undefined frac:", round(float(s_.B1.isna().mean()), 3))
json.dump(res, open(os.path.join(out, "e7.json"), "w"), indent=1, default=float)
