# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""E9: spurious periodicity in naive SPHEREx light curves of stationary stars (consequence of R1).
For each bright star x detector: Lomb-Scargle on (a) naive relative flux, (b) flux / smooth SED (fit on all data).
Count fraction flagged 'periodic' at FAP<1%. Usage: python e9_periodogram.py <e0_phot.npz>"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np
from astropy.timeseries import LombScargle
from core import fit_predict
z = np.load(sys.argv[1]); F, E, LAM, det, mjd = z["F"], z["E"], z["LAM"], z["det"], z["mjd"]
out = {"raw": [], "sed": []}; best = []
for d in range(1, 7):
    m = det == d
    for j in range(F.shape[1]):
        f, e, l, t = F[m, j], E[m, j], LAM[m, j], mjd[m]
        ok = np.isfinite(f) & np.isfinite(l) & (f > 0)
        if ok.sum() < 20 or np.median(f[ok] / e[ok]) < 30: continue
        f, e, l, t = f[ok], e[ok], l[ok], t[ok]
        pr, _ = fit_predict(f[:, None], (e ** 2)[:, None], l[:, None], np.ones(len(f), bool), np.ones(len(f), bool), nk=6, alpha=100.0)
        floor = 0.03
        for key, y in [("raw", f / np.median(f)), ("sed", f / pr[:, 0])]:
            dy = np.sqrt((e / f) ** 2 + floor ** 2)
            ls = LombScargle(t, y, dy)
            freq, pw = ls.autopower(minimum_frequency=1 / 200., maximum_frequency=1 / 0.5, samples_per_peak=10)
            fap = ls.false_alarm_probability(pw.max(), minimum_frequency=1 / 200., maximum_frequency=1 / 0.5)
            out[key].append(fap)
            if key == "raw": best.append(1 / freq[np.argmax(pw)])
for k, v in out.items():
    v = np.array(v); print(f"{k:4s}: N={len(v)}  fraction with FAP<1%: {np.mean(v < 0.01):.2f}   FAP<1e-6: {np.mean(v < 1e-6):.2f}")
b = np.array(best); print("naive best periods [d]: median %.1f, IQR %.1f-%.1f" % (np.median(b), *np.percentile(b, [25, 75])))
