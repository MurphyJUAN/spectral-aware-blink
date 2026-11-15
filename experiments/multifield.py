# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Multi-field summary of the source-level rules. Usage: python multifield.py <results_dir> <out_json> TAG [TAG...]"""
import json
import os
import sys

import numpy as np
import pandas as pd

R, out, tags = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
from astropy.timeseries import LombScargle
from core import fit_predict


def e9(npz):
    z = np.load(npz); F, E, LAM, det, mjd = z["F"], z["E"], z["LAM"], z["det"], z["mjd"]
    fr = {"raw": [], "sed": []}
    for d in range(1, 7):
        m = det == d
        for j in range(F.shape[1]):
            f, e, l, t = F[m, j], E[m, j], LAM[m, j], mjd[m]
            ok = np.isfinite(f) & np.isfinite(l) & (f > 0)
            if ok.sum() < 20 or np.median(f[ok] / e[ok]) < 30:
                continue
            f, e, l, t = f[ok], e[ok], l[ok], t[ok]
            pr, _ = fit_predict(f[:, None], (e ** 2)[:, None], l[:, None], np.ones(len(f), bool), np.ones(len(f), bool), nk=6, alpha=100.0)
            for key, y in [("raw", f / np.median(f)), ("sed", f / pr[:, 0])]:
                ls = LombScargle(t, y, np.sqrt((e / f) ** 2 + 0.03 ** 2))
                fq, pw = ls.autopower(minimum_frequency=1 / 200., maximum_frequency=2., samples_per_peak=10)
                fr[key].append(ls.false_alarm_probability(pw.max(), minimum_frequency=1 / 200., maximum_frequency=2.))
    return len(fr["raw"]), float(np.mean(np.array(fr["raw"]) < 0.01)), float(np.mean(np.array(fr["sed"]) < 0.01))
rms = lambda v: float(np.sqrt(np.nanmean(np.asarray(v, float) ** 2)))
mad = lambda v: float(1.4826 * np.nanmedian(np.abs(np.asarray(v, float))))
rows = []
for t in tags:
    p6 = os.path.join(R, f"e6b_f{t}", "e6.json")
    p7 = os.path.join(R, f"e7b_f{t}", "e7.csv")
    p0 = os.path.join(R, f"e0_f{t}", "e0.json")
    if not (os.path.exists(p6) and os.path.exists(p7)):
        print("missing", t)
        continue
    e0 = json.load(open(p0))
    h1 = pd.DataFrame(e0["h1"])
    al = pd.DataFrame(e0["alias"])
    df = pd.DataFrame(json.load(open(p6)))
    n_bad = int((df.beta.abs() >= 10).sum())
    df = df[df.beta.abs() < 10]  # sanity cut, same for all fields: drop failed SED fits (|beta|>=10 is unphysical)
    df["obs_zp"] = df.obs - df.groupby(["det", "pas"]).obs.transform("median")
    df["bf"] = df.groupby("det").beta.transform("median")
    df["x"] = (df.beta - df.bf) * (df.lin / df.beta.replace(0, np.nan))
    df["ad"] = (df.beta - df.bf).abs()
    pec = df.ad > df.ad.quantile(0.75)
    e7 = pd.read_csv(p7)
    n7 = e7[e7.A == 0].copy()
    n7["pf"] = n7.nb1 / n7.nk
    lo = n7.pf < 0.34
    thr = {m: np.nanpercentile(n7[m].abs(), 95) for m in ["B2", "B1", "C"]}
    a12 = e7[e7.A == 0.12]
    n9, p_raw, p_sed = e9(os.path.join(R, f"e0_f{t}", "e0_phot.npz"))
    rows.append(dict(
        field=t, n_star_det=len(h1), n_frames=int(len(np.load(os.path.join(R, f"e0_f{t}", "e0_phot.npz"))["mjd"])),
        e9_n=n9, e9_periodic_raw=p_raw, e9_periodic_sed=p_sed,
        R2_lam_bright=float(h1[h1.snr > 30].r2_lam.median()), R2_time_bright=float(h1[h1.snr > 30].r2_time.median()),
        chi2_const_bright=float(h1[h1.snr > 30].chi2_const.median()), chi2_spline_bright=float(h1[h1.snr > 30].chi2_spline.median()),
        absrho_D3D6=float(al[al.det.isin([3, 6])].rho.abs().median()) if len(al) else np.nan,
        n_e6=len(df), n_e6_cut=n_bad, B2_raw=rms(df.obs), B2_zp=rms(df.obs_zp), C=rms(df.obs - df.pred), lin=rms(df.obs - df.lin),
        r_pred=float(df[["obs", "pred"]].corr().iloc[0, 1]), r_relcol=float(df[["x", "obs_zp"]].corr().iloc[0, 1]),
        B2zp_peculiar=rms(df.obs_zp[pec]), C_peculiar=rms((df.obs - df.pred)[pec]),
        var_sig_B2=mad(n7.B2), var_sig_B1=mad(n7.B1), var_sig_C=mad(n7.C), B1_undef=float(n7.B1.isna().mean()),
        CoverB1_sparse=mad(n7.C[lo]) / mad(n7.B1[lo]) if lo.sum() > 5 else np.nan, n_sparse=int(lo.sum()),
        compl12_B2=float((a12.B2 > thr["B2"]).mean()), compl12_B1=float((a12.B1 > thr["B1"]).mean()), compl12_C=float((a12.C > thr["C"]).mean()),
    ))
t = pd.DataFrame(rows).set_index("field")
pd.set_option("display.width", 250)
print(t.T.round(3))
t.to_json(out, indent=1)
