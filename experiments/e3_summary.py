# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Aggregate E3 runs: completeness at equal false-alarm probability; flux at 50% completeness per method.
Usage: python e3_summary.py <e3_dir> [fap=0.01] [min_pm=1.0]
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
fap = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
min_pm = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(d, "inj_*.csv")))])
nulls = {}
for f in sorted(glob.glob(os.path.join(d, "null_*.npz"))):
    z = np.load(f)
    for k in z.files:
        nulls.setdefault(k, []).append(z[k])
nulls = {k: np.concatenate(v) for k, v in nulls.items()}
METH = sorted({k.split("|")[0] for k in nulls})
GROUPS = sorted({k.split("|")[1] for k in nulls})
thr = {k: float(np.nanpercentile(v, 100 * (1 - fap))) for k, v in nulls.items()}


def amp50(c):
    a, f = np.array(c.index, float), np.array(c.values, float)
    if f.max() < 0.5:
        return np.inf
    i = np.argmax(f >= 0.5)
    if i == 0:
        return a[0]
    # log-linear interpolation
    la = np.interp(0.5, [f[i - 1], f[i]], np.log([a[i - 1], a[i]]))
    return float(np.exp(la))


res = {"thr": thr, "n_inj": int(len(df)), "fap": fap, "min_pm": min_pm, "amp50": {}, "compl": {}}
sub = df[df.pm >= min_pm]
pd.set_option("display.width", 250)
for g in GROUPS:
    print(f"=== group {g} | thresholds:", {m: round(thr[f'{m}|{g}'], 2) for m in METH})
    rows = []
    for kind in ["star", "tdwarf", "ydwarf"]:
        s = sub[sub.kind == kind]
        row = {"kind": kind}
        for m in METH:
            key = f"{m}|{g}"
            c = (s[key] > thr[key]).groupby(s.amp).mean()
            res["compl"][f"{g}|{kind}|{m}"] = c.to_dict()
            row[m] = amp50(c)
            res["amp50"][f"{g}|{kind}|{m}"] = row[m]
        rows.append(row)
    t = pd.DataFrame(rows).set_index("kind")
    print("flux at 50% completeness [MJy/sr] (lower = better):")
    print(t.round(4))
    best_base = t[[m for m in METH if m.startswith("B")]].min(axis=1)
    print("depth gain of C over best baseline (x):", (best_base / t["C_spectral"]).round(2).to_dict())
print("by pm (all, all kinds, mean completeness):")
print(pd.DataFrame({m: (df[f"{m}|all"] > thr[f"{m}|all"]).groupby(df.pm).mean() for m in METH}).round(3))
json.dump(res, open(os.path.join(d, "summary.json"), "w"), indent=1, default=float)
