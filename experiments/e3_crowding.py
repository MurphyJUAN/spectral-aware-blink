# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""E3b: completeness vs local stationary-sky brightness at the injected position (crowding rule).
Usage: python e3_crowding.py <e3_dir> <grid_npz_from_e1> [fap]
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

d, gpath = sys.argv[1], sys.argv[2]
fap = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(d, "inj_*.csv")))])
nulls = {}
for f in sorted(glob.glob(os.path.join(d, "null_*.npz"))):
    z = np.load(f)
    for k in z.files:
        nulls.setdefault(k, []).append(z[k])
thr = {k: float(np.nanpercentile(np.concatenate(v), 100 * (1 - fap))) for k, v in nulls.items()}
snr = np.load(gpath)["snr"]
n = snr.shape[0]
yy, xx = np.mgrid[0:n, 0:n]


def local(gx, gy, r=1.5):
    m = np.hypot(xx - gx, yy - gy) <= r
    return float(np.nanmax(snr[m])) if m.any() else np.nan


df["sky"] = [local(a, b) for a, b in zip(df.gx, df.gy)]
bins = [-np.inf, 10, 50, 300, np.inf]
labels = ["empty(<10)", "faint star(10-50)", "star(50-300)", "bright star(>300)"]
df["skybin"] = pd.cut(df.sky, bins, labels=labels)
METH = sorted({k.split("|")[0] for k in thr})
df = df[(df.pm >= 1.0) & (df.amp >= 0.02)]
out = {}
pd.set_option("display.width", 250)
for g in ["all", "long"]:
    t = pd.DataFrame({m: (df[f"{m}|{g}"] > thr[f"{m}|{g}"]).groupby(df.skybin, observed=False).mean() for m in METH})
    t["N"] = df.groupby("skybin", observed=False).size()
    print(f"== group {g}: completeness (pm>=1, amp>=0.02) by local stationary-sky SNR at injected position")
    print(t.round(3))
    out[g] = t.to_dict()
print("fraction of sky area per bin:", {l: round(float(np.mean(pd.cut(pd.Series(snr.ravel()), bins, labels=labels) == l)), 3) for l in labels})
json.dump({g: {m: {str(k): v for k, v in d_.items()} for m, d_ in t_.items()} for g, t_ in out.items()}, open(os.path.join(d, "crowding.json"), "w"), indent=1)
