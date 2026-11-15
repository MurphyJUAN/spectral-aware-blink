# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""E3: injection-recovery of moving point sources. Same preprocessing + same detector for all methods;
thresholds set per method from the null (no-injection) maps at equal false-alarm probability.
Usage: python e3_inject.py <cache_dir> <ra0> <dec0> <grid_n> <outdir> <n_real> [seed] [n_src]
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np

from core import load_cache, preprocess, reproject_frames, tan_grid
from bench import run_methods, pass_snr_maps, inject, METHODS

cache, ra0, dec0, n, out, nreal = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
seed = int(sys.argv[7]) if len(sys.argv) > 7 else 0
nsrc = int(sys.argv[8]) if len(sys.argv) > 8 else 8
ETA = float(os.environ.get("ETA", "0.03"))
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(seed)
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
det = np.array([e["det"] for e in fr])
pas = np.array([e["pass"] for e in fr])
mjd = np.array([e["mjd"] for e in fr])
t_ref = float(np.mean(mjd))
passes = np.unique(pas)
tpass = {int(k): float(np.mean(mjd[pas == k])) for k in passes}
N = len(fr)
eta = {m: ETA for m in METHODS}
GROUPS = {"all": None, "long": [4, 5, 6], "short": [1, 2, 3]}


def maps_for(frames_key):
    for e in fr:
        e["cur"] = e[frames_key]
    saved = [e["pre"] for e in fr]
    for e in fr:
        e["pre"] = e["cur"]
    Y, V, L = reproject_frames(fr, ra, dec)
    for e, s in zip(fr, saved):
        e["pre"] = s
    PH = reproject_frames.last_phase.reshape(N, -1, 2)
    res = run_methods(Y.reshape(N, -1), V.reshape(N, -1), L.reshape(N, -1), det, pas, mjd, eta=eta, PH=PH, shape=(n, n))
    Z = {}
    for meth in METHODS:
        R, S = res[meth]
        for g, ds in GROUPS.items():
            Z[(meth, g)] = pass_snr_maps(R, S, pas, (n, n), dets=ds, det=det)
    return Z


def local_max(z, x, y, r=1.0):
    yy, xx = np.mgrid[0:n, 0:n]
    m = (np.hypot(xx - x, yy - y) <= r) & np.isfinite(z)
    return float(np.max(z[m])) if m.any() else np.nan


# ---- null maps and per-method thresholds at equal false-alarm probability ----
Z0 = maps_for("pre")
null = {}
for key, Zk in Z0.items():
    vals = []
    for k, z in Zk.items():
        for _ in range(400):
            x, y = rng.uniform(5, n - 6, 2)
            vals.append(local_max(z, x, y))
    null["|".join(key)] = np.array(vals)
print("null done", time.time() - t0, flush=True)

AMPS = np.array([0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32])  # MJy/sr peak at sed==1
PMS = np.array([0.5, 1.0, 2.0, 4.0, 10.0])  # px/yr
KINDS = ["star", "tdwarf", "ydwarf"]
rows = []
for r in range(nreal):
    srcs = []
    for s in range(nsrc):
        # keep sources apart; positions in inner region
        for _ in range(50):
            x, y = rng.uniform(8, n - 9, 2)
            if all(np.hypot(x - q["gx"], y - q["gy"]) > 6 for q in srcs):
                break
        sra, sdec = w.wcs_pix2world(x, y, 0)
        ang = rng.uniform(0, 2 * np.pi)
        pm = rng.choice(PMS)
        srcs.append(dict(gx=x, gy=y, ra=float(sra), dec=float(sdec), pm=float(pm), pmra=float(pm * np.cos(ang)), pmdec=float(pm * np.sin(ang)),
                         amp=float(rng.choice(AMPS)), kind=str(rng.choice(KINDS))))
    inject(fr, srcs, t_ref)
    Z = maps_for("pre_inj")
    for s in srcs:
        for k in passes:
            dt = (tpass[int(k)] - t_ref) / 365.25
            # true grid position at pass k (grid is north-up, east-left: +RA -> -x)
            gx = s["gx"] - s["pmra"] * dt
            gy = s["gy"] + s["pmdec"] * dt
            row = dict(real=r, kind=s["kind"], amp=s["amp"], pm=s["pm"], pas=int(k), gx=gx, gy=gy)
            for key, Zk in Z.items():
                row["|".join(key)] = local_max(Zk[int(k)], gx, gy) if int(k) in Zk else np.nan
            rows.append(row)
    print(f"real {r+1}/{nreal} {time.time()-t0:.0f}s", flush=True)

np.savez(os.path.join(out, f"null_{seed}.npz"), **null)
import pandas as pd
df = pd.DataFrame(rows)
df.to_csv(os.path.join(out, f"inj_{seed}.csv"), index=False)
# completeness at false-alarm probability 1% per (position, pass)
summ = {}
for key in null:
    thr = np.nanpercentile(null[key], 99)
    df["det_" + key] = df[key] > thr
    summ[key] = float(thr)
json.dump(summ, open(os.path.join(out, f"thr_{seed}.json"), "w"), indent=1)
for g in GROUPS:
    print("== group", g, "thresholds", {m: round(summ[f"{m}|{g}"], 2) for m in METHODS})
    tab = df.groupby(["kind", "amp"])[[f"det_{m}|{g}" for m in METHODS]].mean().round(2)
    tab.columns = METHODS
    print(tab)
print("by pm (all)")
tab = df.groupby(["pm"])[[f"det_{m}|all" for m in METHODS]].mean().round(2)
tab.columns = METHODS
print(tab)
