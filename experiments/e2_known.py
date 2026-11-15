# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""E2: known mover (WISE 0855) in per-pass residual SNR maps for each method.
Usage: python e2_known.py <cache_dir> <ra0> <dec0> <grid_n> <outdir>
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
import astropy.units as u

from core import load_cache, preprocess, reproject_frames, tan_grid
from bench import run_methods, pass_snr_maps, METHODS

cache, ra0, dec0, n, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
os.makedirs(out, exist_ok=True)
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
N = len(fr)
PH = reproject_frames.last_phase.reshape(N, -1, 2)
det = np.array([e["det"] for e in fr])
pas = np.array([e["pass"] for e in fr])
mjd = np.array([e["mjd"] for e in fr])
res = run_methods(Y.reshape(N, -1), V.reshape(N, -1), L.reshape(N, -1), det, pas, mjd, eta={m: 0.03 for m in METHODS}, PH=PH, shape=(n, n))
c = SkyCoord(133.81890825 * u.deg, -7.247048 * u.deg, pm_ra_cosdec=-8123.7 * u.mas / u.yr, pm_dec=673.2 * u.mas / u.yr,
             distance=Distance(parallax=439.0 * u.mas), obstime=Time("J2000"))
outd = {"pos": {}, "z": {}}
for k in np.unique(pas):
    tk = float(np.mean(mjd[pas == k]))
    p = c.apply_space_motion(Time(tk, format="mjd"))
    x, y = w.wcs_world2pix(p.ra.deg, p.dec.deg, 0)
    outd["pos"][int(k)] = (float(x), float(y), tk)
stamps = {}
for meth in METHODS:
    R, S = res[meth]
    for grp, ds in [("D56", [5, 6]), ("D456", [4, 5, 6]), ("D123", [1, 2, 3])]:
        Z = pass_snr_maps(R, S, pas, (n, n), dets=ds, det=det)
        zs = []
        for k, (x, y, tk) in outd["pos"].items():
            yy, xx = np.mgrid[0:n, 0:n]
            m = np.hypot(xx - x, yy - y) <= 1.0
            zs.append(float(np.nanmax(Z[k][m])))
            stamps[f"{meth}|{grp}|{k}"] = Z[k][int(round(y)) - 8:int(round(y)) + 9, int(round(x)) - 8:int(round(x)) + 9]
        outd["z"][f"{meth}|{grp}"] = zs
json.dump(outd, open(os.path.join(out, "e2.json"), "w"), indent=1)
np.savez(os.path.join(out, "e2_stamps.npz"), **stamps)
print(json.dumps(outd["pos"]))
for k, v in outd["z"].items():
    print(f"{k:28s}", [round(x, 1) for x in v])
