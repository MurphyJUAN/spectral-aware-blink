# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""E5: survey-geometry rules over many sky positions (header-only, no pixels).
For random positions inside the SIA footprint: within-pass time<->lambda aliasing, same-colour pairing
usability f(tau) vs a Poisson prediction, B2 effective-wavelength mismatch, predicted SNR gain C vs B1.
Usage: python e5_coverage.py <coverage_npz> <ra0> <dec0> <radius_deg> <n_pos> <outdir>
"""
import json
import os
import sys
import warnings
import logging

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning
from scipy.ndimage import map_coordinates
from scipy.stats import spearmanr

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", AstropyWarning)
logging.getLogger("astropy").setLevel(logging.ERROR)

cov, ra0, dec0, rad, npos, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
os.makedirs(out, exist_ok=True)
z = np.load(cov)
H, mjd, det = z["hdr"], z["mjd"], z["det"]
step = int(z["step"])
order = np.argsort(mjd)
pas = np.zeros(len(mjd), int)
p = 0
for a, b in zip(order[:-1], order[1:]):
    if mjd[b] - mjd[a] > 20:
        p += 1
    pas[b] = p
rng = np.random.default_rng(1)
r = rad * np.sqrt(rng.uniform(0, 1, npos))
th = rng.uniform(0, 2 * np.pi, npos)
pdec = dec0 + r * np.sin(th)
pra = ra0 + r * np.cos(th) / np.cos(np.radians(dec0))
LAM = np.full((len(H), npos), np.nan)
for i, h in enumerate(H):
    w = WCS(fits.Header.fromstring(str(h)))
    x, y = w.all_world2pix(pra, pdec, 0)
    on = (x > 4) & (x < 2035) & (y > 4) & (y < 2035)
    lm = z[f"lam{det[i]}"]
    lam = map_coordinates(lm.astype(float), [y / step, x / step], order=1, mode="nearest")
    LAM[i, on] = lam[on]
TAUS = [0.0025, 0.005, 0.01, 0.02, 0.04]
rows = []
for j in range(npos):
    for d in range(1, 7):
        m = (det == d) & np.isfinite(LAM[:, j])
        if m.sum() < 8:
            continue
        l, t, pp = LAM[m, j], mjd[m], pas[m]
        lo, hi = 0.75 * 0, None
        lnr = np.log(z[f"lam{d}"].max() / z[f"lam{d}"].min())
        rho = []
        for k in np.unique(pp):
            mk = pp == k
            if mk.sum() >= 5:
                rho.append((spearmanr(t[mk], l[mk]).statistic, int(mk.sum()), float(np.ptp(np.log(l[mk])) / lnr)))
        row = dict(pos=j, ra=float(pra[j]), dec=float(pdec[j]), det=d, n=int(m.sum()), npass=int(len(np.unique(pp))),
                   rho=[x[0] for x in rho], rho_n=[x[1] for x in rho], pass_lnspan_frac=[x[2] for x in rho])
        # same-colour pairing usability and predicted gain
        dl = np.abs(np.log(l[:, None] / l[None, :]))
        other = pp[:, None] != pp[None, :]
        for tau in TAUS:
            mm = (dl < tau) & other
            npart = mm.sum(1)
            f = float(np.mean(npart > 0))
            nother = other.sum(1)
            fpred = float(np.mean(1 - np.exp(-nother * 2 * tau / lnr)))
            # noise-limited gain: C uses every frame with ref noise ~ sigma^2/n_other; B1 uses frames with partners, ref noise sigma^2/npart
            infoC = np.sum(1 / (1 + 1 / np.maximum(nother, 1)))
            infoB1 = np.sum(np.where(npart > 0, 1 / (1 + 1 / np.maximum(npart, 1)), 0))
            row[f"f_{tau}"] = f
            row[f"fpred_{tau}"] = fpred
            row[f"gain_{tau}"] = float(np.sqrt(infoC / infoB1)) if infoB1 > 0 else np.inf
        # B2 effective wavelength mismatch between pass stacks (ln-lambda mean per pass)
        mu = [np.mean(np.log(l[pp == k])) for k in np.unique(pp)]
        row["b2_dlnlam_max"] = float(np.ptp(mu)) if len(mu) > 1 else np.nan
        rows.append(row)
json.dump(rows, open(os.path.join(out, "e5_rows.json"), "w"))
import pandas as pd
df = pd.DataFrame(rows)
df["rho_med"] = df["rho"].apply(lambda v: np.median(v) if len(v) else np.nan)
df["absrho_wmean"] = df.apply(lambda r: np.average(np.abs(r["rho"]), weights=r["rho_n"]) if len(r["rho"]) else np.nan, axis=1)
df["span_med"] = df["pass_lnspan_frac"].apply(lambda v: np.median(v) if len(v) else np.nan)
cols = ["n", "npass", "absrho_wmean", "span_med", "b2_dlnlam_max"] + [f"f_{t}" for t in TAUS] + [f"fpred_{t}" for t in TAUS] + [f"gain_{t}" for t in TAUS]
pd.set_option("display.width", 250)
print(df.groupby("det")[cols].median().round(3).T)
print("positions", df.pos.nunique())
