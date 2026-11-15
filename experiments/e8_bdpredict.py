# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""E8: predicted false variability of a lambda-blind detector stack (after per-image zero-point matching to the
field-star population) for different SED classes, over many sky positions, using real SPHEREx wavelength sampling.
FV_T(k) = ln(<s_T(lam)>_k / <s_T(lam)>_all) - ln(<s_field(lam)>_k / <s_field(lam)>_all)
Usage: python e8_bdpredict.py <coverage_npz> <ra0> <dec0> <radius_deg> <n_pos> <outdir>
"""
import json
import os
import sys
import warnings
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning
from scipy.ndimage import map_coordinates

from bench import sed

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
rng = np.random.default_rng(2)
r = rad * np.sqrt(rng.uniform(0, 1, npos))
th = rng.uniform(0, 2 * np.pi, npos)
pdec = dec0 + r * np.sin(th)
pra = ra0 + r * np.cos(th) / np.cos(np.radians(dec0))
LAM = np.full((len(H), npos), np.nan)
for i, h in enumerate(H):
    w = WCS(fits.Header.fromstring(str(h)))
    x, y = w.all_world2pix(pra, pdec, 0)
    on = (x > 4) & (x < 2035) & (y > 4) & (y < 2035)
    lam = map_coordinates(np.nan_to_num(z[f"lam{det[i]}"].astype(float), nan=np.nanmedian(z[f"lam{det[i]}"])), [y / step, x / step], order=1, mode="nearest")
    LAM[i, on] = lam[on]
# field-star slope per detector, measured on real stars (E6, median beta of ln F vs ln lam)
BETA_FIELD = {1: 0.54, 2: 0.19, 3: -1.25, 4: -1.63, 5: -1.89, 6: -1.32}
TEMPL = {
    "field star (beta_field)": None,
    "red star (beta_field+1)": 1.0,
    "blue star (beta_field-1)": -1.0,
    "T dwarf (schematic)": "tdwarf",
    "Y dwarf (schematic)": "ydwarf",
}
rows = []
for j in range(npos):
    for d in range(1, 7):
        m = (det == d) & np.isfinite(LAM[:, j])
        if m.sum() < 8:
            continue
        l, pp = LAM[m, j], pas[m]
        sf = (l / 2.0) ** BETA_FIELD[d]
        for k in np.unique(pp):
            mk = pp == k
            if mk.sum() < 3 or (~mk).sum() < 3:
                continue
            fvf = np.log(sf[mk].mean() / sf.mean())
            row = dict(pos=j, det=d, pas=int(k), dmu=float(np.log(l[mk]).mean() - np.log(l).mean()), fv_field_raw=float(fvf))
            for name, t in TEMPL.items():
                if t is None:
                    s = sf
                elif isinstance(t, float):
                    s = (l / 2.0) ** (BETA_FIELD[d] + t)
                else:
                    s = sed(t, l)
                row[name] = float(np.log(s[mk].mean() / s.mean()) - fvf)
            rows.append(row)
import pandas as pd
df = pd.DataFrame(rows)
df.to_csv(os.path.join(out, "e8.csv"), index=False)
pd.set_option("display.width", 250)
print("positions", df.pos.nunique(), "rows", len(df))
print("raw lambda-blind stack false variability of a FIELD star, median |FV| per detector:")
print(df.groupby("det").fv_field_raw.apply(lambda v: np.median(np.abs(v))).round(3).to_dict())
print("after per-image zero-point matching: median |FV| (and 90th pct) by SED class and detector")
for name in TEMPL:
    med = df.groupby("det")[name].apply(lambda v: np.median(np.abs(v))).round(3).to_dict()
    p90 = df.groupby("det")[name].apply(lambda v: np.percentile(np.abs(v), 90)).round(3).to_dict()
    print(f"  {name:26s} median {med}\n  {'':26s} p90    {p90}")
