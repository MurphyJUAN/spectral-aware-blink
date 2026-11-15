# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Fetch SPHEREx L2 cutouts from the S3 mirror via HTTP range reads.

Usage: python fetch.py <sia_csv> <ra> <dec> <size_px> <outdir> [nthreads]
Writes one npz per image: img, var, zodi, flags, lam, bw, x0, y0, mjd, det, obsid, hdr
"""
import csv
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", AstropyWarning)
import logging
logging.getLogger("astropy").setLevel(logging.ERROR)

IRSA = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"
S3 = "https://nasa-irsa-spherex.s3.amazonaws.com/"


def to_s3(url):
    return url.replace(IRSA, S3)


def fetch_one(row, ra, dec, size, outdir):
    url = to_s3(row["access_url"])
    obsid = row["obs_id"]
    det = int(row["energy_bandpassname"][-1])
    out = os.path.join(outdir, f"{obsid}_D{det}.npz")
    if os.path.exists(out):
        return out, "cached"
    for attempt in range(4):
        try:
            with fits.open(url, use_fsspec=True, fsspec_kwargs={"block_size": 1024 * 1024, "client_kwargs": {"trust_env": True}}) as h:
                hd = h["IMAGE"].header
                w = WCS(hd)
                x, y = w.all_world2pix(ra, dec, 0)
                xc, yc = int(round(float(x))), int(round(float(y)))
                x0, y0 = max(0, xc - size // 2), max(0, yc - size // 2)
                x1, y1 = min(2040, xc + size // 2), min(2040, yc + size // 2)
                if x1 - x0 < 8 or y1 - y0 < 8:
                    return out, "offchip"
                img = np.asarray(h["IMAGE"].section[y0:y1, x0:x1], dtype=np.float32)
                var = np.asarray(h["VARIANCE"].section[y0:y1, x0:x1], dtype=np.float32)
                zodi = np.asarray(h["ZODI"].section[y0:y1, x0:x1], dtype=np.float32)
                flags = np.asarray(h["FLAGS"].section[y0:y1, x0:x1], dtype=np.int32)
                wave = WCS(hd, fobj=h, key="W")
                yy, xx = np.mgrid[y0:y1, x0:x1]
                lam, bw = wave.pixel_to_world_values(xx, yy)
                np.savez_compressed(
                    out, img=img, var=var, zodi=zodi, flags=flags,
                    lam=np.asarray(lam, np.float32), bw=np.asarray(bw, np.float32),
                    x0=x0, y0=y0, xc=float(x), yc=float(y),
                    mjd=float(hd.get("MJD-AVG", float(row["t_min"]))), det=det, obsid=obsid,
                    fwhm=float(hd.get("PSF_FWHM", np.nan)),
                    hdr=hd.tostring(),
                )
            return out, "ok"
        except Exception as e:  # retry with backoff
            err = repr(e)
            time.sleep(2 ** attempt)
    return out, "fail " + err


def main():
    sia, ra, dec, size, outdir = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    nth = int(sys.argv[6]) if len(sys.argv) > 6 else 24
    os.makedirs(outdir, exist_ok=True)
    rows = [r for r in csv.DictReader(open(sia)) if r["dataproduct_type"] == "image"]
    t0 = time.time()
    stats = {}
    with ThreadPoolExecutor(nth) as ex:
        futs = [ex.submit(fetch_one, r, ra, dec, size, outdir) for r in rows]
        for i, f in enumerate(as_completed(futs)):
            out, st = f.result()
            k = st.split()[0]
            stats[k] = stats.get(k, 0) + 1
            if k == "fail":
                print(out, st, flush=True)
            if i % 20 == 0:
                print(f"{i+1}/{len(rows)} {stats} {time.time()-t0:.0f}s", flush=True)
    print("DONE", stats, f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
