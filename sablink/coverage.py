# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""Header-only coverage: IMAGE header of every frame + one full-frame lambda map per detector.
Usage: python coverage.py <sia_csv> <out_npz> [nthreads]
"""
import csv
import sys
import time
import warnings
import logging
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", AstropyWarning)
logging.getLogger("astropy").setLevel(logging.ERROR)
from fetch import to_s3

FK = {"block_size": 64 * 1024, "client_kwargs": {"trust_env": True}}


def hdr(row):
    for a in range(4):
        try:
            with fits.open(to_s3(row["access_url"]), use_fsspec=True, fsspec_kwargs=FK) as h:
                hd = h["IMAGE"].header
                return hd.tostring(), float(hd.get("MJD-AVG", float(row["t_min"]))), int(row["energy_bandpassname"][-1]), row["obs_id"]
        except Exception as e:
            time.sleep(2 ** a)
    return None


def lam_map(row, step=8):
    with fits.open(to_s3(row["access_url"]), use_fsspec=True, fsspec_kwargs=FK) as h:
        wave = WCS(h["IMAGE"].header, fobj=h, key="W")
        g = np.arange(0, 2040, step)
        yy, xx = np.meshgrid(g, g, indexing="ij")
        lam, bw = wave.pixel_to_world_values(xx, yy)
    return np.asarray(lam, np.float32), np.asarray(bw, np.float32)


if __name__ == "__main__":
    sia, out = sys.argv[1], sys.argv[2]
    nth = int(sys.argv[3]) if len(sys.argv) > 3 else 32
    rows = [r for r in csv.DictReader(open(sia)) if r["dataproduct_type"] == "image"]
    t0 = time.time()
    with ThreadPoolExecutor(nth) as ex:
        H = list(ex.map(hdr, rows))
    print("headers", sum(h is not None for h in H), time.time() - t0, flush=True)
    lm = {}
    for d in range(1, 7):
        r = next(r for r in rows if r["energy_bandpassname"].endswith(str(d)))
        lm[d] = lam_map(r)
        print("lam map", d, lm[d][0].min(), lm[d][0].max(), time.time() - t0, flush=True)
    ok = [h for h in H if h is not None]
    np.savez_compressed(out, hdr=np.array([h[0] for h in ok]), mjd=np.array([h[1] for h in ok]), det=np.array([h[2] for h in ok]),
                        obsid=np.array([h[3] for h in ok]), **{f"lam{d}": lm[d][0] for d in lm}, **{f"bw{d}": lm[d][1] for d in lm}, step=8)
    print("DONE", time.time() - t0)
