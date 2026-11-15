# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""Source detection on a common-grid stack + forced aperture photometry in native pixels."""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import maximum_filter, map_coordinates

from core import preprocess, tan_grid, reproject_frames


def detect_sources(frames, ra0, dec0, n, nsig=8.0, dets=(1, 2, 3, 4, 5, 6), edge=3):
    w, ra, dec = tan_grid(ra0, dec0, n)
    for e in frames:
        if "pre" not in e:
            e["pre"] = preprocess(e)
    Y, V, L = reproject_frames(frames, ra, dec)
    det = np.array([e["det"] for e in frames])
    stack = []
    for d in dets:
        s = Y[det == d]
        if len(s) == 0:
            continue
        med = np.nanmedian(s, 0)
        mad = 1.4826 * np.nanmedian(np.abs(med - np.nanmedian(med)))
        stack.append(med / mad)
    snr = np.nanmax(np.array(stack), 0)
    snr = np.nan_to_num(snr)
    pk = (snr == maximum_filter(snr, 3)) & (snr > nsig)
    pk[:edge] = pk[-edge:] = False
    pk[:, :edge] = pk[:, -edge:] = False
    yy, xx = np.nonzero(pk)
    sra, sdec = w.wcs_pix2world(xx, yy, 0)
    return dict(ra=sra, dec=sdec, x=xx, y=yy, snr=snr[yy, xx]), (w, ra, dec, Y, V, L, snr)


def forced_phot(frames, sra, sdec, r_ap=2.0, r_in=4.0, r_out=6.0):
    """Returns F, E, LAM, X, Y arrays (n_img, n_src) in MJy/sr * pix units."""
    n, m = len(frames), len(sra)
    F = np.full((n, m), np.nan)
    E = np.full((n, m), np.nan)
    LAM = np.full((n, m), np.nan)
    XX = np.full((n, m), np.nan)
    YY = np.full((n, m), np.nan)
    for i, e in enumerate(frames):
        if "pre" not in e:
            e["pre"] = preprocess(e)
        img, var, bad = e["pre"]
        w = WCS(fits.Header.fromstring(e["hdr"]))
        x, y = w.all_world2pix(sra, sdec, 0)
        x = x - e["x0"]
        y = y - e["y0"]
        ny, nx = img.shape
        gy, gx = np.mgrid[0:ny, 0:nx]
        for j in range(m):
            if not (r_out < x[j] < nx - r_out - 1 and r_out < y[j] < ny - r_out - 1):
                continue
            rr = np.hypot(gx - x[j], gy - y[j])
            ap = rr < r_ap
            an = (rr >= r_in) & (rr < r_out) & np.isfinite(img)
            if np.any(~np.isfinite(img[ap])) or an.sum() < 10:
                continue
            bkg = np.median(img[an])
            F[i, j] = np.sum(img[ap] - bkg)
            E[i, j] = np.sqrt(np.sum(var[ap]))
            XX[i, j], YY[i, j] = x[j], y[j]
        LAM[i] = map_coordinates(e["lam"].astype(float), [np.clip(y, 0, ny - 1), np.clip(x, 0, nx - 1)], order=1)
    return F, E, LAM, XX, YY


def register(frames, sra, sdec, r=2.5, sig=0.7, iters=3, min_stars=3):
    """Per-frame shift (native px) = median over bright stars of (Gaussian-weighted centroid - WCS position).
    Stores e['dx'], e['dy']; returns array of (dx, dy, n)."""
    out = []
    for e in frames:
        img, var, bad = e["pre"]
        w = WCS(fits.Header.fromstring(e["hdr"]))
        x, y = w.all_world2pix(sra, sdec, 0)
        x = x - e["x0"]; y = y - e["y0"]
        ny, nx = img.shape
        gy, gx = np.mgrid[0:ny, 0:nx]
        dx, dy = [], []
        for j in range(len(x)):
            if not (5 < x[j] < nx - 6 and 5 < y[j] < ny - 6):
                continue
            cx, cy = x[j], y[j]
            ok = True
            for _ in range(iters):
                m = np.hypot(gx - cx, gy - cy) < r
                v = img[m]
                if np.any(~np.isfinite(v)):
                    ok = False; break
                wt = np.clip(v, 0, None) * np.exp(-0.5 * ((gx[m] - cx) ** 2 + (gy[m] - cy) ** 2) / sig ** 2)
                if wt.sum() <= 0:
                    ok = False; break
                cx, cy = (wt * gx[m]).sum() / wt.sum(), (wt * gy[m]).sum() / wt.sum()
            if ok and abs(cx - x[j]) < 1 and abs(cy - y[j]) < 1:
                dx.append(cx - x[j]); dy.append(cy - y[j])
        if len(dx) >= min_stars:
            e["dx"], e["dy"] = float(np.median(dx)), float(np.median(dy))
        else:
            e["dx"], e["dy"] = 0.0, 0.0
        out.append((e["dx"], e["dy"], len(dx)))
    return np.array(out)
