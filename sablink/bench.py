# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""Benchmark harness: methods (B0 naive, B1 same-colour, B2 detector-stack, C spectral LOPO),
injection of moving point sources, per-pass residual SNR maps, detection metrics."""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import convolve

import os as _os
from core import preprocess, reproject_frames, lopo_residuals, samecolor_residuals, tan_grid, phase_features

METHODS = ["B0_naive", "B1_samecolor", "B2_detstack", "C_spectral", "C2_spec_phase"]


# ---------------- SED templates (relative, as function of lambda in micron) ---------------
def sed(kind, lam):
    lam = np.asarray(lam, float)
    if kind == "star":  # Rayleigh-Jeans-ish tail in surface brightness per unit freq ~ lam^-2 ... use F_nu ~ lam^-1
        return (lam / 2.0) ** -1.0
    if kind == "tdwarf":  # T dwarf: peaks J/H/K, deep CH4 + H2O at 1.15, 1.4, 1.6-1.8, 2.2+, 3.3; bright 4-5um
        f = 1.0 + 0 * lam
        for c, wdt, dep in [(1.15, 0.05, 0.6), (1.4, 0.08, 0.85), (1.75, 0.12, 0.8), (2.3, 0.2, 0.7), (3.3, 0.2, 0.9), (2.7, 0.3, 0.8)]:
            f *= 1 - dep * np.exp(-0.5 * ((lam - c) / wdt) ** 2)
        f *= 1 + 1.5 * np.exp(-0.5 * ((lam - 4.4) / 0.4) ** 2)
        return f
    if kind == "ydwarf":  # Y dwarf: essentially only bright at 4-5 um
        return 0.03 + np.exp(-0.5 * ((lam - 4.5) / 0.35) ** 2)
    raise ValueError(kind)


def inject(frames, srcs, t_ref):
    """srcs: list of dict(ra, dec, pmra_pxyr, pmdec_pxyr, amp, kind). amp = peak SB at sed==1 [MJy/sr].
    Adds Gaussian PSF to preprocessed images in native pixels."""
    for e in frames:
        img, var, bad = e["pre"]
        img = img.copy()
        w = WCS(fits.Header.fromstring(e["hdr"]))
        fwhm_px = (e["fwhm"] if np.isfinite(e["fwhm"]) else 7.0) / 6.15
        sig = fwhm_px / 2.3548
        ny, nx = img.shape
        gy, gx = np.mgrid[0:ny, 0:nx]
        dt = (e["mjd"] - t_ref) / 365.25
        for s in srcs:
            dec = s["dec"] + s["pmdec"] * dt * 6.15 / 3600
            ra = s["ra"] + s["pmra"] * dt * 6.15 / 3600 / np.cos(np.radians(s["dec"]))
            x, y = w.all_world2pix(ra, dec, 0)
            x -= e["x0"]
            y -= e["y0"]
            if not (-5 < x < nx + 5 and -5 < y < ny + 5):
                continue
            lam = e["lam"][int(np.clip(round(float(y)), 0, ny - 1)), int(np.clip(round(float(x)), 0, nx - 1))]
            a = s["amp"] * sed(s["kind"], lam)
            img += a * np.exp(-0.5 * ((gx - x) ** 2 + (gy - y) ** 2) / sig ** 2)
        e["pre_inj"] = (img, var, bad)


KMATCH = int(_os.environ.get("KMATCH", "0")) if "_os" in globals() else 0


def kernel_match(R, S, Y, shape, half=2):
    """Per-image PSF/registration matching (Alard-Lupton, delta basis): fit r_i ~ pred_i (*) K_i (K is (2h+1)^2)
    by weighted least squares over the whole field, subtract. pred_i = Y_i - R_i. Applied identically to all methods."""
    ny, nx = shape
    out = R.copy()
    for i in range(R.shape[0]):
        r, s = R[i], S[i]
        pred = Y[i] - r
        ok = np.isfinite(r) & np.isfinite(s) & np.isfinite(pred) & (s > 0)
        if ok.sum() < 100:
            continue
        P2 = np.where(ok, pred, 0.0).reshape(ny, nx)
        cols = []
        for dy in range(-half, half + 1):
            for dx in range(-half, half + 1):
                cols.append(np.roll(np.roll(P2, dy, 0), dx, 1).ravel())
        A = np.stack(cols, 1)
        edge = np.zeros((ny, nx), bool)
        edge[half:ny - half, half:nx - half] = True
        use = ok & edge.ravel() & (np.abs(np.where(ok, pred, 0)) > 5 * np.where(ok, s, np.inf))  # high-SNR template pixels only (avoid errors-in-variables bias)
        if use.sum() < 3 * A.shape[1]:
            continue
        w = 1 / s[use]
        # robust: two passes, clip |z|>5 (movers, cosmic rays)
        sel = np.ones(use.sum(), bool)
        for _ in range(2):
            coef, *_ = np.linalg.lstsq(A[use][sel] * w[sel, None], r[use][sel] * w[sel], rcond=None)
            z = (r[use] - A[use] @ coef) * w
            sel = np.abs(z) < 5
        out[i] = np.where(ok, r - A @ coef, np.nan)
    return out


def run_methods(Y, V, L, det, pas, mjd, eta=None, methods=METHODS, nk=8, alpha=None, tol=0.01, PH=None, shape=None):
    alpha = float(_os.environ.get("PALPHA", "1.0")) if alpha is None else alpha
    """Y,V,L (n, P). Returns {method: (R, S)} residual and sigma arrays (n, P)."""
    eta = eta or {}
    out = {m: (np.full(Y.shape, np.nan), np.full(Y.shape, np.nan)) for m in methods}
    for d in np.unique(det):
        m = det == d
        idx = np.nonzero(m)[0]
        y, v, l, p, t = Y[m], V[m], L[m], pas[m], mjd[m]
        for meth in methods:
            if meth == "C_spectral":
                R, S, _ = lopo_residuals(y, v, l, p, model="spline", eta=eta.get(meth, 0.0), nk=nk, alpha=alpha)
            elif meth == "C2_spec_phase":
                # phase terms scaled by the pixel's typical brightness so one ridge works for all pixels
                scale = np.maximum(np.nanmedian(np.abs(y), 0), np.sqrt(np.nanmedian(v, 0)))
                Xp = phase_features(PH[m], scale)
                R, S, _ = lopo_residuals(y, v, l, p, model="spline", eta=eta.get(meth, 0.0), nk=nk, alpha=alpha, X=Xp, ridge_x=float(np.nanmedian(1 / v)) * 2.0 / 1.0)
            elif meth == "B2_detstack":
                R, S, _ = lopo_residuals(y, v, l, p, model="const", eta=eta.get(meth, 0.0))
            elif meth == "B1_samecolor":
                R, S = samecolor_residuals(y, v, l, p, tol=tol)
                e_ = eta.get(meth, 0.0)
                if e_:
                    S = np.sqrt(S ** 2 + (e_ * (y - R)) ** 2)
            elif meth == "B0_naive":
                R = np.full(y.shape, np.nan)
                S = np.full(y.shape, np.nan)
                for i in range(len(y)):
                    o = np.nonzero(p != p[i])[0]
                    if len(o) == 0:
                        continue
                    j = o[np.argmin(np.abs(t[o] - t[i]))]
                    R[i] = y[i] - y[j]
                    S[i] = np.sqrt(v[i] + v[j] + (eta.get(meth, 0.0) * y[j]) ** 2)
            R[~np.isfinite(S)] = np.nan
            if _os.environ.get("KMATCH", "0") == "1" and shape is not None:
                R = kernel_match(R, S, y, shape, half=int(_os.environ.get("KHALF", "2")))
            out[meth][0][idx] = R
            out[meth][1][idx] = S
    return out


import os as _os
_PSF_SIG = float(_os.environ.get("PSF_SIG", str(float(np.hypot(0.55, float(_os.environ.get("SMOOTH", "0")))))))


def pass_snr_maps(R, S, pas, shape, psf_sigma_px=_PSF_SIG, dets=None, det=None):
    """Per-pass matched-filter SNR maps from per-image residuals (all detectors or subset).
    Returns dict pass -> Z map (ny, nx), and coverage (# images contributing)."""
    ny, nx = shape
    rr = np.arange(-3, 4)
    ker = np.exp(-0.5 * (rr[:, None] ** 2 + rr[None] ** 2) / max(psf_sigma_px, 0.5) ** 2)
    ker /= ker.sum()
    Z = {}
    for k in np.unique(pas):
        m = pas == k
        if dets is not None:
            m &= np.isin(det, dets)
        r, s = R[m], S[m]
        ok = np.isfinite(r) & np.isfinite(s) & (s > 0)
        w = np.where(ok, 1 / np.where(ok, s, 1) ** 2, 0)
        num = (w * np.nan_to_num(r)).sum(0).reshape(ny, nx)
        den = w.sum(0).reshape(ny, nx)
        # matched filter: sum(k*num) / sqrt(sum(k^2*den))
        fn = convolve(num, ker, mode="constant")
        fd = convolve(den, ker ** 2, mode="constant")
        z = np.where(fd > 0, fn / np.sqrt(np.where(fd > 0, fd, 1)), np.nan)
        Z[int(k)] = z
    return Z


def chi2_stats(R, S):
    z = R / S
    return np.nanmean(z ** 2, axis=0), np.isfinite(z).mean(axis=0)
