# SI-generated code: written with Claude Code (Anthropic) from the team's idea and under the team's direction, NASA Space Apps 2026.
"""Core: load cutouts, preprocess, reproject to a common grid, per-pixel spectral models (LOPO)."""
import glob
import logging
import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning
from scipy.ndimage import map_coordinates

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", AstropyWarning)
logging.getLogger("astropy").setLevel(logging.ERROR)

BAD_BITS = [0, 1, 2, 6, 7, 9, 10, 11, 15, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28]
BAD_MASK = sum(1 << b for b in BAD_BITS)
PIX = 6.15 / 3600.0  # deg


def load_cache(d):
    fr = []
    for f in sorted(glob.glob(os.path.join(d, "*.npz"))):
        z = np.load(f, allow_pickle=False)
        e = {k: z[k] for k in z.files}
        e["hdr"] = str(e["hdr"])
        e["obsid"] = str(e["obsid"])
        e["det"] = int(e["det"])
        e["mjd"] = float(e["mjd"])
        e["file"] = f
        fr.append(e)
    fr.sort(key=lambda e: (e["mjd"], e["det"]))
    # passes: split on gaps > 20 days
    t = np.array([e["mjd"] for e in fr])
    order = np.argsort(t)
    p, last = 0, None
    for i in order:
        if last is not None and t[i] - last > 20:
            p += 1
        fr[i]["pass"] = p
        last = t[i]
    return fr


def preprocess(e):
    img = e["img"].astype(np.float64) - e["zodi"]
    var = e["var"].astype(np.float64).copy()
    bad = ((e["flags"] & BAD_MASK) != 0) | ~np.isfinite(img) | ~np.isfinite(var) | (var <= 0)
    img[bad] = np.nan
    img = img - mesh_background(img, BOX)
    if SMOOTH > 0:  # pre-smooth undersampled native frame (normalized convolution) before resampling
        from scipy.ndimage import convolve
        r = int(np.ceil(3 * SMOOTH))
        g = np.exp(-0.5 * (np.arange(-r, r + 1) / SMOOTH) ** 2)
        k = np.outer(g, g); k /= k.sum()
        ok = np.isfinite(img).astype(float)
        num = convolve(np.where(ok > 0, img, 0.0), k, mode="constant")
        den = convolve(ok, k, mode="constant")
        vv = convolve(np.where(ok > 0, var, 0.0), k ** 2, mode="constant")
        good = den > 0.7
        img = np.where(good, num / np.where(good, den, 1), np.nan)
        var = np.where(good, vv / np.where(good, den, 1) ** 2, np.inf)
        bad = ~good
    var[bad] = np.inf
    return img, var, bad


SMOOTH = float(os.environ.get("SMOOTH", "0"))


BOX = int(os.environ.get("BKG_BOX", "16"))


def _clipmed(v):
    v = v[np.isfinite(v)]
    if v.size < 20:
        return np.nan
    m = np.median(v)
    for _ in range(5):
        s = 1.4826 * np.median(np.abs(v - m)) + 1e-12
        v = v[np.abs(v - m) < 2.5 * s]
        if v.size < 10:
            break
        m = np.median(v)
    return m


def mesh_background(img, box=16):
    """SExtractor-like mesh background: clipped median per box, 3x3 median filtered, bilinear to full res."""
    from scipy.ndimage import median_filter, zoom
    ny, nx = img.shape
    if box <= 0:
        return np.full(img.shape, _clipmed(img.ravel()))
    my, mx = max(1, ny // box), max(1, nx // box)
    ey, ex = np.linspace(0, ny, my + 1).astype(int), np.linspace(0, nx, mx + 1).astype(int)
    M = np.array([[_clipmed(img[ey[i]:ey[i + 1], ex[j]:ex[j + 1]].ravel()) for j in range(mx)] for i in range(my)])
    if np.all(~np.isfinite(M)):
        return np.zeros_like(img)
    M = np.where(np.isfinite(M), M, np.nanmedian(M))
    if min(M.shape) >= 3:
        M = median_filter(M, 3, mode="nearest")
    # bilinear interpolation from box centres
    cy = (ey[:-1] + ey[1:]) / 2 - 0.5
    cx = (ex[:-1] + ex[1:]) / 2 - 0.5
    from scipy.interpolate import RegularGridInterpolator
    if len(cy) == 1 or len(cx) == 1:
        return np.full(img.shape, float(np.mean(M)))
    f = RegularGridInterpolator((cy, cx), M, bounds_error=False, fill_value=None)
    yy, xx = np.mgrid[0:ny, 0:nx]
    return f((np.clip(yy, cy[0], cy[-1]), np.clip(xx, cx[0], cx[-1])))


def tan_grid(ra0, dec0, n, scale_deg=PIX):
    """Common tangent-plane grid (north up, east left)."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    w.wcs.cdelt = [-scale_deg, scale_deg]
    yy, xx = np.mgrid[0:n, 0:n]
    ra, dec = w.wcs_pix2world(xx, yy, 0)
    return w, ra, dec


def reproject_frames(frames, ra, dec):
    """Bilinear-interpolate each preprocessed frame onto (ra, dec) grid.
    Returns Y, V, L (n_img, ny, nx), and native pixel coords of the grid in each frame."""
    n = len(frames)
    sh = ra.shape
    Y = np.full((n,) + sh, np.nan)
    V = np.full((n,) + sh, np.inf)
    L = np.full((n,) + sh, np.nan)
    PH = np.zeros((n,) + sh + (2,))
    for i, e in enumerate(frames):
        w = WCS(fits.Header.fromstring(e["hdr"]))
        x, y = w.all_world2pix(ra, dec, 0)
        x = x - e["x0"] + e.get("dx", 0.0)
        y = y - e["y0"] + e.get("dy", 0.0)
        img, var, bad = e["pre"]
        ok = np.isfinite(img)
        im0 = np.where(ok, img, 0.0)
        cov = map_coordinates(ok.astype(float), [y, x], order=1, mode="constant", cval=0.0)
        val = map_coordinates(im0, [y, x], order=1, mode="constant", cval=0.0)
        vv = map_coordinates(np.where(ok, var, 0.0), [y, x], order=1, mode="constant", cval=0.0)
        lam = map_coordinates(e["lam"].astype(float), [y, x], order=1, mode="nearest")
        good = cov > 0.999
        Y[i][good] = val[good]
        V[i][good] = vv[good]
        L[i] = lam
        e["gx"], e["gy"] = x, y
        PH[i, ..., 0] = x - np.floor(x)
        PH[i, ..., 1] = y - np.floor(y)
    reproject_frames.last_phase = PH
    return Y, V, L


# ---------------- spectral models -----------------

def bspline_basis(x, knots_inner, lo, hi, deg=3):
    """Evaluate cubic B-spline basis at x (any shape). Returns (..., K)."""
    from scipy.interpolate import BSpline
    t = np.r_[[lo] * (deg + 1), knots_inner, [hi] * (deg + 1)]
    K = len(t) - deg - 1
    xs = np.clip(x, lo, hi - 1e-9)
    out = np.empty(x.shape + (K,))
    flat = xs.ravel()
    for k in range(K):
        c = np.zeros(K)
        c[k] = 1
        out.reshape(-1, K)[:, k] = BSpline(t, c, deg, extrapolate=False)(flat)
    out = np.nan_to_num(out)
    return out


def fit_predict(Y, V, L, train, test, model="spline", nk=8, alpha=1.0, eta=0.0, lo=None, hi=None, iters=2, clip=4.0, X=None, ridge_x=1.0):
    """Per-pixel weighted fit on 'train' images, predict 'test' images.
    Y,V,L: (n, P) arrays. Returns pred (n_test, P), predvar (n_test, P)."""
    Yt, Vt, Lt = Y[train], V[train], L[train]
    ok = np.isfinite(Yt) & np.isfinite(Vt)
    if model == "const":
        B = np.ones(Yt.shape + (1,))
        Bs = np.ones(L[test].shape + (1,))
        D = np.zeros((1, 1))
    else:
        lo = np.nanmin(L) if lo is None else lo
        hi = np.nanmax(L) if hi is None else hi
        lo, hi = lo - 1e-4, hi + 1e-4
        if model == "linear":
            B = np.stack([np.ones_like(Lt), (Lt - lo) / (hi - lo)], -1)
            Bs = np.stack([np.ones_like(L[test]), (L[test] - lo) / (hi - lo)], -1)
            D = np.zeros((2, 2))
        else:
            kn = np.linspace(lo, hi, nk - 2)[1:-1]
            B = bspline_basis(np.nan_to_num(Lt, nan=lo), kn, lo, hi)
            Bs = bspline_basis(np.nan_to_num(L[test], nan=lo), kn, lo, hi)
            K = B.shape[-1]
            Dm = np.diff(np.eye(K), 2, axis=0)
            D = alpha * Dm.T @ Dm
    if X is not None:  # extra per-(image,pixel) regressors, e.g. sub-pixel phase terms (n, P, E)
        E_ = X.shape[-1]
        B = np.concatenate([B, X[train]], -1)
        Bs = np.concatenate([Bs, X[test]], -1)
        D2 = np.zeros((D.shape[0] + E_, D.shape[1] + E_))
        D2[: D.shape[0], : D.shape[1]] = D
        D = D2
        Dx = np.zeros(D.shape[0]); Dx[-E_:] = ridge_x
        D = D + np.diag(Dx)
    w = np.where(ok, 1.0 / np.where(ok, Vt, 1.0), 0.0)
    Yz = np.where(ok, Yt, 0.0)
    for it in range(iters):
        A = np.einsum("np,npk,npl->pkl", w, B, B) + D[None]
        b = np.einsum("np,npk,np->pk", w, B, Yz)
        A += 1e-10 * np.eye(A.shape[-1])[None] * (1 + np.nanmax(np.abs(A), axis=(1, 2)))[:, None, None]
        c = np.linalg.solve(A, b[..., None])[..., 0]
        if it < iters - 1:  # clip outliers (cosmic rays, movers) and refit
            m = np.einsum("npk,pk->np", B, c)
            s2 = Vt + (eta * m) ** 2
            r = (Yt - m) / np.sqrt(s2)
            w = np.where(ok & (np.abs(r) < clip), 1.0 / np.where(ok, s2, 1.0), 0.0)
    pred = np.einsum("npk,pk->np", Bs, c)
    Ainv = np.linalg.inv(A)
    pv = np.einsum("npk,pkl,npl->np", Bs, Ainv, Bs)
    return pred, pv


def lopo_residuals(Y, V, L, passes, model="spline", eta=0.0, X=None, **kw):
    """Leave-one-pass-out residuals for one detector. Y,V,L: (n, P)."""
    R = np.full(Y.shape, np.nan)
    S = np.full(Y.shape, np.nan)  # sigma
    P = np.full(Y.shape, np.nan)
    lo, hi = np.nanmin(L), np.nanmax(L)
    for k in np.unique(passes):
        te = passes == k
        tr = ~te
        if tr.sum() < 3:
            continue
        pred, pv = fit_predict(Y, V, L, tr, te, model=model, eta=eta, lo=lo, hi=hi, X=X, **kw)
        s2 = V[te] + pv + (eta * pred) ** 2
        R[te] = Y[te] - pred
        S[te] = np.sqrt(s2)
        P[te] = pred
    return R, S, P


def samecolor_residuals(Y, V, L, passes, tol=0.01):
    """B1: predict each pixel from other-pass images with |dlam/lam|<tol (inverse-variance mean)."""
    n, Pn = Y.shape
    R = np.full(Y.shape, np.nan)
    S = np.full(Y.shape, np.nan)
    for i in range(n):
        m = (passes != passes[i])[:, None] & (np.abs(L - L[i][None]) / L[i][None] < tol) & np.isfinite(Y) & np.isfinite(V)
        w = np.where(m, 1.0 / np.where(m, V, 1.0), 0.0)
        sw = w.sum(0)
        good = sw > 0
        pred = np.where(good, (w * np.nan_to_num(Y)).sum(0) / np.where(good, sw, 1), np.nan)
        R[i] = np.where(good, Y[i] - pred, np.nan)
        S[i] = np.where(good, np.sqrt(V[i] + 1.0 / np.where(good, sw, 1)), np.nan)
    return R, S


def phase_features(PH, scale=None):
    """Sub-pixel phase regressors (n, P, 4) from PH (n, P, 2); scaled per pixel by 'scale' (P,) if given."""
    a = 2 * np.pi * PH
    F = np.stack([np.cos(a[..., 0]), np.sin(a[..., 0]), np.cos(a[..., 1]), np.sin(a[..., 1])], -1)
    if scale is not None:
        F = F * np.nan_to_num(scale, nan=0.0)[None, :, None]
    return np.nan_to_num(F)
