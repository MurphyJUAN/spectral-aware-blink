# SI-assisted: written with Claude Code (Anthropic) under the direction of the team, NASA Space Apps 2026.
"""Figures for the report. Usage: python figures.py <results_dir> <out_dir>"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sablink"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from core import fit_predict

R, OUT = sys.argv[1], sys.argv[2]
E7DIR = sys.argv[3] if len(sys.argv) > 3 else "e7_w0855"
E6DIR = sys.argv[4] if len(sys.argv) > 4 else "e6_160"
E3DIR = sys.argv[5] if len(sys.argv) > 5 else "e3_160_s0"
os.makedirs(OUT, exist_ok=True)
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
PASS_C = ["#2a78d6", "#eb6834", "#1baf7a", "#52514e"]
PASS_M = ["o", "s", "^", "D"]
METH = ["B0_naive", "B1_samecolor", "B2_detstack", "C_spectral", "C2_spec_phase"]
MLAB = {"B0_naive": "B0 naive blink", "B1_samecolor": "B1 same-colour pairs", "B2_detstack": "B2 detector stack",
        "C_spectral": "C spectral LOPO (ours)", "C2_spec_phase": "C2 spectral+phase (ours)"}
MCOL = dict(zip(METH, CAT[:5]))
DIV = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#e8e7e3", "#f2a27f", "#b8461b"])
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
                     "lines.linewidth": 2, "legend.frameon": False})


def fig_problem():
    p = os.path.join(R, "e0_160", "e0_phot.npz")
    if not os.path.exists(p):
        return
    z = np.load(p)
    F, E, LAM, det, pas, mjd = z["F"], z["E"], z["LAM"], z["det"], z["pas"], z["mjd"]
    d = 3
    m = det == d
    snr = np.nanmedian(F[m] / E[m], 0)
    cnt = np.isfinite(F[m]).sum(0)
    j = int(np.nanargmax(np.where((cnt > 0.8 * m.sum()) & (snr < 300), snr, -1)))
    f, e, l, p_, t = F[m, j], E[m, j], LAM[m, j], pas[m], mjd[m]
    ok = np.isfinite(f)
    f, e, l, p_, t = f[ok], e[ok], l[ok], p_[ok], t[ok]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
    days = t - t.min()
    for k in np.unique(p_):
        mk = p_ == k
        lab = f"pass {k+1}"
        ax[0].errorbar(days[mk], f[mk], e[mk], fmt=PASS_M[k], ms=5, color=PASS_C[k], label=lab, lw=1)
        ax[1].errorbar(l[mk], f[mk], e[mk], fmt=PASS_M[k], ms=5, color=PASS_C[k], label=lab, lw=1)
    xs = np.linspace(l.min(), l.max(), 200)
    pr, _ = fit_predict(f[:, None], (e ** 2)[:, None], l[:, None], np.ones(len(f), bool), np.ones(len(f), bool), nk=6, alpha=0.1)
    o = np.argsort(l)
    ax[1].plot(l[o], pr[o, 0], color=INK, lw=2, label="spline SED")
    ax[0].set(xlabel="days since first exposure", ylabel="aperture flux [MJy/sr·pix]", title="Same star vs TIME: looks variable")
    ax[1].set(xlabel="wavelength at star [µm]", title="Same star vs WAVELENGTH: smooth SED")
    ax[0].legend(loc="best", fontsize=8)
    # within-pass sweep
    for k in np.unique(p_):
        mk = p_ == k
        if mk.sum() < 3:
            continue
        tt = (t[mk] - t[mk].min())
        ax[2].plot(tt, l[mk], PASS_M[k] + "-", color=PASS_C[k], ms=5, lw=1.2, label=f"pass {k+1}")
    ax[2].set(xlabel="days since start of pass", ylabel="wavelength at star [µm]", title="Within a pass, λ sweeps with time")
    ax[2].legend(fontsize=8)
    fig.suptitle(f"E0 — WISE 0855 field, detector D{d}, star SNR≈{snr[j]:.0f}: frame-to-frame 'flicker' is the LVF reading out the spectrum", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_problem.png"), dpi=150)
    plt.close(fig)


def fig_b2law():
    p = os.path.join(R, E6DIR, "e6.json")
    if not os.path.exists(p):
        return
    rows = json.load(open(p))
    obs = np.array([r["obs"] for r in rows])
    pred = np.array([r["pred"] for r in rows])
    dd = np.array([r["det"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.4))
    lim = np.nanpercentile(np.abs(np.r_[obs, pred]), 99.5) * 1.1
    ax[0].plot([-lim, lim], [-lim, lim], color=INK2, lw=1, ls="--")
    for i, d in enumerate(range(1, 7)):
        m = dd == d
        ax[0].scatter(pred[m], obs[m], s=12, color=CAT[i], label=f"D{d}", edgecolor=SURF, linewidth=0.4)
    ax[0].set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="predicted ln(F_pass/F_all) from SED × λ-sampling (LOPO)",
              ylabel="observed ln(F_pass/F_all), detector stack", title="Detector-stack 'variability' is predictable")
    ax[0].legend(fontsize=8, ncol=2)
    bins = np.linspace(-lim, lim, 41)
    ax[1].hist(obs, bins, color=CAT[1], alpha=0.85, label=f"B2 raw: rms {np.sqrt(np.mean(obs**2)):.3f}")
    ax[1].hist(obs - pred, bins, color=CAT[0], alpha=0.85, label=f"after spectral correction: rms {np.sqrt(np.mean((obs-pred)**2)):.3f}")
    ax[1].set(xlabel="ln flux ratio between pass stacks", ylabel="star × pass × detector", title="Stationary stars: false variability removed")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_b2law.png"), dpi=150)
    plt.close(fig)


def fig_residual_maps():
    d = os.path.join(R, "e1_160_small")
    if not os.path.exists(os.path.join(d, "Zp.npz")):
        return
    Zs = dict(np.load(os.path.join(d, "Zp.npz")))
    g = np.load(os.path.join(d, "grid.npz"))
    k = 1
    show = ["B0_naive", "B1_samecolor", "B2_detstack", "C_spectral"]
    fig, ax = plt.subplots(1, len(show) + 1, figsize=(3.1 * (len(show) + 1), 3.5))
    ax[0].imshow(np.arcsinh(g["snr"] / 10), origin="lower", cmap="Greys_r")
    ax[0].set_title("sky (median SNR stack)", fontsize=9)
    for a, m in zip(ax[1:], show):
        Z = Zs[m]
        im = a.imshow(Z[k], origin="lower", cmap=DIV, vmin=-6, vmax=6)
        zz = np.concatenate([z[3:-3, 3:-3].ravel() for z in Z])
        a.set_title(f"{MLAB[m]}\n|Z|>5 px (4 passes): {int(np.nansum(np.abs(zz) > 5))}", fontsize=9)
    for a in ax:
        a.set_xticks([]); a.set_yticks([]); a.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label="per-pass residual SNR")
    fig.suptitle("E1 — pass-2 residual SNR maps of the stationary sky (17′×17′, no injections). Ideal = pure noise", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig3_residual_maps.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_e7():
    import pandas as pd
    p = os.path.join(R, E7DIR, "e7.csv")
    if not os.path.exists(p):
        return
    df = pd.read_csv(p)
    null = df[df.A == 0].copy()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    fap = 0.05
    for m, lab, c in [("B2", "B2 detector stack (λ-blind)", CAT[2]), ("B1", "B1 same-colour pairs", CAT[1]), ("C", "C spectral LOPO (ours)", CAT[3])]:
        thr = np.nanpercentile(null[m].abs(), 100 * (1 - fap))
        comp = (df[m] > thr).groupby(df.A).mean()
        ax[0].plot(comp.index[1:] * 100, comp.values[1:], "o-", color=c, label=lab, ms=6)
    ax[0].set(xscale="log", xlabel="injected brightening of one pass [%]", ylabel=f"completeness at FAP={fap:.0%}",
              title="Real stars: detecting true variability")
    ax[0].legend(fontsize=8)
    null["pf"] = null.nb1 / null.nk
    bins = [-0.01, 0.34, 0.67, 0.99, 1.01]
    labs = ["<1/3", "1/3–2/3", "2/3–<1", "all"]
    null["pfb"] = pd.cut(null.pf, bins, labels=labs)
    xs = np.arange(len(labs))
    for i, (lo, hi, nm) in enumerate([(10, 30, "SNR 10–30"), (30, 1e9, "SNR >30")]):
        s_ = null[(null.snr >= lo) & (null.snr < hi)]
        g = s_.groupby("pfb", observed=False).apply(lambda g: 1.4826 * np.nanmedian(np.abs(g.C)) / (1.4826 * np.nanmedian(np.abs(g.B1))))
        ax[1].plot(xs, g.values, "o-", color=CAT[i * 2 + 0 if i == 0 else 6], label=nm, ms=7)
    ax[1].axhline(1, color=INK2, lw=1, ls="--")
    ax[1].set_xticks(xs, labs)
    ax[1].set(xlabel="fraction of pass-k exposures that have a same-colour partner (|Δλ/λ|<1%)",
              ylabel="noise ratio  σ_C / σ_B1", title="When does the spectral model beat same-colour pairing?")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig6_variability_rule.png"), dpi=150)
    plt.close(fig)


def fig_completeness():
    import glob
    import pandas as pd
    fs = sorted(glob.glob(os.path.join(R, E3DIR, "inj_*.csv")))
    if not fs:
        return
    df = pd.concat([pd.read_csv(f) for f in fs])
    thr = {}
    for f in sorted(glob.glob(os.path.join(R, E3DIR, "null_*.npz"))):
        z = np.load(f)
        for k in z.files:
            thr.setdefault(k, []).append(z[k])
    thr = {k: np.nanpercentile(np.concatenate(v), 99) for k, v in thr.items()}
    out = {}
    for grp in ["all", "long"]:
        fig, ax = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
        for a, kind in zip(ax, ["star", "tdwarf", "ydwarf"]):
            sub = df[(df.kind == kind) & (df.pm >= 1.0)]
            for m in METH:
                key = f"{m}|{grp}"
                if key not in sub:
                    continue
                c = sub.assign(det=sub[key] > thr[key]).groupby("amp")["det"].mean()
                out[(grp, kind, m)] = c.to_dict()
                a.plot(c.index, c.values, "o-", color=MCOL[m], ms=6, label=MLAB[m])
            a.set_xscale("log")
            a.set(title=f"{kind} SED, μ≥1 px/yr", xlabel="injected peak surface brightness [MJy/sr]")
        ax[0].set_ylabel("completeness at equal false-alarm prob. (1%)")
        ax[0].legend(fontsize=8, loc="upper left")
        fig.suptitle(f"E3 — injection-recovery of moving point sources, detectors: {grp}", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"fig4_completeness_{grp}.png"), dpi=150)
        plt.close(fig)
    json.dump({"|".join(k): v for k, v in out.items()}, open(os.path.join(OUT, "completeness.json"), "w"), indent=1)


def fig_coverage():
    p = os.path.join(R, "e5", "e5_rows.json")
    if not os.path.exists(p):
        return
    import pandas as pd
    df = pd.DataFrame(json.load(open(p)))
    taus = [0.0025, 0.005, 0.01, 0.02, 0.04]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, d in enumerate(range(1, 7)):
        g = df[df.det == d]
        ax[0].plot(taus, [g[f"f_{t}"].median() for t in taus], "o-", color=CAT[i], label=f"D{d} measured")
        ax[0].plot(taus, [g[f"fpred_{t}"].median() for t in taus], "--", color=CAT[i], lw=1.2)
        ax[1].plot(taus, [np.median(np.clip(g[f"gain_{t}"], 0, 10)) for t in taus], "o-", color=CAT[i], label=f"D{d}")
    for a in ax:
        a.set_xscale("log"); a.set_xticks(taus); a.set_xticklabels(["0.25%", "0.5%", "1%", "2%", "4%"]); a.minorticks_off()
    ax[1].set_ylim(0.9, 5.2)
    ax[0].set(xlabel="same-colour tolerance τ = |Δλ/λ|", ylabel="fraction of exposures with a partner",
              title="Same-colour pairing usability (dashed: random-λ Poisson)")
    ax[1].set(xlabel="same-colour tolerance τ", ylabel="noise-limited SNR gain, C vs B1 (clipped at 5)", title="Predicted gain from using all exposures")
    ax[0].legend(fontsize=7, ncol=2)
    ax[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig5_coverage.png"), dpi=150)
    plt.close(fig)


def fig_relcolour():
    import pandas as pd
    p = os.path.join(R, E6DIR, "e6.json")
    if not os.path.exists(p):
        return
    df = pd.DataFrame(json.load(open(p)))
    df["obs_zp"] = df.obs - df.groupby(["det", "pas"]).obs.transform("median")
    df["bf"] = df.groupby("det").beta.transform("median")
    df["dmu"] = df.lin / df.beta.replace(0, np.nan)
    df["x"] = (df.beta - df.bf) * df.dmu
    rms = lambda v: float(np.sqrt(np.nanmean(np.asarray(v) ** 2)))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    lim = np.nanpercentile(np.abs(np.r_[df.obs_zp, df.x]), 99.5) * 1.1
    ax[0].plot([-lim, lim], [-lim, lim], color=INK2, lw=1, ls="--")
    for i, d in enumerate(range(1, 7)):
        m = df.det == d
        ax[0].scatter(df.x[m], df.obs_zp[m], s=12, color=CAT[i], label=f"D{d}", edgecolor=SURF, linewidth=0.4)
    r = df[["x", "obs_zp"]].corr().iloc[0, 1]
    ax[0].set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="(β_star − β_field) × Δ⟨ln λ⟩_pass", ylabel="stack false variability after zero-point match",
              title=f"Relative-colour law (r = {r:.2f})")
    ax[0].legend(fontsize=8, ncol=2)
    df["ad"] = (df.beta - df.bf).abs()
    q = pd.qcut(df.ad, 4)
    g = df.groupby(q, observed=False).apply(lambda g: pd.Series(dict(x=g.ad.median(), B2zp=rms(g.obs_zp), C=rms(g.obs - g.pred), B2=rms(g.obs))))
    ax[1].plot(g.x, g.B2, "o-", color=CAT[2], label="B2 stack, raw")
    ax[1].plot(g.x, g.B2zp, "o-", color=CAT[1], label="B2 stack + zero-point match")
    ax[1].plot(g.x, g.C, "o-", color=CAT[3], label="C spectral LOPO (ours)")
    ax[1].set(xlabel="|β_star − β_field|  (quartile medians)", ylabel="rms false variability (ln flux)",
              title="Spectral model pays off for colour-peculiar sources")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig7_relative_colour_law.png"), dpi=150)
    plt.close(fig)


def fig_e2():
    p = os.path.join(R, "e2_160", "e2_stamps.npz")
    if not os.path.exists(p):
        return
    st = np.load(p)
    info = json.load(open(os.path.join(R, "e2_160", "e2.json")))
    rowsdef = [("B2_detstack", "D123"), ("B2_detstack", "D56"), ("C_spectral", "D123"), ("C_spectral", "D56")]
    fig, ax = plt.subplots(len(rowsdef), 4, figsize=(8.4, 8.8))
    for i, (m, g) in enumerate(rowsdef):
        zs = info["z"][f"{m}|{g}"]
        for k in range(4):
            a = ax[i, k]
            key = f"{m}|{g}|{k}"
            if key in st.files:
                a.imshow(st[key], origin="lower", cmap=DIV, vmin=-6, vmax=6)
            a.set_xticks([]); a.set_yticks([]); a.grid(False)
            a.plot(8, 8, "+", color=INK, ms=10, mew=1)
            a.set_title(f"pass {k+1}: Z={zs[k]:.1f}", fontsize=8)
            if k == 0:
                a.set_ylabel(f"{MLAB[m].split(' (')[0]}\n{g}", fontsize=8)
    fig.suptitle("E2 — WISE 0855 (Y dwarf, 1.3 px/yr) in per-pass residuals; + = predicted position\n"
                 "D5/D6: dipole flips sign between pass 1 and pass 4 as the dwarf crosses its mean position (motion, not variability).\n"
                 "D1–D3 (where a Y dwarf is dark): spectral-model Z = 1.0 / 2.3 / 3.5 / 0.1 vs detector-stack 0.2 / 4.5 / 4.8 / 1.3", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig8_wise0855.png"), dpi=150)
    plt.close(fig)


def fig_multifield():
    p = os.path.join(R, "multifield.json")
    if not os.path.exists(p):
        return
    import pandas as pd
    t = pd.read_json(p)
    if "field" in t.columns:
        t = t.set_index("field")
    f = list(t.index)
    x = np.arange(len(f))
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for key, lab, c, mk in [("B2_raw", "B2 λ-blind stack", CAT[2], "o"), ("B2_zp", "B2 + zero-point", CAT[1], "s"), ("C", "C spectral (ours)", CAT[3], "D")]:
        ax[0].plot(x, 100 * t[key], mk, color=c, ms=7, label=lab)
    ax[0].set(xticks=x, xticklabels=f, ylabel="rms pass-to-pass false variability [%]", title="Stationary stars, all detectors (R2)")
    ax[0].legend(fontsize=8)
    for key, lab, c, mk in [("B2zp_peculiar", "B2 + zero-point", CAT[1], "s"), ("C_peculiar", "C spectral (ours)", CAT[3], "D")]:
        ax[1].plot(x, 100 * t[key], mk, color=c, ms=7, label=lab)
    ax[1].set(xticks=x, xticklabels=f, ylabel="rms false variability [%]", title="Colour-peculiar quartile |β−β_field| (R3)")
    ax[1].legend(fontsize=8)
    if "e9_periodic_raw" in t:
        ax[2].plot(x, 100 * t["e9_periodic_raw"], "o", color=CAT[1], ms=7, label="naive light curve")
        ax[2].plot(x, 100 * t["e9_periodic_sed"], "D", color=CAT[3], ms=7, label="÷ spectral model")
        ax[2].axhline(1, color=INK2, ls="--", lw=1)
        ax[2].set(xticks=x, xticklabels=f, ylabel="stationary stars flagged periodic [%]", title="Spurious periods at FAP<1% (R1b)")
        ax[2].legend(fontsize=8)
    for a in ax:
        a.set_xlabel("field")
    fig.suptitle(f"Rules tested across {len(f)} independent 17′ fields in the WISE 0855 footprint (uniform cut |β|<10)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig9_multifield.png"), dpi=150)
    plt.close(fig)


for f in [fig_problem, fig_b2law, fig_residual_maps, fig_completeness, fig_coverage, fig_e7, fig_relcolour, fig_e2, fig_multifield]:
    try:
        f()
        print("ok", f.__name__)
    except Exception as ex:
        import traceback
        traceback.print_exc()
