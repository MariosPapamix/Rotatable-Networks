"""Experiment 3: real data (Fig. 3 + Table 1).

(a) HCP group-average functional connectome (Schaefer-200 parcellation,
    ENIGMA toolbox): Haar-conditional null densities vs observed values;
    conditional residual z-scores are also saved (obs - null mean)/null sd.
(b) S&P 500 rolling-window return-correlation matrices: conditional
    p-values across windows x statistics, raw and after Gaussianization.
    If data/all_stocks_5yr.csv is absent (see data/README.md), cached
    values from results/exp3.json are reused and flagged.
(c) Air-traffic network for Table 1; the Gaussianized analysis is
    repeated over 5 tie-breaking seeds (many tied structural zeros).
"""
import json, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import (mpl_style, conditional_pvalues, fisher_z, gaussianize,
                    all_stats, top_vecs, PAL, SERIES, save_json,
                    excess_kurtosis, hill_estimator)

mpl_style()
rng = np.random.default_rng(20260812)
res = {}
WHICH = ["ipr", "clust", "mod", "cv", "kurt"]
LAB = {"ipr": "IPR", "clust": "clustering", "mod": "modularity",
       "cv": "strength CV", "kurt": "kurtosis"}

# ------------------------------------------------------------- (a) HCP FC
FC = np.loadtxt("../data/funcMatrix_ctx_schaefer_200.csv", delimiter=",")
Z = fisher_z(FC)
pv_h, obs_h, null_h = conditional_pvalues(Z, 499, rng, which=WHICH,
                                          return_null=True, hollow=True)
res["hcp_fc"] = {"p": pv_h, "obs": obs_h,
                 "null_mean": {s: float(np.mean(null_h[s])) for s in WHICH},
                 "null_sd": {s: float(np.std(null_h[s])) for s in WHICH},
                 "z": {s: float((obs_h[s] - np.mean(null_h[s]))
                                / max(np.std(null_h[s]), 1e-12))
                       for s in WHICH}}
print("HCP FC p-values:", pv_h)
print("HCP FC residual z:", res["hcp_fc"]["z"])

# Gaussianized HCP (no ties generically; seed for reproducibility)
Zg = gaussianize(Z, rng=np.random.default_rng(101))
pv_hg, _ = conditional_pvalues(Zg, 499, rng, which=WHICH, hollow=True)
res["hcp_fc_gauss"] = {"p": pv_hg}

# ------------------------------------------------------- (b) S&P 500 rolls
Bwin = 99
WIN_WHICH = ["ipr", "clust", "cv", "kurt"]
csv_path = "../data/all_stocks_5yr.csv"
if os.path.exists(csv_path):
    px = pd.read_csv(csv_path)
    piv = px.pivot_table(index="date", columns="Name", values="close")
    piv = piv.dropna(axis=1)                  # complete histories only
    ret = np.log(piv).diff().dropna()
    names = ret.columns
    win, stride = 120, 42
    starts = range(0, len(ret) - win + 1, stride)
    sub = rng.choice(len(names), size=min(300, len(names)), replace=False)
    R = ret.values[:, sub]
    pmat, pmat_g, dates = [], [], []
    for s0 in starts:
        Rw = R[s0:s0 + win]
        Rw = (Rw - Rw.mean(0)) / Rw.std(0)
        C = np.corrcoef(Rw.T)
        Zw = fisher_z(C)
        pv, _ = conditional_pvalues(Zw, Bwin, rng, which=WIN_WHICH,
                                    hollow=True)
        pv_g, _ = conditional_pvalues(gaussianize(Zw, rng=rng), Bwin, rng,
                                      which=WIN_WHICH, hollow=True)
        pmat.append([pv[s] for s in WIN_WHICH])
        pmat_g.append([pv_g[s] for s in WIN_WHICH])
        dates.append(str(ret.index[s0 + win - 1])[:10])
    pmat, pmat_g = np.array(pmat), np.array(pmat_g)
    res["sp500"] = {"which": WIN_WHICH, "dates": dates, "p": pmat,
                    "p_gauss": pmat_g, "n": int(len(sub)),
                    "n_windows": int(pmat.shape[0]), "cached": False}
else:
    print("S&P csv absent; reusing cached results/exp3.json sp500 block")
    cached = json.load(open("../results/exp3.json"))["sp500"]
    WIN_WHICH = cached["which"]
    pmat = np.array(cached["p"]); pmat_g = np.array(cached["p_gauss"])
    dates = cached["dates"]
    res["sp500"] = dict(cached, cached=True)
print("S&P windows:", pmat.shape, "median p per stat:", np.median(pmat, 0))

# ------------------------------------------- (c) air-traffic network stats
fl = pd.read_csv("../data/flights_airport.csv")
airports = pd.unique(pd.concat([fl.origin, fl.destination]))
amap = {a: i for i, a in enumerate(airports)}
nA = len(airports)
A = np.zeros((nA, nA))
for o, d, c in zip(fl.origin, fl.destination, fl["count"]):
    A[amap[o], amap[d]] += c
A = A + A.T                       # symmetrize flows
deg = (A > 0).sum(1)
keep = deg >= 5                   # trim isolated airfields
A = A[np.ix_(keep, keep)]
nA = A.shape[0]
np.fill_diagonal(A, 0.0)
Alog = np.log1p(A)
pv_f, obs_f = conditional_pvalues(Alog, 199, rng, which=WHICH, hollow=True)
pv_fg_seeds = []
for s in range(5):
    Ag = gaussianize(Alog, rng=np.random.default_rng(200 + s))
    pv_fg, _ = conditional_pvalues(Ag, 199, rng, which=WHICH, hollow=True)
    pv_fg_seeds.append(pv_fg)
res["flights"] = {"n": int(nA), "p_log": pv_f,
                  "p_gauss_seeds": pv_fg_seeds,
                  "p_gauss": pv_fg_seeds[0],
                  "hill": hill_estimator(A[np.triu_indices(nA, 1)]),
                  "kurt": excess_kurtosis(Alog)}
print("flights n:", nA, "p:", pv_f)
print("flights Gaussianized (5 tie seeds):", pv_fg_seeds)

# --------------------------------------------------------------------- fig
fig = plt.figure(figsize=(7.0, 4.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.55, wspace=0.35)

show = ["ipr", "clust", "cv"]
for j, s in enumerate(show):
    ax = fig.add_subplot(gs[0, j])
    x = null_h[s]
    ax.hist(x, bins=24, color=SERIES[0], alpha=0.55, density=True,
            edgecolor="white", linewidth=0.3)
    ax.axvline(obs_h[s], color=SERIES[1], lw=1.6)
    lo, hi = min(x.min(), obs_h[s]), max(x.max(), obs_h[s])
    pad = 0.08 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_title(f"({'abc'[j]}) HCP: {LAB[s]}  $p={pv_h[s]:.3f}$", loc="left")
    ax.set_yticks([])
    if j == 0:
        ax.set_ylabel("Haar-conditional null")
    ax.annotate("observed", xy=(obs_h[s], ax.get_ylim()[1] * 0.9),
                color=SERIES[1], fontsize=7,
                ha="left" if obs_h[s] < (lo + hi) / 2 else "right")

ax = fig.add_subplot(gs[1, :])
M = np.vstack([pmat.T, pmat_g.T])
im = ax.imshow(M, aspect="auto", cmap="Blues_r", vmin=0, vmax=1,
               interpolation="nearest")
ax.axhline(len(WIN_WHICH) - 0.5, color="white", lw=2.0)
ax.set_yticks(range(2 * len(WIN_WHICH)),
              [LAB[s] for s in WIN_WHICH] + [LAB[s] + " (g)" for s in WIN_WHICH],
              fontsize=7)
step = max(1, len(dates) // 9)
ax.set_xticks(range(0, len(dates), step),
              [dates[i] for i in range(0, len(dates), step)],
              fontsize=6.5, rotation=0)
ax.set_title("(d) S\\&P 500 rolling windows: conditional $p$-values "
             "(raw / after Gaussianization ``(g)'')", loc="left")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.set_label("$p$", rotation=0, labelpad=6)
ax.grid(False)

fig.savefig("../figures/fig3_real.pdf", bbox_inches="tight")
save_json(res, "../results/exp3_new.json")
print("exp3 done")
