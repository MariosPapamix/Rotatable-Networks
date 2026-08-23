"""Experiment 5: heavy tails and the tail-versus-structure diagnostic
(Fig. 5).

Four networks are pushed through the before/after-Gaussianization
readout of Section 6:
  (1) spiked GOE (Gaussian control);
  (2) the Theorem 4 spiked p-stable law itself, p = 1.5: SaS Wigner
      noise with the canonical diagonal plus SaS-factor outer products
      (the factors are heavy tailed, hence l4-localized);
  (3) stable Wigner noise plus delocalized (Haar-vector) spikes: heavy
      marginal, generic eigenvector geometry;
  (4) the real air-traffic flow network (log(1+flow)); its Gaussianized
      analysis is repeated over 5 tie-breaking seeds.
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import (mpl_style, spiked_goe, stable_wigner, spiked_stable,
                    haar, gaussianize, conditional_pvalues, hill_estimator,
                    offdiag, PAL, SERIES, save_json)

mpl_style()
rng = np.random.default_rng(20260812)
res = {}
n = 200
alpha = 1.5

# four networks ---------------------------------------------------------------
Xg = spiked_goe(n, [2.5, 1.6], 1.0, rng)

Xs_fac = spiked_stable(n, [6.0, 3.0], 1.0, alpha, rng) / n ** (1 / alpha)

Zs = stable_wigner(n, alpha, rng) / n ** (1 / alpha)
U = haar(n, rng)[:, :2]
Xs_del = Zs + (U * np.array([3.0, 2.0])) @ U.T      # delocalized spikes

fl = pd.read_csv("../data/flights_airport.csv")
airports = pd.unique(pd.concat([fl.origin, fl.destination]))
amap = {a: i for i, a in enumerate(airports)}
nA = len(airports)
A = np.zeros((nA, nA))
for o, d, c in zip(fl.origin, fl.destination, fl["count"]):
    A[amap[o], amap[d]] += c
A = A + A.T
keep = (A > 0).sum(1) >= 5
A = A[np.ix_(keep, keep)]
np.fill_diagonal(A, 0.0)

# (a) tails ------------------------------------------------------------------
def ccdf(x):
    x = np.sort(np.abs(x[x != 0]))
    return x, 1.0 - np.arange(1, len(x) + 1) / (len(x) + 1)

nets = {
    "Gaussian (spiked GOE)": (offdiag(Xg) * np.sqrt(n), SERIES[0]),
    rf"spiked $1.5$-stable (Thm 4)": (offdiag(Xs_fac), SERIES[1]),
    "air traffic (flows)": (A[np.triu_indices(A.shape[0], 1)], SERIES[2]),
}
hills = {k: hill_estimator(v) for k, (v, c) in nets.items()}
res["hill"] = hills
res["hill_delocalized"] = hill_estimator(offdiag(Xs_del))

# (b) diagnostics ------------------------------------------------------------
diag = {}
cases = [("Gaussian", Xg), ("stable, SaS factors", Xs_fac),
         ("stable, deloc. spikes", Xs_del), ("air traffic", np.log1p(A))]
for key, X in cases:
    # synthetic ensembles carry their model diagonals (Theorem 4 specifies
    # them); only the genuinely hollow air-traffic data use the hollow
    # pipeline.  For extremely heavy spectra, zeroing the diagonal of a
    # replicate itself localizes its small-eigenvalue eigenvectors, so the
    # full-matrix convention keeps the IPR readout interpretable here.
    hol = (key == "air traffic")
    pv_raw, _ = conditional_pvalues(X, 199, rng, which=["kurt", "ipr"],
                                    hollow=hol)
    if key == "air traffic":
        seeds = []
        for s in range(5):
            Xg_ = gaussianize(X, rng=np.random.default_rng(300 + s))
            pv_g, _ = conditional_pvalues(Xg_, 199, rng,
                                          which=["kurt", "ipr"], hollow=hol)
            seeds.append(pv_g)
        diag[key] = {"raw": pv_raw, "gauss": seeds[0], "gauss_seeds": seeds}
    else:
        Xg_ = gaussianize(X, rng=np.random.default_rng(300))
        pv_g, _ = conditional_pvalues(Xg_, 199, rng,
                                      which=["kurt", "ipr"], hollow=hol)
        diag[key] = {"raw": pv_raw, "gauss": pv_g}
    print(key, "raw:", pv_raw, "| gaussianized:", diag[key]["gauss"])
res["diag"] = diag

# fig -------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.7))
ax = axes[0]
for k, (v, c) in nets.items():
    x, s = ccdf(v / np.median(np.abs(v[v != 0])))
    ax.loglog(x, s, color=c, lw=1.4,
              label=f"{k}  ($\\hat\\alpha_{{\\rm Hill}}={hills[k]:.2f}$)")
ax.set_xlim(left=3e-3)
ax.set_xlabel("normalized $|$weight$|$")
ax.set_ylabel("CCDF")
ax.set_title("(a) edge-weight tails", loc="left")
ax.legend(loc="lower left", fontsize=6.6)

ax = axes[1]
keys = [k for k, _ in cases]
xpos = np.arange(len(keys))
w = 0.19
for i, (stat, off) in enumerate([("kurt", -1.5), ("ipr", -0.5)]):
    ax.bar(xpos + off * w, [diag[k]["raw"][stat] for k in keys], w,
           color=SERIES[i], label=f"{stat.upper()} raw")
for i, (stat, off) in enumerate([("kurt", 0.5), ("ipr", 1.5)]):
    ax.bar(xpos + off * w, [diag[k]["gauss"][stat] for k in keys], w,
           color=SERIES[i], alpha=0.45, hatch="//",
           label=f"{stat.upper()} Gaussianized")
ax.axhline(0.05, color=PAL["muted"], lw=0.8, ls="--")
ax.annotate(r"$\alpha=0.05$", xy=(2.52, 0.065), fontsize=7, color=PAL["ink2"])
ax.set_xticks(xpos, keys, fontsize=6.4)
ax.set_ylabel("conditional $p$-value")
ax.set_ylim(0, 1.04)
ax.set_title("(b) tail-driven vs structural violations", loc="left")
ax.legend(loc="upper right", fontsize=6.2, ncol=2)

fig.tight_layout()
fig.savefig("../figures/fig5_heavy.pdf", bbox_inches="tight")
save_json(res, "../results/exp5.json")
print("exp5 done")
