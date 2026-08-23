"""Re-render Fig. 2 from results/exp2.json (no re-simulation)."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import mpl_style, PAL, SERIES

mpl_style()
res = json.load(open("../results/exp2.json"))
WHICH = ["ipr", "evks", "clust", "cv", "kurt"]
LAB = {"ipr": "localization (IPR)", "evks": "eigenvector KS",
       "clust": "weighted clustering", "cv": "strength CV",
       "kurt": "entry kurtosis"}
pnull = {s: np.array(v) for s, v in res["null_p"].items()}
power = res["power"]

fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
ax = axes[0]
u = np.linspace(0, 1, 200)
ax.plot(u, u, color=PAL["axis"], lw=0.9, ls=":")
for i, s in enumerate(WHICH):
    p = np.sort(pnull[s])
    ax.step(np.concatenate([[0], p]), np.arange(len(p) + 1) / len(p),
            where="post", color=SERIES[i], label=LAB[s], lw=1.3)
ax.set_xlabel("conditional $p$-value")
ax.set_ylabel("empirical CDF")
ax.set_title("(a) validity under $H_0$, hollow pipeline (250 draws)", loc="left")
ax.legend(loc="lower right")

ax = axes[1]
order = ["localized spike", "heteroskedastic", "graphon signal"]
for i, name in enumerate(order):
    d = power[name]
    ax.plot(d["etas"], d["power"], "o-", color=SERIES[i], label=name, ms=3.4)
ax.axhline(0.05, color=PAL["muted"], lw=0.8, ls="--")
ax.annotate(r"$\alpha=0.05$", xy=(0.02, 0.085), fontsize=7.5, color=PAL["ink2"])
ax.annotate("smooth graphon $\\approx$ rotatable\n(near size; see text)",
            xy=(0.72, 0.13), fontsize=7, color=SERIES[2], ha="center")
ax.set_xlabel(r"departure from rotatability $\eta$")
ax.set_ylabel("rejection rate")
ax.set_ylim(-0.03, 1.05)
ax.set_title("(b) power of the combined test ($n=120$, $B=99$)", loc="left")
ax.legend(loc="center left", fontsize=7)

fig.tight_layout()
fig.savefig("../figures/fig2_test.pdf", bbox_inches="tight")
print("fig2 replotted")
