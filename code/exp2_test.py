"""Experiment 2: the exact Haar-conditional test (Fig. 2).

All p-values use hollow=True: the observed matrix and every Haar replicate
are diagonal-zeroed before statistics are evaluated, matching the pipeline
applied to the real (hollow) networks in exp3, so the level simulation
validates the procedure actually used on data.

Stages (checkpointed for constrained environments):
  level0 level1 | loc0 loc1 | graphon0 graphon1 | hetero0 hetero1 | plot
Chunk seeds: base + 10 * chunk_index (documented for reproducibility).

Panel (a): validity: ECDF of conditional p-values over 250 null draws
           (spiked GOE, random spikes).
Panel (b): power against three exchangeable, non-rotatable alternatives.
"""
import sys, json, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import (mpl_style, spiked_goe, goe, haar, conditional_pvalues,
                    PAL, SERIES, save_json)

mpl_style()
STAGE = sys.argv[1] if len(sys.argv) > 1 else "plot"
if STAGE == "all":
    import subprocess
    for _st in ['level0', 'level1', 'loc0', 'loc1', 'graphon0', 'graphon1', 'hetero0', 'hetero1', 'plot']:
        subprocess.run([sys.executable, __file__, _st], check=True)
    sys.exit(0)
PART = "../results/exp2_parts"
os.makedirs(PART, exist_ok=True)
res = {}
n, B = 120, 99
WHICH = ["ipr", "evks", "clust", "cv", "kurt"]
LAB = {"ipr": "localization (IPR)", "evks": "eigenvector KS",
       "clust": "weighted clustering", "cv": "strength CV",
       "kurt": "entry kurtosis"}
R = 250          # null replicates (2 chunks of 125)
Rpow = 50        # replicates per power point (2 chunks of 25)

# ---------------------------------------------------------------- level
if STAGE.startswith("level"):
    ci = int(STAGE[-1]); m = R // 2
    rng = np.random.default_rng(20260812 + 10 * ci)
    pn = {s: np.empty(m) for s in WHICH}
    for r in range(m):
        k = rng.integers(0, 4)
        th = np.sort(rng.uniform(1.2, 3.0, size=k))[::-1] if k else []
        X = spiked_goe(n, th, 1.0, rng)
        pv, _ = conditional_pvalues(X, B, rng, which=WHICH, hollow=True)
        for s in WHICH:
            pn[s][r] = pv[s]
    save_json(pn, f"{PART}/level_c{ci}.json")
    print(f"level chunk {ci} done")
    sys.exit(0)

# ---------------------------------------------------------------- power
def alt_localized(n, eta, rng):
    m = int(np.ceil(n ** (2 / 3)))
    v_loc = np.zeros(n); idx = rng.choice(n, m, replace=False)
    v_loc[idx] = rng.standard_normal(m); v_loc /= np.linalg.norm(v_loc)
    v_del = rng.standard_normal(n); v_del /= np.linalg.norm(v_del)
    v = (1 - eta) * v_del + eta * v_loc; v /= np.linalg.norm(v)
    return spiked_goe(n, [1.8], 1.0, rng, vecs=v[:, None])

def alt_graphon(n, rho, rng):
    U = rng.uniform(0, 1, n)
    F = np.stack([np.cos(np.pi * U), np.cos(2 * np.pi * U)], axis=1)
    Wg = F @ np.diag([1.0, 0.5]) @ F.T
    return (1.0 / np.sqrt(n)) * goe(n, rng) + 3.6 * rho * Wg / n

def alt_hetero(n, eta, rng):
    s = np.exp(eta * rng.standard_normal(n))
    s /= np.sqrt(np.mean(s ** 2))
    S = np.sqrt(np.outer(s, s))
    W = goe(n, rng) * S / np.sqrt(n)
    return W + spiked_goe(n, [1.8], 0.0, rng)

ALTS = {"loc": ("localized spike", alt_localized,
                np.linspace(0, 1.0, 5), 20260814),
        "graphon": ("graphon signal", alt_graphon,
                    np.linspace(0, 1.0, 5), 20260815),
        "hetero": ("heteroskedastic", alt_hetero,
                   np.linspace(0, 0.8, 5), 20260816)}

for stub, (name, fn, etas, seed) in ALTS.items():
    if STAGE.startswith(stub):
        ci = int(STAGE[-1]); m = Rpow // 2
        rng = np.random.default_rng(seed + 10 * ci)
        rej = np.zeros(len(etas))
        for i, eta in enumerate(etas):
            for r in range(m):
                X = fn(n, eta, rng)
                pv, _ = conditional_pvalues(X, B, rng, which=WHICH,
                                            hollow=True)
                # Bonferroni min-p across the five statistics
                if min(pv.values()) <= 0.05 / len(WHICH):
                    rej[i] += 1
        save_json({"etas": etas, "rej": rej, "m": m},
                  f"{PART}/{stub}_c{ci}.json")
        print(f"{name} chunk {ci} rejections:", rej)
        sys.exit(0)

# ---------------------------------------------------------------- plot
pnull = {s: np.concatenate([
    np.array(json.load(open(f"{PART}/level_c{c}.json"))[s])
    for c in range(2)]) for s in WHICH}
res["null_p"] = pnull
res["hollow"] = True
power = {}
for stub, (name, fn, etas, seed) in ALTS.items():
    parts = [json.load(open(f"{PART}/{stub}_c{c}.json")) for c in range(2)]
    rej = np.array(parts[0]["rej"]) + np.array(parts[1]["rej"])
    mm = parts[0]["m"] + parts[1]["m"]
    power[name] = {"etas": np.array(parts[0]["etas"]), "power": rej / mm}
res["power"] = power

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
ax.set_title("(a) validity under $H_0$, hollow pipeline (250 draws)",
             loc="left")
ax.legend(loc="lower right")

ax = axes[1]
order = ["localized spike", "heteroskedastic", "graphon signal"]
for i, name in enumerate(order):
    d = power[name]
    ax.plot(d["etas"], d["power"], "o-", color=SERIES[i], label=name, ms=3.4)
ax.axhline(0.05, color=PAL["muted"], lw=0.8, ls="--")
ax.annotate(r"$\alpha=0.05$", xy=(0.02, 0.085), fontsize=7.5,
            color=PAL["ink2"])
ax.annotate("low power against weak\nsmooth graphons (see text)",
            xy=(0.72, 0.16), fontsize=7, color=SERIES[2], ha="center")
ax.set_xlabel(r"departure from rotatability $\eta$")
ax.set_ylabel("rejection rate")
ax.set_ylim(-0.03, 1.05)
ax.set_title("(b) power of the combined test ($n=120$, $B=99$)", loc="left")
ax.legend(loc="center left", fontsize=7)

fig.tight_layout()
fig.savefig("../figures/fig2_test.pdf", bbox_inches="tight")
save_json(res, "../results/exp2.json")
print("exp2 assembled; combined null rejection at .05:",
      np.mean(np.minimum.reduce([pnull[s] for s in WHICH])
              <= 0.05 / len(WHICH)))
