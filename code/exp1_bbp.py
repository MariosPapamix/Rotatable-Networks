"""Experiment 1: finite-n calibration of the ergodic rotatable law (Fig. 1).

Panel (a): top eigenvalue of the spiked ensemble vs the BBP prediction
           theta + sigma^2/theta (theta > sigma), bulk edge 2 sigma below.
Panel (b): nested corners of single fixed realizations of the ergodic
           array (fixed lambda in l^2, sigma = 1): one infinite-array draw
           X_ij = W_ij + sum_k lam_k (xi_ki xi_kj - delta_ij) is generated
           at n_max = 1600 and its leading principal n x n submatrices,
           divided by sqrt(n), are diagonalized as n grows.  The effective
           spikes s_k = sqrt(n) lam_k cross the detection threshold one by
           one; curves average 8 independent array realizations.
"""
import sys, os, json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import mpl_style, spiked_goe, goe, PAL, SERIES, save_json

mpl_style()
STAGE = sys.argv[1] if len(sys.argv) > 1 else "all"
if STAGE == "all":
    import subprocess
    for _st in ['a0', 'a1', 'a2', 'b', 'plot']:
        subprocess.run([sys.executable, __file__, _st], check=True)
    sys.exit(0)
PART = "../results/exp1_parts"
os.makedirs(PART, exist_ok=True)
res = {}

# ---------------------------------------------------------------- panel (a)
n = 2000
sigma = 1.0
thetas = np.linspace(0.1, 3.0, 25)
reps = 6
CHUNKS = [(0, 9), (9, 17), (17, 25)]
if STAGE.startswith("a"):
    ci = int(STAGE[1])
    lo, hi = CHUNKS[ci]
    rng = np.random.default_rng(20260812 + ci)
    emp = np.zeros((hi - lo, reps))
    for i, th in enumerate(thetas[lo:hi]):
        for r in range(reps):
            X = spiked_goe(n, [th], sigma, rng)
            emp[i, r] = np.linalg.eigvalsh(X)[-1]
    save_json({"emp": emp}, f"{PART}/a{ci}.json")
    print(f"panel a chunk {ci} done")
    sys.exit(0)
if STAGE in ("all", "plot"):
    emp = np.vstack([np.array(json.load(open(f"{PART}/a{c}.json"))["emp"])
                     for c in range(3)])
    grid = np.linspace(0.05, 3.0, 400)
    pred = np.where(grid > sigma,
                    grid + sigma ** 2 / np.maximum(grid, 1e-9), 2 * sigma)
    res["panel_a"] = {"thetas": thetas, "emp_mean": emp.mean(1),
                     "emp_sd": emp.std(1)}

# ---------------------------------------------------------------- panel (b)
# TRUE nested corners: one fixed array realization per replicate.
lams = np.array([0.28, 0.14, 0.07])          # fixed array-level spectrum
sigma = 1.0
nmax = 1600
ns = np.array([25, 50, 100, 200, 400, 800, 1600])
reps = 8
if STAGE == "b":
    top4 = np.zeros((len(ns), 4, reps))
    rng = np.random.default_rng(20260830)
    for r in range(reps):
        W = goe(nmax, rng)                                # infinite-GOE corner
        Xi = rng.standard_normal((len(lams), nmax))       # Gaussian factors
        A = sigma * W + (Xi.T * lams) @ Xi - lams.sum() * np.eye(nmax)
        for j, nn in enumerate(ns):
            ev = np.linalg.eigvalsh(A[:nn, :nn] / np.sqrt(nn))
            top4[j, :, r] = ev[-4:][::-1]
    save_json({"top4": top4}, f"{PART}/b.json")
    print("panel b done")
    sys.exit(0)
top4 = np.array(json.load(open(f"{PART}/b.json"))["top4"])
res["panel_b"] = {"ns": ns, "lams": lams, "mean": top4.mean(2),
                  "sd": top4.std(2), "nested": True, "reps": reps}

# Panel (a): standard BBP check -> Supplement Figure S2.
figA, ax = plt.subplots(figsize=(3.6, 2.6))
ax.plot(grid, pred, color=PAL["ink"], lw=1.4, label="BBP prediction")
ax.axhline(2 * 1.0, color=PAL["axis"], lw=0.8, ls=":")
ax.errorbar(thetas, emp.mean(1), yerr=emp.std(1), fmt="o", ms=3.2,
            color=SERIES[0], ecolor=SERIES[0], elinewidth=0.9, capsize=0,
            label=r"simulation ($n=2000$)")
ax.axvline(1.0, color=PAL["muted"], lw=0.8, ls="--")
ax.annotate(r"$\theta=\sigma$", xy=(1.0, 3.35), fontsize=7.5,
            color=PAL["ink2"], ha="left", xytext=(1.07, 3.35))
ax.annotate("bulk edge $2\\sigma$", xy=(2.45, 2.04), fontsize=7.5,
            color=PAL["ink2"])
ax.set_xlabel(r"spike strength $\theta$")
ax.set_ylabel(r"largest eigenvalue of $X_n$")
ax.set_title("BBP prediction in the ergodic rotatable law", loc="left",
             fontsize=8.5)
ax.legend(loc="upper left")
figA.tight_layout()
figA.savefig("../figures/figS2_bbp_check.pdf", bbox_inches="tight")

# Panel (b): the fixed-array reading -> main-text Figure 1.
fig, ax = plt.subplots(figsize=(4.6, 2.8))
cols = [SERIES[0], SERIES[1], SERIES[2]]
for k in range(3):
    s = np.sqrt(ns) * lams[k]
    th_pred = np.where(s > 1, s + 1 / np.maximum(s, 1e-9), 2.0)
    ax.plot(ns, th_pred, color=cols[k], lw=1.0, ls="--")
    ax.plot(ns, top4[:, k, :].mean(1), "o-", ms=3.2, color=cols[k], lw=1.4,
            label=rf"$\lambda_{k+1}={lams[k]:.2f}$")
ax.plot(ns, top4[:, 3, :].mean(1), "o-", ms=2.6, color=PAL["muted"],
        lw=1.0, label="bulk")
ax.axhline(2.0, color=PAL["axis"], lw=0.8, ls=":")
ax.set_xscale("log")
ax.set_xlabel(r"observed nodes $n$")
ax.set_ylabel(r"top eigenvalues of $X_n/\sqrt{n}$")
ax.set_title("nested corners of a fixed array", loc="left", fontsize=8.5)
ax.legend(loc="upper left", ncol=1)

fig.tight_layout()
fig.savefig("../figures/fig1_bbp.pdf", bbox_inches="tight")
save_json(res, "../results/exp1.json")
print("exp1 done (nested corners)")
