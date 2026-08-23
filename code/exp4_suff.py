"""Experiment 4: spectral sufficiency, conditional residuals, and the
cohort protocol on synthetic cohorts calibrated to connectomic dimensions
(Fig. 4).  Cohorts are synthetic; the ABIDE replication protocol is in
Appendix B.

Stages (checkpointed): for K in {cohortA, cohortB}:
  K_gen | K_res0 | K_res1 | K_test ; then: finish
"""
import json, os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.linalg import eigh
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from common import (mpl_style, spiked_goe, all_stats, top_vecs, fisher_z,
                    haar, conditional_pvalues, PAL, SERIES, save_json)

mpl_style()
res = {}
n, S = 200, 240
NS = ["ipr", "clust", "mod", "cv", "kurt"]
STAGE = sys.argv[1] if len(sys.argv) > 1 else "finish"
if STAGE == "all":
    import subprocess
    for _st in ['cohortA_gen', 'cohortA_res0', 'cohortA_res1', 'cohortA_test', 'cohortB_gen', 'cohortB_res0', 'cohortB_res1', 'cohortB_test', 'finish']:
        subprocess.run([sys.executable, __file__, _st], check=True)
    sys.exit(0)
PART = "../results/exp4_parts"
os.makedirs(PART, exist_ok=True)
BRES = 32                       # Haar draws per subject for residual means
SEEDS = {"cohortA": {"gen": 20260812, "res": 20260813, "test": 20260840},
         "cohortB": {"gen": 20260822, "res": 20260823, "test": 20260841}}

def spectral_features(X, k=10):
    w = np.sort(eigh(X)[0])
    top = w[-k:][::-1]; bot = w[:3]
    return np.concatenate([top, bot,
                           [np.mean(w ** 2), np.mean(w ** 3), np.mean(w ** 4)]])

def nonspectral_features(X):
    _, V = top_vecs(X)
    st = all_stats(X, V=V, which=NS)
    return np.array([st[s] for s in NS])

def make_subject_A(g, rng):
    d = 0.55 * g
    th = np.array([2.4 + d + 0.25 * rng.standard_normal(),
                   1.7 + 0.20 * rng.standard_normal(),
                   1.25 - 0.5 * d + 0.15 * rng.standard_normal()])
    th = np.maximum(th, 1.05)
    return spiked_goe(n, th, 1.0, rng)

def make_subject_B(g, rng):
    th = np.array([2.4 + 0.25 * rng.standard_normal(),
                   1.7 + 0.20 * rng.standard_normal()])
    th = np.maximum(th, 1.05)
    if g == 0:
        return spiked_goe(n, th, 1.0, rng)
    m = 20
    v1 = np.zeros(n); idx = rng.choice(n, m, replace=False)
    v1[idx] = rng.standard_normal(m); v1 /= np.linalg.norm(v1)
    v2 = rng.standard_normal(n); v2 -= v1 * (v1 @ v2)
    v2 /= np.linalg.norm(v2)
    return spiked_goe(n, th, 1.0, rng, vecs=np.stack([v1, v2], 1))

MAKERS = {"cohortA": make_subject_A, "cohortB": make_subject_B}
y = np.repeat([0, 1], S // 2)

if STAGE != "finish":
    key, sub = STAGE.split("_", 1)
    if sub == "gen":
        rng = np.random.default_rng(SEEDS[key]["gen"])
        mats = np.stack([MAKERS[key](g, rng) for g in y])
        np.save(f"{PART}/{key}_mats.npy", mats)
        Xs = np.array([spectral_features(M) for M in mats])
        Xn = np.array([nonspectral_features(M) for M in mats])
        np.save(f"{PART}/{key}_Xs.npy", Xs)
        np.save(f"{PART}/{key}_Xn.npy", Xn)
        print(key, "generated; feature shapes", Xs.shape, Xn.shape)
    elif sub.startswith("res"):
        ci = int(sub[-1]); lo, hi = (0, S // 2) if ci == 0 else (S // 2, S)
        mats = np.load(f"{PART}/{key}_mats.npy")
        Xn = np.load(f"{PART}/{key}_Xn.npy")
        rng = np.random.default_rng(SEEDS[key]["res"] + 100 * ci)
        Xr = np.empty((hi - lo, len(NS)))
        for i in range(lo, hi):
            w = eigh(mats[i])[0]
            nulls = np.empty((BRES, len(NS)))
            for b in range(BRES):
                U = haar(n, rng)
                Xb = (U * w) @ U.T
                nulls[b] = nonspectral_features(Xb)
            Xr[i - lo] = Xn[i] - nulls.mean(0)
        np.save(f"{PART}/{key}_Xr{ci}.npy", Xr)
        print(key, f"residual chunk {ci} done")
    elif sub == "test":
        mats = np.load(f"{PART}/{key}_mats.npy")
        rng = np.random.default_rng(SEEDS[key]["test"])
        sel = np.concatenate([np.arange(20), np.arange(S // 2, S // 2 + 20)])
        rej = {0: [], 1: []}
        for i in sel:
            pv, _ = conditional_pvalues(mats[i], 99, rng, which=NS)
            rej[int(y[i])].append(min(pv.values()) <= 0.05 / len(NS))
        out = {"rej_g0": float(np.mean(rej[0])),
               "rej_g1": float(np.mean(rej[1]))}
        save_json(out, f"{PART}/{key}_test.json")
        print(key, "per-subject rejection rates:", out)
    sys.exit(0)

# ------------------------------------------------------------------ finish
def cv_auc(F, yy):
    skf = StratifiedKFold(5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(F, yy):
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000, C=1.0))
        clf.fit(F[tr], yy[tr])
        aucs.append(roc_auc_score(yy[te], clf.predict_proba(F[te])[:, 1]))
    return (float(np.mean(aucs)), float(np.std(aucs) / np.sqrt(5)))

FEATSETS = ["spectral", "non-spectral", "residual", "both"]
for key in ("cohortA", "cohortB"):
    Xs = np.load(f"{PART}/{key}_Xs.npy")
    Xn = np.load(f"{PART}/{key}_Xn.npy")
    Xr = np.vstack([np.load(f"{PART}/{key}_Xr0.npy"),
                    np.load(f"{PART}/{key}_Xr1.npy")])
    res[key] = {"spectral": cv_auc(Xs, y), "non-spectral": cv_auc(Xn, y),
                "residual": cv_auc(Xr, y),
                "both": cv_auc(np.hstack([Xs, Xn]), y)}
    res[key + "_test"] = json.load(open(f"{PART}/{key}_test.json"))
    print(key, res[key]); print(key, "tests:", res[key + "_test"])

# ------------------------------------------------------- (c) S&P windows
csv_path = "../data/all_stocks_5yr.csv"
if os.path.exists(csv_path):
    rng = np.random.default_rng(20260812)
    px = pd.read_csv(csv_path)
    piv = px.pivot_table(index="date", columns="Name",
                         values="close").dropna(axis=1)
    ret = np.log(piv).diff().dropna()
    sub = rng.choice(ret.shape[1], size=min(300, ret.shape[1]), replace=False)
    R = ret.values[:, sub]
    win, stride = 120, 21
    feats, targ = [], []
    for s0 in range(0, len(R) - win + 1, stride):
        Rw = R[s0:s0 + win]
        Rw = (Rw - Rw.mean(0)) / Rw.std(0)
        Zw = fisher_z(np.corrcoef(Rw.T))
        w = np.sort(eigh(Zw)[0])
        feats.append(np.concatenate([w[-5:][::-1], w[:1], [np.mean(w ** 2)]]))
        _, V = top_vecs(Zw)
        st = all_stats(Zw, V=V, which=NS)
        targ.append([st[s] for s in NS])
    feats, targ = np.array(feats), np.array(targ)
    r2 = {}
    loo = LeaveOneOut()
    for j, s in enumerate(NS):
        pred = np.empty(len(feats))
        for tr, te in loo.split(feats):
            m = make_pipeline(StandardScaler(), LinearRegression())
            m.fit(feats[tr], targ[tr, j])
            pred[te] = m.predict(feats[te])
        ss = 1 - np.sum((targ[:, j] - pred) ** 2) / np.sum(
            (targ[:, j] - targ[:, j].mean()) ** 2)
        r2[s] = float(ss)
    res["sp500_r2"] = r2
    res["sp500_r2_cached"] = False
else:
    print("S&P csv absent; reusing cached results/exp4.json sp500_r2")
    r2 = json.load(open("../results/exp4.json"))["sp500_r2"]
    res["sp500_r2"] = r2
    res["sp500_r2_cached"] = True
print("S&P LOO R2:", r2)

# --------------------------------------------------------------------- fig
LABS = {"ipr": "IPR", "clust": "clust.", "mod": "modul.", "cv": "str. CV",
        "kurt": "kurt."}
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
for ax, key, ttl in [(axes[0], "cohortA",
                      "(a) rotatable cohort: signal in spikes"),
                     (axes[1], "cohortB",
                      "(b) non-rotatable cohort: signal in eigenvectors")]:
    mu = [res[key][k][0] for k in FEATSETS]
    se = [res[key][k][1] for k in FEATSETS]
    bars = ax.bar(range(len(FEATSETS)), mu, yerr=se, width=0.62,
                  color=[SERIES[0], SERIES[1], SERIES[3], SERIES[2]],
                  error_kw=dict(lw=0.9, ecolor=PAL["ink2"]))
    ax.axhline(0.5, color=PAL["muted"], lw=0.8, ls="--")
    ax.set_xticks(range(len(FEATSETS)),
                  ["spectral", "summaries", "residual", "both"], fontsize=7)
    ax.set_ylim(0.35, 1.05)
    ax.set_ylabel("CV AUC" if key == "cohortA" else "")
    ax.set_title(ttl, loc="left", fontsize=8.2)
    for b, m in zip(bars, mu):
        ax.annotate(f"{m:.2f}", (b.get_x() + b.get_width() / 2, m + 0.025),
                    ha="center", fontsize=6.6, color=PAL["ink2"])

ax = axes[2]
vals = [r2[s] for s in NS]
bars = ax.bar(range(len(NS)), vals, width=0.62, color=SERIES[0])
ax.set_xticks(range(len(NS)), [LABS[s] for s in NS], fontsize=7.5)
ax.set_ylim(0, 1.05)
ax.set_ylabel("LOO $R^2$ from eigenvalues")
ax.set_title("(c) S\\&P 500: summaries from eigenvalues", loc="left",
             fontsize=8.2)
for b, m in zip(bars, vals):
    ax.annotate(f"{m:.2f}", (b.get_x() + b.get_width() / 2, m + 0.02),
                ha="center", fontsize=7, color=PAL["ink2"])

fig.tight_layout()
fig.savefig("../figures/fig4_suff.pdf", bbox_inches="tight")
save_json(res, "../results/exp4_new.json")
print("exp4 done")
