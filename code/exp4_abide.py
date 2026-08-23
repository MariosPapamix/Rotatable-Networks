"""[TO RUN] Real-cohort replication of experiment 4 on ABIDE.

This script was NOT executed in the environment that produced the paper
(no network access to the ABIDE bucket); every number it would produce is
therefore absent from the manuscript, and Section 7.4 uses the labeled
synthetic cohorts instead. Running it requires:

    pip install nilearn
    python fetch_real_data.py abide      # downloads ~1 GB via nilearn
    python exp4_abide.py

It then applies the identical pipeline to real subjects: Fisher-z
per-subject functional connectivity (CC200 parcellation), four feature
sets (spectral, non-spectral summaries, conditional residuals, spectral
plus summaries), stratified 5-fold logistic AUC for DX_GROUP, and
per-subject Haar rotatability tests, writing
../results/exp4_abide.json and ../figures/fig4_abide.pdf. If adopted,
these results would replace the synthetic Section 7.4 after the text is
updated accordingly.
"""
import json, os, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.linalg import eigh
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from common import (mpl_style, all_stats, top_vecs, fisher_z, haar,
                    conditional_pvalues, PAL, SERIES, save_json)

mpl_style()
FC = "../data/abide_fc/fc_cc200.npy"          # (S, n, n) correlations
DX = "../data/abide_fc/dx_group.npy"          # (S,) 1 = ASD, 2 = control
if not (os.path.exists(FC) and os.path.exists(DX)):
    sys.exit("ABIDE arrays not found. Run: python fetch_real_data.py abide")

NS = ["ipr", "clust", "mod", "cv", "kurt"]
BRES, BTEST = 32, 99
rng = np.random.default_rng(20260812)

def spectral_features(X, k=10):
    w = np.sort(eigh(X)[0])
    return np.concatenate([w[-k:][::-1], w[:3],
                           [np.mean(w**2), np.mean(w**3), np.mean(w**4)]])

def nonspectral_features(X):
    _, V = top_vecs(X)
    st = all_stats(X, V=V, which=NS)
    return np.array([st[s] for s in NS])

mats = fisher_z(np.load(FC))
y = (np.load(DX) == 1).astype(int)
S, n, _ = mats.shape
print(f"ABIDE: {S} subjects, n={n}")

Xs = np.array([spectral_features(M) for M in mats])
Xn = np.array([nonspectral_features(M) for M in mats])
Xr = np.empty_like(Xn)
for i, M in enumerate(mats):
    w = eigh(M)[0]
    nulls = np.empty((BRES, len(NS)))
    for b in range(BRES):
        U = haar(n, rng)
        nulls[b] = nonspectral_features((U * w) @ U.T)
    Xr[i] = Xn[i] - nulls.mean(0)

def cv_auc(F):
    skf = StratifiedKFold(5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(F, y):
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000))
        clf.fit(F[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(F[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs) / np.sqrt(5))

res = {"n_subjects": int(S), "n": int(n),
       "spectral": cv_auc(Xs), "non-spectral": cv_auc(Xn),
       "residual": cv_auc(Xr), "both": cv_auc(np.hstack([Xs, Xn]))}

sel = rng.choice(S, size=min(40, S), replace=False)
rej = {0: [], 1: []}
for i in sel:
    pv, _ = conditional_pvalues(mats[i], BTEST, rng, which=NS, hollow=True)
    rej[int(y[i])].append(min(pv.values()) <= 0.05 / len(NS))
res["per_subject_rej"] = {g: float(np.mean(v)) for g, v in rej.items() if v}
save_json(res, "../results/exp4_abide.json")
print(json.dumps(res, indent=2))

fig, ax = plt.subplots(figsize=(4.2, 2.6))
keys = ["spectral", "non-spectral", "residual", "both"]
mu = [res[k][0] for k in keys]; se = [res[k][1] for k in keys]
ax.bar(range(4), mu, yerr=se, width=0.62,
       color=[SERIES[0], SERIES[1], SERIES[3], SERIES[2]],
       error_kw=dict(lw=0.9, ecolor=PAL["ink2"]))
ax.axhline(0.5, color=PAL["muted"], lw=0.8, ls="--")
ax.set_xticks(range(4), ["spectral", "summaries", "residual", "both"],
              fontsize=7)
ax.set_ylabel("CV AUC (DX)")
ax.set_title("ABIDE per-subject cohort (real data)", loc="left", fontsize=8.5)
fig.tight_layout()
fig.savefig("../figures/fig4_abide.pdf", bbox_inches="tight")
print("wrote figures/fig4_abide.pdf")
