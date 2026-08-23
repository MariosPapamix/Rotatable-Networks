"""Experiment 6: calibration of competing null models under rotatable truth.

Under the SAME null ensemble as the level study of experiment 2 (spiked
GOE, 0-3 random spikes theta ~ U(1.2, 3.0), sigma=1, n=120, hollow
pipeline, B=99, Bonferroni-combined battery at nominal 0.05), measure the
rejection rate of two null models used in practice:

  shuffle : weight-permutation null (off-diagonal weights of the observed
            hollow matrix are randomly permuted, symmetrically);
  pboot   : parametric bootstrap from the fitted spiked model
            (sigma-hat from the bulk quartile, lambda-hat from
            psi^{-1} of outliers, replicates simulated from the fit).

The Haar-conditional test's rate on this ensemble is 17/250 = 0.068
(results/exp2_test.json). Stages: shuffle0 shuffle1 pboot0 pboot1
assemble | all.
"""
import json, os, sys
import numpy as np
from numpy.linalg import eigh
from common import (spiked_goe, haar, all_stats, top_vecs, save_json,
                    STAT_NAMES)

STAGE = sys.argv[1] if len(sys.argv) > 1 else "assemble"
if STAGE == "all":
    import subprocess
    for _st in ["shuffle0", "shuffle1", "pboot0", "pboot1", "assemble"]:
        subprocess.run([sys.executable, __file__, _st], check=True)
    sys.exit(0)

WHICH = ["ipr", "evks", "clust", "cv", "kurt"]
n, B, R = 120, 99, 100
ALPHA = 0.05
PART = "../results/exp6_parts"
os.makedirs(PART, exist_ok=True)
Q75_SC = 0.8079  # upper quartile of the radius-2 semicircle

def observed_stats(X0):
    _, V = top_vecs(X0)
    return all_stats(X0, V=V, which=WHICH)

def pvals_against(null, obs):
    return {s: float((1 + np.sum(null[s] >= obs[s] - 1e-12)) / (B + 1))
            for s in WHICH}

def shuffle_null(X0, rng):
    iu = np.triu_indices(n, 1)
    vals = X0[iu]
    null = {s: np.empty(B) for s in WHICH}
    for b in range(B):
        Xb = np.zeros_like(X0)
        Xb[iu] = rng.permutation(vals)
        Xb = Xb + Xb.T
        _, Vb = top_vecs(Xb)
        sb = all_stats(Xb, V=Vb, which=WHICH)
        for s in WHICH:
            null[s][b] = sb[s]
    return null

def fit_spiked(X0):
    """sigma-hat from the bulk quartile (semicircle radius 2*sigma),
    theta-hat = sigma*s from outliers mu = w/sigma with |mu| > 2.1,
    s = (|mu| + sqrt(mu^2 - 4))/2, matching Proposition S2 at the
    order-1 scale of spiked_goe."""
    w = np.sort(eigh(X0)[0])
    sig = np.quantile(w, 0.75) / Q75_SC
    sig = max(sig, 1e-3)
    mu = w / sig
    th = []
    for m in mu:
        if abs(m) > 2.1:
            s = (abs(m) + np.sqrt(m * m - 4)) / 2
            th.append(np.sign(m) * sig * s)
    return sig, sorted(th, key=abs, reverse=True)

def pboot_null(X0, rng):
    sig, th = fit_spiked(X0)
    null = {s: np.empty(B) for s in WHICH}
    for b in range(B):
        Xb = spiked_goe(n, th, sig, rng)
        np.fill_diagonal(Xb, 0.0)
        _, Vb = top_vecs(Xb)
        sb = all_stats(Xb, V=Vb, which=WHICH)
        for s in WHICH:
            null[s][b] = sb[s]
    return null

if STAGE.startswith(("shuffle", "pboot")):
    arm = "shuffle" if STAGE.startswith("shuffle") else "pboot"
    ci = int(STAGE[-1]); m = R // 2
    rng = np.random.default_rng((20260850 if arm == "shuffle" else 20260860)
                                + 10 * ci)
    pmat = {s: np.empty(m) for s in WHICH}
    for r in range(m):
        k = rng.integers(0, 4)
        th = np.sort(rng.uniform(1.2, 3.0, size=k))[::-1] if k else []
        X = spiked_goe(n, th, 1.0, rng)
        np.fill_diagonal(X, 0.0)
        obs = observed_stats(X)
        null = shuffle_null(X, rng) if arm == "shuffle" else pboot_null(X, rng)
        pv = pvals_against(null, obs)
        for s in WHICH:
            pmat[s][r] = pv[s]
    save_json(pmat, f"{PART}/{arm}_c{ci}.json")
    print(f"{arm} chunk {ci} done")
    sys.exit(0)

if STAGE == "assemble":
    out = {"n": n, "B": B, "R": R, "alpha": ALPHA,
           "haar_reference": {"combined": 17 / 250,
                              "note": "from exp2 level study, same ensemble"}}
    for arm in ["shuffle", "pboot"]:
        ps = {s: np.concatenate([
            np.asarray(json.load(open(f"{PART}/{arm}_c{c}.json"))[s])
            for c in (0, 1)]) for s in WHICH}
        comb = np.min(np.column_stack([ps[s] for s in WHICH]), 1) <= ALPHA / len(WHICH)
        out[arm] = {"combined": float(np.mean(comb)),
                    "per_stat": {s: float(np.mean(ps[s] <= ALPHA))
                                 for s in WHICH}}
    save_json(out, "../results/exp6_baselines.json")
    print(json.dumps(out, indent=2))
