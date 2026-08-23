"""Experiment 2c: empirical size of the hollow test on sample-correlation
networks built from isotropic Gaussian panels (finite-T pipeline check).

Y is a T x n panel of iid N(0,1) entries; the analyzed network is the
Fisher-z transformed sample correlation matrix with zeroed diagonal,
exactly the S&P 500 pipeline of exp3/exp4.  The unit-diagonal constraint
means this matrix is not the hollow part of an exactly rotatable array,
so Theorem 3 does not apply verbatim; this experiment measures the
resulting size distortion of the hollow test at nominal 0.05.
"""
import sys, json, os
import numpy as np
from common import conditional_pvalues, fisher_z, save_json

n, B, R = 120, 99, 100
WHICH = ["ipr", "evks", "clust", "cv", "kurt"]
PART = "../results/exp2c_parts"
os.makedirs(PART, exist_ok=True)
MODE = sys.argv[1] if len(sys.argv) > 1 else "assemble"
if MODE == "all":
    import subprocess
    for _st in ['120', '360', 'assemble']:
        subprocess.run([sys.executable, __file__, _st], check=True)
    sys.exit(0)
res = {}
Ts = [int(MODE)] if MODE != "assemble" else []
for T in Ts:
    rng = np.random.default_rng(20260812 + T)
    pmat = {s: np.empty(R) for s in WHICH}
    for r in range(R):
        Y = rng.standard_normal((T, n))
        Z = fisher_z(np.corrcoef(Y.T))
        pv, _ = conditional_pvalues(Z, B, rng, which=WHICH, hollow=True)
        for s in WHICH:
            pmat[s][r] = pv[s]
    comb = np.mean(np.minimum.reduce([pmat[s] for s in WHICH])
                   <= 0.05 / len(WHICH))
    out = {"per_stat_rej05": {s: float(np.mean(pmat[s] <= 0.05))
                              for s in WHICH},
           "combined_rej05": float(comb),
           "p": {s: pmat[s] for s in WHICH}}
    save_json(out, f"{PART}/T{T}.json")
    print(f"T={T}: combined rej @.05 = {comb:.3f}; per-stat:",
          {s: round(float(np.mean(pmat[s] <= 0.05)), 3) for s in WHICH})
if MODE == "assemble":
    for T in (120, 360):
        res[f"T{T}"] = json.load(open(f"{PART}/T{T}.json"))
    save_json(res, "../results/exp2c.json")
    print("exp2c assembled")
