"""exp7: nodewise attribution with exact familywise control.

Per-node statistic: energy of node i in the K=10 leading eigenvectors
(by |eigenvalue|) of the hollow matrix, l_i = sum_k v_k(i)^2.  Single-step
max adjustment: p_i = (1 + #{b: max_j l_j(X^(b)) >= l_i(A)}) / (B+1).
Under rotatability, P(any node flagged at level a) <= a exactly, because
min_i p_i is the Theorem-3 p-value of the statistic max_i l_i.

Outputs results/exp7_nodewise.json.  Seeds fixed; ~1 minute.
"""
import json
import numpy as np
import pandas as pd
from common import fisher_z

K = 10
B = 499
ALPHA = 0.05
SEED = 20260870


def node_energy(M, K=K):
    w, V = np.linalg.eigh(M)
    idx = np.argsort(-np.abs(w))[:K]
    return (V[:, idx] ** 2).sum(axis=1)


def nodewise(A, B, rng):
    A = A.copy()
    np.fill_diagonal(A, 0.0)
    n = A.shape[0]
    d = np.linalg.eigvalsh(A)
    obs = node_energy(A)
    maxes = np.empty(B)
    for b in range(B):
        Z = rng.standard_normal((n, n))
        Q, _ = np.linalg.qr(Z)
        Xb = (Q * d) @ Q.T
        np.fill_diagonal(Xb, 0.0)
        maxes[b] = node_energy(Xb).max()
    padj = (1 + (maxes[None, :] >= obs[:, None]).sum(1)) / (B + 1)
    return obs, padj


res = {}
rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------- HCP
FC = np.loadtxt("../data/funcMatrix_ctx_schaefer_200.csv", delimiter=",")
Zfc = fisher_z(FC)          # match the Table-1 pipeline
obs, padj = nodewise(Zfc, B, rng)
flag = np.where(padj <= ALPHA)[0]
order = np.argsort(padj)
res["hcp"] = {
    "n": int(Zfc.shape[0]), "B": B, "alpha": ALPHA,
    "n_flagged": int(len(flag)),
    "global_p": float(padj.min()),
    "top_nodes": [[int(i), float(obs[i]), float(padj[i])] for i in order[:8]],
}
print("HCP flagged:", len(flag), "of", FC.shape[0],
      "| global p:", padj.min())

# --------------------------------------------------------- air traffic
fl = pd.read_csv("../data/flights_airport.csv")
airports = pd.unique(pd.concat([fl.origin, fl.destination]))
amap = {a: i for i, a in enumerate(airports)}
nA = len(airports)
A = np.zeros((nA, nA))
for o, dd, c in zip(fl.origin, fl.destination, fl["count"]):
    A[amap[o], amap[dd]] += c
A = A + A.T
deg = (A > 0).sum(1)
keep = deg >= 5
A = A[np.ix_(keep, keep)]
names = np.array(airports)[keep]
np.fill_diagonal(A, 0.0)
Alog = np.log1p(A)

obs, padj = nodewise(Alog, B, rng)
flag = np.where(padj <= ALPHA)[0]
order = np.argsort(-obs)
res["flights"] = {
    "n": int(Alog.shape[0]), "B": B, "alpha": ALPHA,
    "n_flagged": int(len(flag)),
    "global_p": float(padj.min()),
    "flagged_airports": [[str(names[i]), float(obs[i]), float(padj[i])]
                         for i in order if padj[i] <= ALPHA][:12],
}
print("flights flagged:", len(flag), "of", Alog.shape[0],
      "| global p:", padj.min())
print("top flagged:", [x[0] for x in res["flights"]["flagged_airports"]][:8])

with open("../results/exp7_nodewise.json", "w") as f:
    json.dump(res, f, indent=1)
print("saved results/exp7_nodewise.json")
