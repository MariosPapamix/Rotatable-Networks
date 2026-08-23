"""Appendix pilot: combined-test power against a strong smooth graphon
(spikes 5.0, 2.5; n=200) -- quantifies 'weak smooth graphons are
effectively rotatable' (Section 7.2)."""
import json, numpy as np
from common import goe, conditional_pvalues
rng = np.random.default_rng(7)
n, c, R = 200, 10.0, 50
rej = 0
for r in range(R):
    U = rng.uniform(0, 1, n)
    F = np.stack([np.cos(np.pi*U), np.cos(2*np.pi*U)], axis=1)
    X = goe(n, rng)/np.sqrt(n) + c*(F @ np.diag([1.0, 0.5]) @ F.T)/n
    pv, _ = conditional_pvalues(X, 99, rng,
                                which=["ipr","evks","clust","cv","kurt"])
    rej += min(pv.values()) <= 0.01
json.dump({"n": n, "spikes": [5.0, 2.5], "reps": R, "power": rej/R},
          open("../results/exp2b.json", "w"))
print("strong-graphon power:", rej/R)
