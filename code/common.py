"""Shared utilities for the rotatable-networks paper experiments.

Model conventions
-----------------
Finite-n spiked ensemble (the ergodic rotatable law, Theorem 1):
    X = (sigma/sqrt(n)) * W + sum_k theta_k u_k u_k',
with W a GOE-type Wigner matrix (W_ij = W_ji ~ N(0,1) off-diagonal,
W_ii ~ N(0,2)), u_k orthonormal (Haar), theta_k spike strengths.
Bulk edge at 2*sigma; BBP: lambda_1 -> theta + sigma^2/theta iff theta > sigma.

Infinite-array (Kallenberg) coordinates: X_ij = sigma*zeta_ij
+ sum_k lam_k (xi_ki xi_kj - delta_ij); the n x n corner, divided by sqrt(n),
is the ensemble above with theta_k = sqrt(n) * lam_k  (Proposition BBP).
"""

import json
import numpy as np
from numpy.linalg import eigh, qr

# ---------------------------------------------------------------- palette ---
PAL = {
    "blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a",
    "yellow": "#eda100", "magenta": "#e87ba4", "green": "#008300",
    "violet": "#4a3aa7", "red": "#e34948",
    "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781",
    "grid": "#e1e0d9", "axis": "#c3c2b7",
}
SERIES = [PAL["blue"], PAL["orange"], PAL["aqua"], PAL["yellow"],
          PAL["magenta"], PAL["green"], PAL["violet"], PAL["red"]]


def mpl_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": PAL["axis"], "axes.linewidth": 0.8,
        "xtick.color": PAL["muted"], "ytick.color": PAL["muted"],
        "axes.labelcolor": PAL["ink2"], "text.color": PAL["ink"],
        "axes.grid": True, "grid.color": PAL["grid"], "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "lines.linewidth": 1.6, "lines.markersize": 4.0,
        "axes.prop_cycle": mpl.cycler(color=SERIES),
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


# ---------------------------------------------------------------- sampling ---
def goe(n, rng):
    """GOE-type symmetric matrix: off-diag N(0,1), diag N(0,2)."""
    A = rng.standard_normal((n, n))
    W = (A + A.T) / np.sqrt(2.0)
    return W


def haar(n, rng):
    """Haar-distributed orthogonal matrix (QR with sign correction)."""
    Z = rng.standard_normal((n, n))
    Q, R = qr(Z)
    return Q * np.sign(np.diag(R))


def spiked_goe(n, thetas, sigma, rng, vecs=None):
    """X = (sigma/sqrt n) W + sum theta_k u_k u_k^T, u Haar orthonormal."""
    X = (sigma / np.sqrt(n)) * goe(n, rng)
    if len(thetas):
        if vecs is None:
            U = haar(n, rng)[:, :len(thetas)]
        else:
            U = vecs
        X = X + (U * np.asarray(thetas)) @ U.T
    return X


def sas(alpha, size, rng):
    """Symmetric alpha-stable via Chambers-Mallows-Stuck (scale 1)."""
    U = rng.uniform(-np.pi / 2, np.pi / 2, size=size)
    Wexp = rng.exponential(1.0, size=size)
    if abs(alpha - 1.0) < 1e-12:
        return np.tan(U)
    return (np.sin(alpha * U) / np.cos(U) ** (1 / alpha)
            * (np.cos((1 - alpha) * U) / Wexp) ** ((1 - alpha) / alpha))


def stable_wigner(n, alpha, rng):
    A = sas(alpha, (n, n), rng)
    return np.triu(A, 1) + np.triu(A, 1).T + np.diag(np.diag(A))


# ------------------------------------------------------------- statistics ---
def offdiag(X):
    n = X.shape[0]
    return X[~np.eye(n, dtype=bool)]


def ipr_stat(V, k=10):
    """Max inverse participation ratio over the k leading eigenvectors
    (by |eigenvalue|); Haar null ~ 3/n. Larger = localized."""
    return float((V ** 4).sum(axis=0).max())


def top_vecs(X, k=10):
    w, V = eigh(X)
    idx = np.argsort(-np.abs(w))[:k]
    return w, V[:, idx]


def onnela_clustering(X):
    """Mean Onnela weighted clustering on the positive part, hollow."""
    A = np.maximum(X, 0.0)
    np.fill_diagonal(A, 0.0)
    m = A.max()
    if m <= 0:
        return 0.0
    W3 = np.cbrt(A / m)
    num = np.einsum('ii->i', W3 @ W3 @ W3)
    deg = (A > 0).sum(axis=1).astype(float)
    den = deg * (deg - 1)
    ok = den > 0
    return float(np.mean(num[ok] / den[ok])) if ok.any() else 0.0


def spectral_modularity(X):
    """Modularity value of the leading-eigenvector bipartition (Newman),
    computed on the positive part."""
    A = np.maximum(X, 0.0)
    np.fill_diagonal(A, 0.0)
    k = A.sum(axis=1)
    two_m = k.sum()
    if two_m <= 0:
        return 0.0
    B = A - np.outer(k, k) / two_m
    w, V = eigh(B)
    s = np.sign(V[:, -1])
    s[s == 0] = 1.0
    return float(s @ B @ s / two_m)


def strength_cv(X):
    s = np.abs(X).sum(axis=1)
    return float(np.std(s) / np.mean(s))


def excess_kurtosis(X):
    x = offdiag(X)
    x = x - x.mean()
    m2 = np.mean(x ** 2)
    return float(np.mean(x ** 4) / m2 ** 2 - 3.0)


def evec_ks(V, k=3):
    """Max Kolmogorov-Smirnov distance of sqrt(n) * (leading eigenvector
    entries) from N(0,1); Haar null ~ delocalized Gaussian entries.
    Detects non-Gaussian eigenvector geometry (graphon eigenfunctions).
    Computed directly (identical value to scipy.stats.kstest's two-sided
    statistic, vectorized for speed)."""
    from scipy.stats import norm
    n = V.shape[0]
    grid = np.arange(1, n + 1) / n
    best = 0.0
    for j in range(min(k, V.shape[1])):
        c = norm.cdf(np.sort(np.sqrt(n) * V[:, j]))
        d = max(np.max(grid - c), np.max(c - (grid - 1.0 / n)))
        best = max(best, float(d))
    return best


STAT_NAMES = ["ipr", "evks", "clust", "mod", "cv", "kurt"]


def all_stats(X, V=None, which=STAT_NAMES):
    if V is None and ("ipr" in which or "evks" in which):
        _, V = top_vecs(X)
    out = {}
    if "ipr" in which:
        out["ipr"] = ipr_stat(V)
    if "evks" in which:
        out["evks"] = evec_ks(V)
    if "clust" in which:
        out["clust"] = onnela_clustering(X)
    if "mod" in which:
        out["mod"] = spectral_modularity(X)
    if "cv" in which:
        out["cv"] = strength_cv(X)
    if "kurt" in which:
        out["kurt"] = excess_kurtosis(X)
    return out


# ------------------------------------------------- exact conditional test ---
def haar_null_stats(X, B, rng, which=STAT_NAMES, return_null=False,
                    hollow=False):
    """Draw B Haar-conditional replicates X' = U diag(spec(X)) U' and
    compute statistics.  With hollow=True every statistic is computed on
    the diagonal-zeroed matrix (matching hollow observed networks); this
    pipeline is approximate (see Remark 3 of the paper: exactness holds
    for fully observed matrices).  With hollow=False, eigenvectors of X'
    are the columns of
    U, so the IPR statistic needs no re-diagonalisation."""
    n = X.shape[0]
    w = eigh(X)[0]
    idx = np.argsort(-np.abs(w))[:10]
    null = {s: np.empty(B) for s in which}
    for b in range(B):
        U = haar(n, rng)
        Xb = (U * w) @ U.T
        if hollow:
            np.fill_diagonal(Xb, 0.0)
            Vb = top_vecs(Xb)[1] if ("ipr" in which or "evks" in which) else None
        else:
            Vb = U[:, idx]
        sb = all_stats(Xb, V=Vb, which=which)
        for s in which:
            null[s][b] = sb[s]
    if return_null:
        return null
    return null


def conditional_pvalues(X, B, rng, which=STAT_NAMES, return_null=False,
                        hollow=False):
    """Exact Monte-Carlo conditional p-values (upper tail):
    p = (1 + #{b : T_b >= T_obs}) / (B + 1)."""
    X0 = X.copy()
    if hollow:
        np.fill_diagonal(X0, 0.0)
    _, V = top_vecs(X0)
    obs = all_stats(X0, V=V, which=which)
    null = haar_null_stats(X, B, rng, which=which, return_null=True,
                           hollow=hollow)
    pvals = {s: float((1 + np.sum(null[s] >= obs[s] - 1e-12)) / (B + 1))
             for s in which}
    if return_null:
        return pvals, obs, null
    return pvals, obs


# ----------------------------------------------------------------- helpers ---
def fisher_z(C, clip=0.999):
    C = np.clip(C, -clip, clip)
    Z = np.arctanh(C)
    np.fill_diagonal(Z, 0.0)
    return Z


def hill_estimator(x, k=None):
    """Hill tail-index estimate from the top-k order statistics of |x|."""
    x = np.sort(np.abs(np.asarray(x, dtype=float)))
    x = x[x > 0]
    if k is None:
        k = max(20, int(0.05 * len(x)))
    tail = x[-k - 1:]
    logs = np.log(tail[1:]) - np.log(tail[0])
    return float(1.0 / np.mean(logs))


def gaussianize(X, rng=None):
    """Normal-scores (rank) transform of the off-diagonal entries,
    preserving symmetry and the zero diagonal.  Ties are broken uniformly
    at random when an rng is supplied (deterministically, by input order,
    otherwise).  With many tied entries -- e.g. structural zeros in sparse
    flow networks -- conclusions should be checked across tie-breaking
    seeds; exp3/exp5 do this."""
    from scipy.stats import norm
    n = X.shape[0]
    iu = np.triu_indices(n, 1)
    v = X[iu]
    if rng is None:
        order = np.argsort(v, kind="stable")
    else:
        order = np.lexsort((rng.random(v.size), v))
    r = np.empty(v.size, dtype=float)
    r[order] = np.arange(v.size)
    q = norm.ppf((r + 0.5) / v.size) * np.std(v)
    Y = np.zeros_like(X)
    Y[iu] = q
    return Y + Y.T


def spiked_stable(n, lams, sigma, alpha, rng, return_factors=False):
    """Theorem 4 construction (the p-rotatable ergodic law), p = alpha:
        X_ij = sigma * Z_ij + sum_k lam_k * xi_ki * xi_kj   (i != j),
    with Z the SaS(alpha) Wigner array with the canonical diagonal
    Z_ii ~ SaS(2^{(alpha-1)/alpha}) forced by Theorem 4(ii), xi_ki iid
    SaS(alpha) factors, and the diagonal completed as
    X_ii = sigma * Z_ii + sum_k lam_k * xi_ki^2."""
    A = sas(alpha, (n, n), rng)
    Z = np.triu(A, 1) + np.triu(A, 1).T
    np.fill_diagonal(Z, 2.0 ** ((alpha - 1.0) / alpha) * sas(alpha, n, rng))
    Xi = (sas(alpha, (len(lams), n), rng) if len(lams)
          else np.zeros((0, n)))
    X = sigma * Z
    for k in range(len(lams)):
        X = X + lams[k] * np.outer(Xi[k], Xi[k])
    return (X, Xi) if return_factors else X


def save_json(obj, path):
    def enc(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=enc)

def minp_pvalue(obs, null):
    """Exact single-step min-p combination using the same Haar draws.

    obs: dict statistic -> observed value; null: dict -> array of B null
    values.  Returns the Monte Carlo p-value of the combined statistic
    q = min_j p_j, calibrated jointly (exact by exchangeability of the
    B+1 joint statistic vectors; cf. Hemerik & Goeman 2018).
    """
    import numpy as _np
    stats = list(obs)
    B = len(null[stats[0]])
    allv = {s: _np.concatenate([[obs[s]], null[s]]) for s in stats}
    # p-value of each of the B+1 draws for each statistic, within the pool
    q = _np.ones(B + 1)
    for s in stats:
        v = allv[s]
        ge = (v[None, :] >= v[:, None] - 1e-12).sum(1) / (B + 1.0)
        q = _np.minimum(q, ge)
    return float((q <= q[0] + 1e-12).sum() / (B + 1.0))


def rotatable_test(X, B=499, seed=0, which=STAT_NAMES, hollow=True):
    """One-call interface: Haar-conditional rotation test + residuals.

    Returns a dict with per-statistic observed value, Haar-conditional
    mean and sd, standardized residual z, upper-tail p, and the joint
    min-p combination.
    """
    import numpy as _np
    rng = _np.random.default_rng(seed)
    pv, obs, null = conditional_pvalues(X, B, rng, which=which,
                                        return_null=True, hollow=hollow)
    out = {}
    for s in which:
        mu, sd = float(_np.mean(null[s])), float(_np.std(null[s]))
        out[s] = {"observed": float(obs[s]), "null_mean": mu,
                  "null_sd": sd,
                  "z": float((obs[s] - mu) / sd) if sd > 0 else float("nan"),
                  "p_upper": pv[s]}
    out["minp_combined"] = minp_pvalue(obs, null)
    return out
