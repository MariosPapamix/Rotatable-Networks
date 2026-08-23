# Notation ledger (symbol | meaning | first defined)

| Symbol | Meaning | First defined |
|---|---|---|
| Sym_n, Sym_inf | real symmetric n x n matrices; infinite symmetric arrays | Sec 1.3 Notation |
| O(n), Haar_n | orthogonal group; Haar probability on it | Sec 1.3 |
| X, X_n | random symmetric array; leading n x n corner | Sec 1.3 |
| hollow(A) | A with diagonal zeroed | Sec 1.3 |
| W, GOE(n) | infinite GOE pattern (N(0,1) off-diag, N(0,2) diag); its corner | Sec 1.3 |
| p_m(lambda) | power sum sum_k lambda_k^m, m >= 2 | Sec 1.3 |
| D(A) | eigenvalues in decreasing order | Sec 1.3 |
| V(A) | sign-normalized eigenvector matrix (simple spectrum) | Sec 1.3 |
| K_loc = 10 | eigenpairs by largest \|eigenvalue\| for localization statistics | Sec 1.3 |
| s_k = sqrt(n) lambda_k / sigma | effective spike at sample size n | Sec 1.3; used Prop S2 |
| eq (1) = eq:rotdef | joint rotatability U X U^T =d X | Sec 1 |
| tilde U = U + I | embedded orthogonal matrix | Sec 2.1 |
| Def 1 def:rot | joint rotatability (array form) | Sec 2.1 |
| Def 2 def:consistent | orthogonal invariance; sampling consistency | Sec 2.2 |
| R_{rho,sigma,lambda} | ergodic law, eq (2) = eq:rep | Thm 1 |
| P^net_{sigma,lambda} | ergodic rotatable network law (hollow part, rho=0) | Sec 3.1 |
| eq (3) = eq:cycles | cycle moments and second moment | Sec 3.2 |
| kappa(D, .) | Haar orbit law given spectrum | Thm 2 (i) |
| p_j, tilde p_j | upper-/lower-tail Monte Carlo p-values | Test box; Sec 5 |
| r_S, z_S | conditional residual; standardized residual | Def 3 def:residual |
| epsilon | TV gap to Haar conditional | Prop 1 prop:robust |
| T_IPR, T_KS | localization statistics | Sec 5.1 |
| Q_a(X) | quadratic form sum a_i a_j X_ij | Sec 6.1 |
| Def 4 def:prot | joint p-rotatability | Sec 6.1 |
| SaS(s) | symmetric p-stable, scale s | Sec 6.2 |
| Z^(p), c_p = 2^{(p-1)/p} | p-stable Wigner array; canonical diagonal scale | Sec 6.2; Thm 4 |
| eq (5) = eq:stablerep | spiked p-stable ensemble | Thm 4 |
| G(X) | normal-scores (Gaussianization) transform | Sec 6.3 |
| Supplement S-numbering | Lemma S1; Props S1-S3; Cor S1; Tables S1-S2; Figs S1-S2 | Supplement |

Checked: no symbol used before definition; no symbol defined twice with
different meanings; "network/node" fixed as the terminology ("graph"
reserved for graphons, binary graphs, and graph-theoretic objects F).

## Additions from the corrected theorem set (this revision)
- `\succeq_*`, `\ell_{2,*}`: canonical order on spike sequences; parameter space Theta = R x [0,inf) x l_{2,*}. Defined once in Section 2; makes the mixing measure exactly unique.
- `\nu_n`: sign-normalized Haar law of the eigenframe (Theorem 2(iii)). Not Haar measure on O(n); stated as such in main text.
- `\kappa(d,A)` with `\kappa_d := \kappa(d,.)`: Haar-orbit kernel; subscript form used in Proposition 1 and its proof.
- `\zeta^{(p)}_{ki}`: SaS(1) factors in Section 6, distinguishing them from the Gaussian factors `\xi_{ki}` of Theorem 1.
- `\varsigma(a)`: sign factor local to the Theorem 1 proof (Supplement).
- `m_\theta`: directing finite measure in the Theorem 1 proof (renamed from the source's nu_theta to avoid collision with nu_n).
- `\mathcal L_{j,d}`: law of T_j under the Haar orbit at spectrum d (Theorem 3 proof; renamed from nu_{j,d}).
- `s_*(a)`: stable scale local to the Theorem 4 proof (renamed from s(a) to avoid the eigenframe sign map s(v)).
- `s_o, s_d`: generic off-diagonal/diagonal stable scales, local to Theorem 4(ii).
- `osc(g)`, `d_TV`: defined once in Section 2.
- `\Sigma`: generic sigma-field in the total-variation convention (Supplement; the source's mathscr S).
- `\mu_i` in Proposition S2 remains the BBP eigenvalue notation, unrelated to any measure.
- `\rho`: diagonal ergodic parameter (Theorem 1); unobserved for networks.
- `Q_d(eps)`: concentration modulus of the replicate-statistic law (Proposition 2); rational-endpoint supremum, Borel in d.
- `Delta = max_i |X_ii|`: unobserved diagonal magnitude of the completed matrix (Proposition 2).
- `ell_i`: energy of node i in the K=10 leading eigenvectors; `p_i^adj`: single-step max-adjusted nodewise p-value (Section 5, exp7).
- `q_*`: shorthand E Q_{D(X)}(2 L Delta) inside the Proposition 2 proof (supplement-local).
