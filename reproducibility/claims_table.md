# Claims table (claim | type | evidence | status)

| Claim | Type | Evidence | Status |
|---|---|---|---|
| Rotatable + sampling-consistent = mixtures of spiked GOE | theory | Thm 1 (representation incl. symmetric case due to Kallenberg; model-class equivalence via Lemma S1 is ours); proofs S2 | proved |
| Exactness scope: full matrices exact; hollow/unit-diagonal degenerate, pipeline approximate | theory + empirical | Remark 3 (degeneracy argument); Prop 1; measured sizes 0.068 / 0.060 / 0.040 | proved + measured |
| Competing nulls miscalibrate under rotatable truth (shuffle 0.51; fitted bootstrap 0.03, per-stat to 0.10) | empirical | Sec 7.2; Table S2; exp6_baselines.py | run |
| Cycle moments are power sums; identifiability | theory | Prop S1 (corrected scope: vertex-disjoint cycles; bridges vanish; counterexamples stated) | proved |
| Spectrum sufficient; Haar conditional; converse | theory | Thm 2 | proved |
| Non-spectral summaries conditionally parameter-free; Rao-Blackwell; frame ancillarity | theory | Cor S1 | proved |
| Exact finite-sample composite test | theory | Thm 3 | proved |
| Size <= alpha + epsilon under approximate rotatability; bounded-loss risk bound | theory | Prop 1 | proved |
| IPR consistency, supercritical l4-nondegenerate spikes; null O_P(log n / n) | theory | Prop S3 | proved |
| BBP calibration; n^{-1/2} boundary; plug-in consistency | theory (transfer) | Prop S2, citing Benaych-Georges & Nadakuditi; Onatski; Perry et al. | proved (transfer) |
| Spiked p-stable ensembles p-rotatable; forced diagonal 2^{(p-1)/p} | theory | Thm 4 | proved |
| Completeness of p-stable construction | conjecture | Conj 1 | open (stated as such) |
| Test level correct in hollow pipeline | empirical | Fig 2a; 0.068 (Wilson [0.043,0.106]) | run |
| Correlation-pipeline size approx nominal | empirical | Table S1; 0.060 / 0.040 | run |
| Power: localization 0.82, hetero 1.0; graphon blind spot 0.02-0.08, strong 0.56 | empirical | Fig 2b; exp2b | run |
| HCP directions (IPR+, kurt+, mod-, cv-, clust null; G readout) | empirical | Table 1; Fig 3 | run |
| S&P shares and lower-tail findings | empirical | Table 1 (cached; rerun step 1) | run (cache) |
| Air-traffic persistence under G across tie seeds | empirical | Table 1; Sec 7.3 | run |
| Cohort A/B AUC pattern; residuals at chance / at 1.00 | empirical (synthetic, labeled) | Fig 4a,b | run |
| S&P LOO R^2 per summary | empirical (descriptive) | Fig 4c (cached; rerun step 1) | run (cache) |
| Heavy-tail Hill indices; diagnostic cases | empirical | Fig S1 | run |
| "To our knowledge, new" reduction + network use | novelty claim | Related work para; literature check | qualified, cited |

Exaggeration audit: no "first/only/always/guarantees" outside proved
statements; "exact" used only for Thm 3-type guarantees; universal
quantifiers restricted to theorem scopes.

## Updates from the corrected theorem set (this revision)
| Claim | Status |
|---|---|
| Theorem 1: mixing measure exactly unique on the canonically ordered parameter space; R_theta exactly the ergodic (extreme) laws; theta -> R_theta a Borel kernel; a.s. convergence simultaneous over entries; finite-family equivalence internal to the theorem | Proved in Supplement, self-contained from Kallenberg (1988, Thm 5.1, eqs. (17)-(18)); Lemma S1 retired |
| Theorem 2: common-kernel sufficiency without domination; measurable eigenframe V constructed explicitly; V ~ nu_n (sign-normalized, not Haar); converse (iv) | Proved in Supplement |
| Theorem 3: floor bound floor((B+1)u)/(B+1); tie robustness; conditional uniformity under atomless orbit laws; conditional Bonferroni | Proved in Supplement (weak-upper-rank argument) |
| Proposition 1: preprocessing clause (T o h with D from the unprocessed matrix); floor-refined size bounds; osc-loss risk bounds | Proved in Supplement (product-TV, single charge of epsilon) |
| Hollow/unit-diagonal conditional laws vs Haar orbits: mutually singular, TV distance = 1; Prop 1 vacuous for those pipelines | Supplement Remarks S1-S2; main Remark 3, intro, and discussion reworded accordingly; empirical sizes remain the calibration evidence |
| Theorem 4: q in (0,p/2); joint disjoint-support quadratic forms; conditional row CF; p=2 endpoint corrected: centered variant = full ergodic class (sigma_rep = sqrt2 tau, lam_rep = 2 gamma), uncentered construction = strict subfamily (rho_rep = sum lam_rep, lam in union of l_q, q<1) | Proved in Supplement; prior 'coincides at p=2' claim replaced |
| Proposition 2 (statistic-level size of the zero-completion pipeline): size <= u + 2 sqrt(E Q_{D(X)}(2 L Delta)) + 1/(B+1); rate n^{-1/4} (log n)^{1/2} under Lipschitz + bounded density | Proved in Supplement (coupling to the exact comparator via Weyl + window count + Markov); E max|X_ii| = O(log n) recorded in Section S2 |
| Corollary S2 (corrected to a trichotomy): spherical-direction spikes are orthogonally invariant, hence inside the composite null with exact level at every strength; localized (l4-nondegenerate) spikes: consistent detection at every fixed s>1; below s=1 no consistent test (contiguity, Perry et al. 2018) | Proved in Supplement; earlier boundary-attainment phrasing retracted |
| Nodewise attribution, exact FWER via max statistic; air traffic flags the 8 major hubs (global p = 0.006), HCP (Fisher-z, matching Table 1) flags 8 of 200 parcels, global p = 0.008; raw-FC variant flags none, preprocessing-alignment fix recorded | exp7_nodewise.json, seeds fixed |

| Proposition S2 restated: sigma>0; ratio consistency mu_k/s_k -> 1 with O_P(1) location gap from the realized factor norm; outlier count over [-2-eps,2+eps]; eps_n rate n^{2/3}eps_n -> infinity; mirrored negative-spike estimator; hollow robustness clause | Statement corrected per referee audit; proof addenda in supplement |
| Lemma S1 (E max diag = O(log n)) promoted from heuristic text | Proved in Supplement |
