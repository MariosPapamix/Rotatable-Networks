# Exactly what is left to run before submission

Everything in the manuscript is backed by executed code EXCEPT the items
below. Numbers currently in the paper that depend on step 1 are marked.

## 1. S&P 500 raw data (REQUIRED to regenerate two result blocks fresh)
The per-window S&P analyses were computed in an earlier session and are
loaded from caches (`results/exp3.json`, `results/exp4.json`) because the
29 MB csv could not be re-downloaded in the build environment.
- Download: see `data/README.md` (Kaggle "S&P 500 stock data",
  `all_stocks_5yr.csv`), place at `data/all_stocks_5yr.csv`.
- Run: `cd code && python exp3_real.py && python exp4_suff.py all`
- Effect: regenerates the S&P panel of Figure 3, the S&P rows of
  Table 1 (shares 0.25/0.75/0.00/0.25 raw and 0.89/1.00/0.00/0.00 after
  G; lower-tail shares 1.00 and 0.68), and the Figure 4(c) LOO R^2
  values (0.99/0.93/0.82/0.25/0.03). Values should reproduce the caches
  up to Monte Carlo noise; update the manuscript if any digit moves.

## 2. Full clean-machine reproduction (RECOMMENDED for provenance)
- `pip install -r requirements.txt` on a fresh machine,
  then `cd code && bash make_all.sh` (~2 to 4 hours single-core;
  now includes `exp6_baselines.py`, the competing-null calibration of
  Supplement Table S2).
- Confirms every figure and table regenerates from one command and
  refreshes `results/session_info.txt`.

## 3. Real ABIDE cohort (OPTIONAL upgrade to Section 7.4)
- `pip install nilearn && python fetch_real_data.py abide` (~1 GB), then
  `python exp4_abide.py` (header marks it [TO RUN]; it was never
  executed here). If adopted, rewrite Section 7.4 around
  `results/exp4_abide.json` and drop the synthetic-cohort labeling.

## 4. Recompile both documents
- `pdflatex main; bibtex main; pdflatex main; pdflatex main;`
  then the same cycle for `supplement` (xr reads `main.aux`).
- Produce the unblinded PDF by setting `\newcommand{\blind}{0}` in
  `main.tex`; fill the Acknowledgments and the AI-disclosure
  [model version] placeholder there.

## 5. Gates that REQUIRE a networked machine (not runnable in this build)
- Citation verification: resolve every entry in `refs.bib` against
  DOI/Crossref/arXiv. IMPORTANT: five entries were added from the
  referee reports with UNVERIFIED metadata (marked in refs.bib):
  langsrud2005rotation, perry2010rotation, hoff2008eigenmodel,
  dobriban2022consistency, koning2024more. Verify volume/pages/venues
  or correct them. Additionally verify and insert, at the marked
  positions in Section 1.2: Mastrandrea et al. (2017)
  spectrum-preserving rotation null for connectomes (heuristics
  sentence, next to vasa2022null); Pickrell (1991) invariant-matrix
  classification (attribution sentence, Section 3.1); Chiu, Sharp &
  Bloem-Reddy (invariance-testing sentence); Koussis et al.
  eigenstrapping (heuristics sentence); MacMahon & Garlaschelli RMT
  modularity (S&P modularity finding, Section 7.3); Bun, Bouchaud &
  Potters (finance line).
- Fetch the official ASA/T&F JASA LaTeX template and its agsm-style
  .bst; move the section files in unchanged, switch
  `\bibliographystyle`, delete the stand-in preamble comment, and
  re-verify the page ceiling in the official template.
- Check the live JASA Instructions for Authors and the ACC
  reproducibility form against `reproducibility/ACC_form_draft.md`.
- Verify the exact chapter/section numbering of the cited Kallenberg
  (2005) representation theorem against a library copy.

## 6. Author-only items
- Fill author block, ORCID, funding, Acknowledgments (unblinded).
- Record the AI model version in the disclosure section.
- Sign off on the cover letter (`reproducibility/cover_letter.md`).

## Note on the corrected proofs (this revision)
The Supplement (50 pp) now contains full self-contained proofs of Theorems 1-4 and Proposition 1, plus Remarks S1-S2 (mutual singularity for hollow/fixed-diagonal observations). Main text is exactly 35 pages under the documented double-spacing stand-in with references set \small at single spacing; re-measure both against the official JASA class and restore the class's own reference styling there.
- exp7_nodewise.py is part of make_all.sh (about one minute; seed 20260870); rerun with everything else on the clean machine.
- Stand-in page count after the Proposition 2 / Corollary S2 / nodewise additions: main 37 (33 text + references), supplement 53. Re-measure in the official JASA class; if over there, the ranked cuts are the Section 6 diagnostic-readings paragraph, related-work citation lists, Section 7.1 narration.

## Referee-driven items queued (require compute or network)
1. Two-sided p-values: rerun exp3/exp5 with the lower-tail counter added (one flag in conditional_pvalues) and report 2min(p,ptilde) alongside upper-tail values.
2. Scale up calibration: exp2 with R>=2000 null replicates and B=999 for the multiplicity analyses; report Wilson intervals (script stages exist).
3. Real-data benchmark suite on HCP/S&P/air: rewiring, eigenstrapping-style frame randomization, fitted spiked bootstrap (exp6 machinery reusable).
4. Regenerate Figure 3 (stray panel text) and Figure 4 (uncertainty bars) via exp3/exp4 plot stages.
5. Hill sensitivity analysis over thresholds for exp5.
6. Per-statistic verification of Proposition 2 hypotheses, or a redesigned provably Lipschitz battery.
7. Real ABIDE run as previously specified.
