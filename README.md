# Rotatable Networks: code, data, and manuscript

Reproduction package for "Rotatable Networks" (JASA Theory & Methods
submission).
Structure follows the jasa-acs reproducibility template: one script per
figure/table, pinned environment, recorded session info, fixed seeds.

## Layout
- `main.tex`, `sec_*.tex`, `refs.bib`  -- main manuscript (blinded via `\blind`)
- `supplement.tex`, `supp_*.tex`       -- supplement: deferred statements,
  all proofs, Table S1, Figure S1, data appendix
- `figures/`   -- vector PDFs used by both documents
- `code/`      -- experiment scripts (Python 3.12; see requirements.txt)
- `data/`      -- shipped datasets + `data/README.md` access instructions
- `results/`   -- JSON outputs backing every reported number + session_info.txt
- `reproducibility/` -- RUN_BEFORE_SUBMISSION.md, ACC form draft, cover
  letter, notation ledger, claims table

## One-command reproduction
    pip install -r requirements.txt
    # place data/all_stocks_5yr.csv (see data/README.md) for the S&P blocks
    cd code && bash make_all.sh        # ~2-4 h single core

Script -> output mapping: `exp1_bbp.py` -> Fig 1; `exp2_test.py` +
`exp2b_graphon_strong.py` -> Fig 2 and the strong-graphon number;
`exp2c_corr_level.py` -> Table S1; `exp3_real.py` -> Fig 3 + Table 1;
`exp4_suff.py` -> Fig 4; `exp5_heavy.py` -> Fig S1; `exp6_baselines.py` -> Table S2 (competing
null calibration). `common.rotatable_test(X)` is the one-call
interface: observed value, Haar-conditional mean/sd, residual z,
upper-tail p per statistic, and the exact joint min-p combination. Staged scripts
(`exp1`, `exp2`, `exp2c`, `exp4`) accept a stage argument for
constrained environments, or `all` to run every stage in sequence; the
stage lists are at the top of each script. All seeds are fixed
in-script; `results/session_info.txt` records the environment that
produced the shipped results.

`code/exp4_abide.py` is marked [TO RUN]: it replicates the cohort
analysis on real ABIDE data after `python fetch_real_data.py abide` and
was never executed in the build environment.

## Compiling the manuscript
    pdflatex main && bibtex main && pdflatex main && pdflatex main
    pdflatex supplement && bibtex supplement && pdflatex supplement && pdflatex supplement
Compile `main` first: the supplement resolves cross-references through
`main.aux` (xr). `\newcommand{\blind}{0}` in `main.tex` produces the unblinded
version. At template adoption, move the section files into the
official ASA/T&F JASA class and switch the bibliography style to its
agsm-style .bst; `apalike` is the stand-in author-year style in this
build. Manuscript sources are maintained separately; this repository
contains the code, data, and results.

## Before submitting
Work through `reproducibility/RUN_BEFORE_SUBMISSION.md`; it lists the
exact remaining runs (S&P csv re-run, optional ABIDE, clean-machine
make_all), the networked verification gates (citations, official
template, live ACC form), and the author-only items.
