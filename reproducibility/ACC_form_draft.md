# ACC (Author Contributions Checklist) form -- DRAFT
# Verify field names against the live JASA/ACC form before submission.

## Data
- Abstract: Three data sources. (1) HCP group-average cortical functional
  connectome (Schaefer-200), shipped as csv from the ENIGMA toolbox
  release. (2) S&P 500 daily closes 2013-2018 (Kaggle "S&P 500 stock
  data", all_stocks_5yr.csv, ~29 MB, user-downloaded). (3) U.S. air
  traffic flow network (airport-pair flight counts), shipped as csv.
  ABIDE subject-level data are fetched by script (nilearn) and used only
  by the optional, not-yet-run exp4_abide.py.
- Availability: (1) and (3) included in the repository under data/ with
  provenance in data/README.md; (2) public download, link and checksum
  instructions in data/README.md; ABIDE public via nilearn fetcher.
- Format: csv (matrices), npy (fetched ABIDE arrays).

## Code
- Description: Python 3.12; numpy/scipy/pandas/scikit-learn/matplotlib,
  versions pinned in requirements.txt and recorded in
  results/session_info.txt. All Monte Carlo seeds fixed in-script.
- Structure: one script per figure/table (mapping in README.md);
  code/make_all.sh reproduces everything in order; staged execution
  supported for constrained environments.
- License: [choose: MIT / GPL-3].

## Instructions for reproduction
- pip install -r requirements.txt
- Place data/all_stocks_5yr.csv (see data/README.md)
- cd code && bash make_all.sh   (~2-4 h single core)
- Compile main.tex then supplement.tex (pdflatex + bibtex cycles).

## Notes
- results/exp3.json and results/exp4.json contain cached S&P blocks
  produced by the same scripts in an earlier session; step 1 of
  RUN_BEFORE_SUBMISSION.md regenerates them from the raw csv.
