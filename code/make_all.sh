#!/usr/bin/env bash
# Master reproduction script: one command per figure/table, in paper order.
# Runs from code/. On a clean machine: pip install -r ../requirements.txt
# Stage-by-stage invocations (for constrained environments) are documented
# in ../README.md; each staged script also accepts a single "all" argument.
set -euo pipefail

python session_info.py

python exp1_bbp.py all            # Figure 1
python exp2_test.py all           # Figure 2 (level + power)
python exp2b_graphon_strong.py    # strong-graphon power quoted in Sec 7.2
python exp2c_corr_level.py all    # finite-T correlation size (Table S1)
python exp3_real.py               # Figure 3 + Table 1 (S&P block needs
                                  #   ../data/all_stocks_5yr.csv; see data/README.md)
python exp4_suff.py all           # Figure 4 (S&P panel needs the same csv)
python exp5_heavy.py              # Supplement Figure S1

echo "All experiments finished; results in ../results, figures in ../figures."
python exp7_nodewise.py
