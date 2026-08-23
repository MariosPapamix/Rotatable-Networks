# Data files

Included (small, redistributable):
- funcMatrix_ctx_schaefer_200.csv, funcMatrix_ctx.csv,
  strucMatrix_ctx_schaefer_200.csv — HCP group-average connectivity
  (ENIGMA toolbox, MICA-MNI/ENIGMA repository)
- flights_airport.csv — US airport-pair flight counts (vega-datasets)

Not included (size): the S&P 500 daily price panel. Download with

    curl -o all_stocks_5yr.csv \
      https://raw.githubusercontent.com/plotly/datasets/master/all_stocks_5yr.csv

before running exp3_real.py / exp4_suff.py.

Subject-level ABIDE and other full datasets: see ../code/fetch_real_data.py.
