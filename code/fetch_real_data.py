"""Fetch the full real datasets used (or referenced) in the paper.

Run this on an unrestricted machine; the analysis sandbox used for the
included figures could reach GitHub-hosted mirrors only (see README).

1. ABIDE preprocessed connectomes (subject-level; no DUA required)
   -> per-subject CC200 ROI time series, then Fisher-z correlation matrices.
2. HCP group-average FC/SC (ENIGMA toolbox mirrors; already included).
3. Kenneth French 49-industry daily returns (alternative market panel).
4. BACI / CEPII world-trade flows (heavy-tailed; registration-free bulk file).

Usage:  python fetch_real_data.py abide|french|baci
"""
import sys, os, io, zipfile, urllib.request

import numpy as np


def fetch_abide(n_subjects=200, out="abide_fc"):
    from nilearn.datasets import fetch_abide_pcp
    from nilearn.connectome import ConnectivityMeasure
    os.makedirs(out, exist_ok=True)
    d = fetch_abide_pcp(derivatives=["rois_cc200"], pipeline="cpac",
                        band_pass_filtering=True, global_signal_regression=True,
                        n_subjects=n_subjects)
    cm = ConnectivityMeasure(kind="correlation")
    mats = cm.fit_transform(d.rois_cc200)
    np.save(os.path.join(out, "fc_cc200.npy"), mats)
    np.save(os.path.join(out, "dx_group.npy"),
            np.array(d.phenotypic["DX_GROUP"]))
    print("saved", mats.shape, "->", out)


def fetch_french(out="french49.csv"):
    url = ("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/"
           "ftp/49_Industry_Portfolios_daily_CSV.zip")
    z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(url).read()))
    name = z.namelist()[0]
    open(out, "wb").write(z.read(name))
    print("saved", out)


def fetch_baci(out="baci"):
    print("Download BACI (HS17, latest release) from "
          "https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37 "
          "then aggregate bilateral flows: exporter x importer total value.")


if __name__ == "__main__":
    {"abide": fetch_abide, "french": fetch_french,
     "baci": fetch_baci}[sys.argv[1]]()
