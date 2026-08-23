# Overleaf project: Rotatable Networks

Two compile targets in one project:
1. `main.tex`  (compile first: pdfLaTeX + BibTeX; blinded via `\blind`=1)
2. `supplement.tex` (compile second; it reads `main.aux` through the xr
   package, so run main to completion once before the supplement)

`figures/` and `refs.bib` are shared. Precompiled `main.bbl` and
`supplement.bbl` are included so the project also builds where BibTeX
runs are restricted. The bibliography style is `apalike` as a stand-in;
switch to the official ASA/T&F JASA class and its agsm-style .bst at
template adoption (one line each in the two preambles).
