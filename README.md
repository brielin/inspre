# Inverse Sparse Regression (inspre)

inspre is an R package for calculating a sparse approximate matrix inverse with
per-entry weights. Applications include directed and undirected graphical
model inference. For more information on the method please consult the paper
[Phenome-scale causal network discovery with bidirectional mediated Mendelian
randomization](https://www.biorxiv.org/content/10.1101/2020.06.18.160176v2).

## Installation
At the moment the package must be installed using `devtools::install_github()`.
Eventually we intend to make the package available on CRAN.

## Requirements
To run the main function `inspre`,
- R
- Rcpp
- RcppEigen
- foreach
- doMC
- Rlinsolve

To run tests, run `GGM.Rmd` and make plots,
- huge
- egg
- devtools
- testthat
