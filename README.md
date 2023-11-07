# Inverse Sparse Regression (`inspre`)

`inspre` is an R package for calculating a sparse approximate matrix inverse with
per-entry weights. Applications include directed and undirected graphical
model inference. For more information on the method please consult the paper
[Large-scale causal discovery using interventional data sheds light on the regulatory network structure of blood traits](https://www.biorxiv.org/content/10.1101/2023.10.13.562293v1).

## Installation
At the moment the package must be installed using `devtools::install_github("brielin/inspre")`.
Eventually we intend to make the package available on CRAN.

If you want to run tests, you will need to add the install option to download
the tests as well. 

```
> devtools::install_github("brielin/inspre", INSTALL_opts="--install-tests")
> library(testthat)
> library(inspre)
> test_package("inspre")
```

## Requirements
`inspre` was developed using:

- R (4.2.1)
- Rcpp (1.0.11)
- RcppEigen (0.3.3.9.3)
- RcppArmadillo (0.12.4.1.0)
- foreach (1.5.2)
- doMC (1.3.8)
- dplyr (1.1.2)
- purrr (1.0.2)

To use functions dependent on h5, we require

- hdf5r (1.3.8)

To use cross-validation with sample-level data, we require

- caret (6.0.94)

To run and reproduce the GWPS analysis in the manuscript (`Rmd/run_gwps_analysis.R`, `Rmd/analyze_inspre_gwps.Rmd`), we require

 - ggplot2 (3.4.2)
 - tidyr (1.3.0)
 - tibble (3.2.1)
 - igraph (1.2.11)
 - ggupset (0.3.0)
 - betareg (3.1.4)
 - readr (2.1.4)

The simulation code is located at `Rmd/simulate_*.Rmd`. It is designed to be run in parallel on a slurm cluster using `rslurm`. To run and reproduce the simulations in the manuscript, we require

- pcalg (2.7.8)
- rslurm (0.6.2)
- ggplot2 (3.4.2)

Several of the comparison methods need to be run in Python. Running these requires

- python (3.7.8)
- numpy
- pandas
- sklearn
- tensorflow

The parenthetical version indicates the version used in development. This is not necessarily the required version but one that is guaranteed to work.

## Use
The primary interface to the inverse sparse regression method is the funcion `inspre::inspre`. We provide several helper functions for data in different formats, performing cross-validation, etc.

- `inspre::fit_inspre_from_R` if you already have ACE data in a `features x features` matrix and just need to fit a network.
- `inspre::fit_inspre_from_X` if you have a `samples x features` data matrix with a set of `samples x interventions` target indicators and you need to learn the effects of the interventions and the network.
- `inspre::fit_inspre_from_h5X` similar to the above but with the data matrix stored in an hdf5 file as is commonly distributed with Perturb-seq data. See `Rmd/run_gwps_analysis.R` for an example using this interface.

For now see the function documentation for additional usage information. This `README` will be updated with more details as the manuscript and analyses are finalized.
