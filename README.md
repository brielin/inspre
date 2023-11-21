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

Some noteable parameters you could consider setting in your analysis that are common across the fit functions
- `verbose` - `0`, `1`, or `2`. Sets the verbosity level. The default (`1`) prints the convergence statistics after every `lambda`. When assessing whether convergence was successful, you may consider setting this to `2` and noting whether convergence was very slow (hitting the iteration limit) or too
fast (diverging very quickly for many `lambda`) values.
- `weighted` and `max_med_ratio` - Set `weighted = TRUE` to use a weighted analysis (default). We generally recommend using weights with real data as standard errors can vary hugely per entry in the R matrix. On the other hand, sometimes some standard errors are so tiny that inverse variance weighting puts huge weights on just a handful of entries. You can set `max_med_ratio` to circumvent this and ensure that the maximum weight is only this times the median weight. By default this is set to `NULL`, or no limit. In our GWPS analysis we set this to `50`.
- `lambda`, `nlambda`, `lambda_min_ratio`. By default, `inspre` tests 20 `lambda` values ranging from the maximum off-diagonal element of the `R` matrix down to this number times `lambda_min_ratio`. For example, if the max off diagonal is 0.8 and `lambda_min_ratio` is 0.01 (default), the max tested `lambda`  will be 0.8 and the min will be 0.008. If you want to set your own `lambda` sequence you can supply it via the `lambda=` parameter and the other two parameters will be ignored.
- `cv_folds`. Set to a number `k` > 0 to do k-fold cross-validation. If so, several cross-validation performance metrics will be returned in the result object. `eps_hat_G`, `eps_hat_I`, and `eps_hat_Rte` can all be used to evaluate held-out sample performance. They represent 3 different ways of computing cross validation error and will not always necessarily agree on which `lambda` is the best.
- `DAG`. Set to `TRUE` to enforce an approximate DAG constraint. Note that this is generally only recommended to improve convergence behaviors if your fit is quickly diverging when DAG is set to `FALSE`. Our tests show that even if the underlying model is a DAG, the unrestricted fit generally produces better results. However, the restricted model is easier to fit, so this is useful if your data is challenging.
- `rho`. This is the primary convergence parameter. This is set very high by default (`100`) so that the model focuses on keeping the inverse constraint on difficult datasets, but this can be slow if the data is easier. If you are seeing very slow convergence behavior (hitting the iteration limit well before converging), then you can try reducing this. In our GWP analysis we use `100`, in simulated data we use `10`.

## Example.
You will find an example data matrix, targets list, true G and true R under the folder `tests/example_data/`. To run this minimal example in `R`:

```
library(inspre)
X <- as.matrix(read.table('tests/example_data/example_dataset.txt', header=TRUE))
targets <- scan('tests/example_data/example_targets.txt', character())
inspre_res <- inspre::fit_inspre_from_X(X, targets)
```

This `README` will be updated with more details as the manuscript and analyses are finalized.
