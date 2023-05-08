## usethis namespace: start
#' @useDynLib inspre, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
## usethis namespace: end
NULL


off_diagonal <- function(X) {
  D <- ncol(X)
  return(X[!diag(D)])
}


inspre_primal <- function(X, U, V, lambda, gamma) {
  x_dist <- .5 * sum((X - U)^2, na.rm = TRUE)
  l1_term <- lambda * sum(abs(off_diagonal(V)))
  det_term <- -1 * gamma * determinant(V)$modulus
  return(c(x_dist, l1_term, det_term))
}


inspre_lagrangian <- function(X, U, V, theta, rho, lambda, gamma) {
  VU_minus_I <- V %*% U - diag(nrow(U))
  VU_term1 <- sum(theta * VU_minus_I)
  VU_term2 <- .5 * rho * sum(VU_minus_I^2)
  return(c(inspre_primal(X, U, V, lambda, gamma), VU_term1, VU_term2))
}


#' Helper function to do basic filtering of the TCE matrix.
#'
#' Large values of R, entries with a high SE, and row/columns with many nans
#' can be removed.
#'
#' @param R_tce Matrix or data.frame. Estimates of TCE.
#' @param SE_tce Matrix or data.frame. Standard errors of the entries in R_tce.
#' @param max_R Float. Set all entries where `abs(R_tce) > max_R` to `NA`.
#' @param max_SE Float. Set all entries whwere `SE > max_SE` tp `NA`.
#' @param max_nan_perc Float. Remove columns and rows that are more than
#'   `max_nan_perc` NAs.
#' @export
filter_tce <- function(R_tce, SE_tce, max_R = 1, max_SE = 0.5, max_nan_perc = 0.5) {
  R_tce[is.nan(SE_tce)] <- NA
  SE_tce[is.nan(SE_tce)] <- NA

  R_too_large <- abs(R_tce) > max_R
  R_tce[R_too_large] <- NA
  SE_tce[R_too_large] <- NA
  SE_too_large <- SE_tce > max_SE
  R_tce[SE_too_large] <- NA
  SE_tce[SE_too_large] <- NA

  row_nan_perc <- rowMeans(is.na(R_tce))
  col_nan_perc <- colMeans(is.na(R_tce))
  max_row_nan = max(row_nan_perc)
  max_col_nan = max(col_nan_perc)
  while((max_row_nan > max_nan_perc) | (max_col_nan > max_nan_perc)){
    if(max_row_nan >= max_col_nan){
      which_max_row_nan <- which.max(row_nan_perc)
      R_tce <- R_tce[-which_max_row_nan, -which_max_row_nan]
      SE_tce <- SE_tce[-which_max_row_nan, -which_max_row_nan]
    } else{
      which_max_col_nan <- which.max(col_nan_perc)
      R_tce <- R_tce[-which_max_col_nan, -which_max_col_nan]
      SE_tce <- SE_tce[-which_max_col_nan, -which_max_col_nan]
    }
    row_nan_perc <- rowMeans(is.na(R_tce))
    col_nan_perc <- colMeans(is.na(R_tce))
    max_row_nan = max(row_nan_perc)
    max_col_nan = max(col_nan_perc)
  }
  return(list("R" = R_tce, "SE" = SE_tce))
}


#' Helper function to make weights for inspre.
#'
#' @param SE DxD matrix of standard errors.
#' @param max_med_ratio Float > 1. Ratio of maximum weight to minimum non-zero
#'   weight. Improves conditioning when some SEs are very small.
#' @export
make_weights <- function(SE, max_med_ratio = NULL) {
  weights <- 1 / SE^2
  weights[is.na(SE)] <- 0

  if (is.null(max_med_ratio)) {
    infs <- is.infinite(weights)
    weights[infs] <- 0
    max_weight <- max(weights)
    weights[infs] <- max_weight
  } else {
    median_weight <- stats::median(weights)
    max_weight <- median_weight * max_med_ratio
    weights[weights > max_weight] <- max_weight
  }

  weights <- weights / mean(weights)
  return(weights)
}


#' Worker function to fit inspre for a single value of lambda.
#'
#' @importFrom foreach %dopar%
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param W DxD Matrix of weights.
#' @param rho Float, learning rate for ADMM.
#' @param lambda Float. L1 regularization strength on inverse of X.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param warm_start List with elements "U" and "V" containing starting point.
#' @param symmetrize True to force the output to be symmetric. If the input
#'   is symmetric, the output isn't always perfectly symmetric due to numerical
#'   issues.
#' @param verbose Bool. True to print solution progress.
#' @param gamma Float. Determinant penalty strength.
#' @param mu rho modification parameter for ADMM. Rho will be
#'   increased/decreased when the dual constrant and primal constraint are off
#'   by a factor of > mu.
#' @param tau rho modification parameter for ADMM. When called for, rho will be
#'   increased/decreased by the factor tau.
#' @param solve_its Integer, number of iterations of bicgstab/lasso to run
#'   for each U and V update.
#' @param ncores Integer, number of cores to use.
inspre_worker <- function(X, W = NULL, rho = 1.0, lambda = 0.01,
                          its = 100, delta_target = 1e-4, warm_start = NULL,
                          symmetrize = FALSE, verbose = FALSE, gamma = 0.0,
                          mu = 5, tau = 1.5, solve_its = 3, ncores = 1) {
  D <- nrow(X)
  if(is.null(W)){
    W <- matrix(1L, nrow = D, ncol = D)
  }
  theta <- matrix(0.0, D, D)
  if (is.null(warm_start)) {
    V <- diag(D)
    U <- diag(D)
    U[is.na(U)] <- 0
  } else {
    V <- warm_start$V
    U <- warm_start$U
  }
  L <- sum(inspre_lagrangian(X, U, V, theta, rho, lambda, gamma))

  min_delta_met <- 2
  delta_met <- 0

  doMC::registerDoMC(cores = ncores)
  start_time <- Sys.time()
  for (iter in 1:its) {
    rvtv <- rho * t(V) %*% V
    WX <- W * X
    WX[is.na(WX)] <- 0
    rhs <- WX + rho * t(V) - t(V) %*% theta

    d <- NULL
    U_next <- foreach::foreach (d = 1:D, .combine = cbind) %dopar% {
      Rlinsolve::lsolve.bicgstab(
        A = diag(W[, d]) + rvtv,
        B = rhs[, d],
        xinit = U[, d],
        maxiter = solve_its,
        verbose = FALSE)$x
    }

    glm_X <- sqrt(rho) * t(U_next)
    glm_Y <- ((rho + gamma) * diag(D) - t(theta)) / sqrt(rho)

    d <- NULL
    V_next <- foreach::foreach (d = 1:D, .combine = rbind) %dopar% {
      penalty_factor <- rep(1, D)
      penalty_factor[d] <- 0
      V_row <- V[d, ]
      lasso(glm_X, glm_Y[, d], V_row, lambda = penalty_factor * lambda,
            niter = solve_its)
      V_row
    }

    if (symmetrize) {
      V_next <- 0.5 * (V_next + t(V_next)) *
        (abs(V_next) > 1e-10 & abs(t(V_next)) > 1e-10)
    }

    constraint_resid <- V_next %*% U_next - diag(D)
    dual_resid <- rho * t(V) %*% (V_next - V) %*% U_next
    theta <- theta + rho * constraint_resid

    constraint_resid_norm <- sqrt(mean(constraint_resid^2))
    dual_resid_norm <- sqrt(mean(dual_resid^2))
    if(constraint_resid_norm > mu*dual_resid_norm){
      rho <- tau*rho
    } else if(dual_resid_norm > mu*constraint_resid_norm){
      rho <- rho/tau
    }

    rmsd_u <- sqrt(mean((U - U_next)^2))
    rmsd_v <- sqrt(mean((V - V_next)^2))
    U <- U_next
    V <- V_next
    lagrangian <- inspre_lagrangian(X, U, V, theta, rho, lambda, gamma)
    L_next <- sum(lagrangian)
    L_delta <- abs(L_next - L) / L
    L <- L_next
    if (L_delta > delta_target) {
      delta_met <- 0
    } else {
      delta_met <- delta_met + 1
    }

    if (verbose) {
      cat(iter, L, rho, constraint_resid_norm, dual_resid_norm, rmsd_u, rmsd_v,
          L_delta, Sys.time() - start_time, "\n")
    }

    # Converged if some iterations in a row have delta < target.
    if (delta_met >= min_delta_met) {
      break
    }
  }
  end_time <- Sys.time()
  return(list(
    "V" = V, "U" = U, "L" = L, "F_term" = lagrangian[1],
    "l1_term" = lagrangian[2], "det_term" = lagrangian[3], "rho" = rho,
    "L_delta" = L_delta, "constraint_resid" = constraint_resid_norm,
    "dual_resid" = dual_resid_norm, "iter" = iter,
    "time" = end_time - start_time))
}


#' Fits inspre model for sequence of lambda values.
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param W DxD Matrix of weights.
#' @param rho Float. Initial learning rate for ADMM.
#' @param lambda Float or sequence of floats. Path of L1 regularization strength
#'   on inverse of X.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param symmetrize True to force the output to be symmetric. If the input
#'   is symmetric, the output isn't always perfectly symmetric due to numerical
#'   issues.
#' @param verbose 0, 1 or 2. 2 to print convergence progress for each lambda,
#'   1 to print convergence result for each lambda, 0 for no output.
#' @param gamma Float or sequence of floats. Determinant penalty strength.
#' @param train_prop Proportion of the dataset to use to train the model. Set
#'   to 1.0 to use all data.
#' @param mu rho modification parameter for ADMM. Rho will be
#'   increased/decreased when the dual constrant and primal constraint are off
#'   by a factor of > mu.
#' @param tau rho modification parameter for ADMM. When called for, rho will be
#'   increased/decreased by the factor tau.
#' @param solve_its Integer, number of iterations of bicgstab/lasso to run
#'   for each U and V update.
#' @param ncores Integer, number of cores to use.
fit_inspre_sequence <- function(X, lambda, W = NULL, rho = 1.0,
                                its = 100, delta_target = 1e-4,
                                symmetrize = FALSE, verbose = 1,
                                gamma = 0.0, train_prop = 1.0,
                                mu = 5, tau = 2.5, solve_its = 3, ncores = 1,
                                warm_start = TRUE) {
  # Break the matrix into training and test sets, equally and at random.
  D <- nrow(X)
  if(is.null(gamma)){
    gamma <- 0.0
  }
  if(train_prop == 1.0){
    train_W <- W
    test_W <- numeric(0)
  } else {
    if(is.null(W)){
      W <- matrix(1L, nrow=D, ncol=D)
    }
    rand <- matrix(stats::runif(D * D), D)
    train <- rand < train_prop
    test <- !train & (W > 0)
    diag(train) = TRUE
    diag(test) = TRUE
    train_W <- W
    train_W[!train] <- 0
    train_W <- train_W/mean(train_W)
    test_W <- W
    test_W[!test] <- 0
    test_W <- test_W/mean(test_W)
  }

  V_all <- array(dim = c(D, D, length(lambda)))
  U_all <- array(dim = c(D, D, length(lambda)))
  j <- 1
  rho_used <- vector("numeric", length = length(lambda))
  L_path <- vector("numeric", length = length(lambda))
  F_term_path <- vector("numeric", length = length(lambda))
  train_error <- vector("numeric", length = length(lambda))
  test_error <- vector("numeric", length = length(lambda))
  warm_start_res = NULL
  for (i in 1:length(lambda)) {
    lambda_i <- lambda[i]
    if(length(gamma) > 1){
      gamma_i = gamma[i]
    } else {
      gamma_i = gamma
    }
    inspre_res <- inspre_worker(
      X = X, W = train_W, lambda = lambda_i, rho = rho, warm_start = warm_start_res,
      its = its, delta_target = delta_target, symmetrize = symmetrize,
      verbose = (verbose == 2), gamma = gamma_i, mu = mu, tau = tau,
      solve_its = solve_its, ncores = ncores)
    if (verbose) {
      cat(sprintf(paste("Converged to L = %f in %d iterations for lambda = %f.",
                        " Time %f.",
                        "\n L_delta = %e, constraint_resid = %e, dual_resid = ",
                        "%e.\n", sep = ""),
                  inspre_res$L, inspre_res$iter, lambda_i, inspre_res$time, inspre_res$L_delta,
                  inspre_res$constraint_resid, inspre_res$dual_resid))
    }
    if(warm_start){
      warm_start_res <- inspre_res
    }
    V_all[, , i] <- inspre_res$V
    U_all[, , i] <- inspre_res$U
    rho_used[i] <- inspre_res$rho
    L_path[i] <- inspre_res$L
    if(is.null(train_W)){
      train_error[i] <- sqrt(
        mean(off_diagonal(X - inspre_res$U)^2, na.rm = TRUE))
    } else{
      train_error[i] <- sqrt(
        mean(off_diagonal(train_W * (X - inspre_res$U))^2, na.rm = TRUE))
    }
    test_error[i] <- sqrt(
      mean(off_diagonal(test_W * (X - inspre_res$U))^2, na.rm = TRUE))
  }
  return(list("V" = V_all, "U" = U_all, "lambda" = lambda, "gamma" = gamma, "rho" = rho_used,
              "L" = L_path, train_error = train_error, test_error = test_error))
}


#' Finds the inverse of X using Inverse Sparse Regression.
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param W DxD Matrix of weights.
#' @param rho Float. Initial learning rate for ADMM.
#' @param lambda Float, sequence of floats of NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmicallly spaced set of values between
#'   the maximimum absolute off diagonal element of X and lambda_min_ratio
#'   times this value will be used.
#' @param lambda_min_ratio Float, ratio of maximum lambda to minimum lambda.
#' @param nlambda Integer. Number of lambda values to try.
#' @param alpha Float between 0 and 1 or NULL. If > 0, the model will be fit
#'   once with gamma = 0 to find L0, then all subsequent fits will use
#'   gamma = alpha * L0 / D. Set to NULL to provide gamma directly.
#' @param gamma Float or sequence of nlambda floats or NULL. Determinant
#'   regularization strength to use (for each lambda value). It is recommended
#'   to set alpha rather than setting this directly.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param symmetrize True to force the output to be symmetric. If the input
#'   is symmetric, the output isn't always perfectly symmetric due to numerical
#'   issues.
#' @param verbose 0, 1 or 2. 2 to print convergence progress for each lambda,
#'   1 to print convergence result for each lambda, 0 for no output.
#' @param train_prop Float between 0 and 1. Proportion of data to use for
#'   training in cross-validation.
#' @param cv_folds Integer. Number of cross-validation folds to perform.
#' @param mu rho modification parameter for ADMM. Rho will be
#'   increased/decreased when the dual constrant and primal constraint are off
#'   by a factor of > mu.
#' @param tau rho modification parameter for ADMM. When called for, rho will be
#'   increased/decreased by the factor tau.
#' @param solve_its Integer, number of iterations of bicgstab/lasso to run
#'   for each U and V update.
#' @param ncores Integer, number of cores to use.
#' @export
inspre <- function(X, W = NULL, rho = 10.0, lambda = NULL,
                   lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                   gamma = NULL, its = 100, delta_target = 1e-4,
                   symmetrize = FALSE, verbose = 1, train_prop = 0.8,
                   cv_folds = 0, mu = 10, tau = 2, solve_its = 3, ncores = 1,
                   warm_start = TRUE, min_nz = 1e-5) {
  D <- ncol(X)

  if(any(is.na(X))){
    if(is.null(W)){
      W <- matrix(1L, nrow = D, ncol = D)
    }
    W[is.na(X)] <- 0
  }

  if (is.null(lambda)) {
    lambda_max <- min(max(abs(off_diagonal(X)), na.rm = TRUE), 1.0)
    lambda_min <- lambda_min_ratio * lambda_max
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  }

  if(verbose){
    cat("Fitting model with full dataset. \n")
  }
  if(!is.null(alpha)){
    gamma <- 0
    if(alpha > 0){
      if(verbose){
        cat("Fitting with gamma = 0. \n")
      }
      L0 <- fit_inspre_sequence(
        X = X, W = W, rho = rho, lambda = lambda, delta_target = delta_target,
        symmetrize = symmetrize, verbose = verbose, gamma = 0, its = its,
        train_prop = 1.0, mu = mu, tau = tau, solve_its = solve_its,
        ncores = ncores, warm_start = warm_start)$L

      if(verbose){
        cat("Refitting with gamma = alpha * L / D. \n")
      }
      gamma <- alpha * L0 / D
    }
  }
  full_res <- fit_inspre_sequence(
    X = X, W = W, rho = rho, lambda = lambda, delta_target = delta_target,
    symmetrize = symmetrize, verbose = verbose, gamma = gamma, its = its,
    train_prop = 1.0, mu = mu, tau = tau, solve_its = solve_its,
    ncores = ncores, warm_start = warm_start)

  if (cv_folds > 0) {
    xi_mat <- array(0, dim = c(D, D, length(lambda)))
    error_matrix <- array(dim = c(cv_folds, length(lambda)))
    for (i in 1:cv_folds){
      if (verbose){
        cat(sprintf("Cross-validation iteration %d.\n", i))
      }
      cv_res <- fit_inspre_sequence(
        X = X, W = W, rho = rho, lambda = lambda, delta_target = delta_target,
        symmetrize = symmetrize, verbose = verbose, gamma = gamma,
        train_prop = train_prop, its = its, mu = mu, tau = tau,
        solve_its = solve_its, ncores = ncores, warm_start = warm_start)
      error_matrix[i, ] = cv_res$test_error
      V_nz <- abs(cv_res$V) > min_nz
      xi_mat <- xi_mat + V_nz
    }
    xi_mat <- xi_mat/cv_folds
    D_hat <- 2 * xi_mat * (1-xi_mat)
    D_hat <- apply(D_hat, 3, mean)
    D_hat_se <- apply(
      xi_mat, 3, function(x){ stats::sd(x)/sqrt(length(x)) })

    test_error <- apply(error_matrix, 2, mean)
    test_error_se <- apply(
      error_matrix, 2, function(x){ stats::sd(x)/sqrt(length(x))})
    full_res$test_error_se <- test_error_se
    full_res$test_error <- test_error
    full_res$D_hat <- D_hat
    full_res$D_hat_se <- D_hat_se
    full_res$xi_mat <- xi_mat
  }
  return(full_res)
}


#' Fits inverse sparse regression model.
#'
#' See also inspre::inspre() for more details.
#'
#' @param R_tce D x D matrix of  "total causal effects".
#' @param W DxD Matrix of weights.
#' @param rho Float. Initial learning rate for ADMM.
#' @param lambda Float, sequence of floats of NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmicallly spaced set of values between
#'   the maximimum absolute off diagonal element of X and lambda_min_ratio
#'   times this value will be used.
#' @param lambda_min_ratio Float, ratio of maximum lambda to minimum lambda.
#' @param nlambda Integer. Number of lambda values to try.
#' @param alpha Float between 0 and 1 or NULL. If > 0, the model will be fit
#'   once with gamma = 0 to find L0, then all subsequent fits will use
#'   gamma = alpha * L0 / D. Set to NULL to provide gamma directly.
#' @param gamma Float or sequence of nlambda floats or NULL. Determinant
#'   regularization strength to use (for each lambda value). It is recommended
#'   to set alpha rather than setting this directly.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param verbose 0, 1 or 2. 2 to print convergence progress for each lambda,
#'   1 to print convergence result for each lambda, 0 for no output.
#' @param train_prop Float between 0 and 1. Proportion of data to use for
#'   training in cross-validation.
#' @param cv_folds Integer. Number of cross-validation folds to perform.
#' @param mu rho modification parameter for ADMM. Rho will be
#'   increased/decreased when the dual constrant and primal constraint are off
#'   by a factor of > mu.
#' @param tau rho modification parameter for ADMM. When called for, rho will be
#'   increased/decreased by the factor tau.
#' @param solve_its Integer, number of iterations of bicgstab/lasso to run
#'   for each U and V update.
#' @param ncores Integer, number of cores to use.
#' @param warm_start Logical. Whether to use previous lambda value result as
#'   starting point for next fit.
#' @export
fit_inspre <- function(R_tce, W = NULL, rho = 10.0, lambda = NULL,
                       lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                       gamma = NULL, its = 100, delta_target = 1e-4,
                       verbose = 1, train_prop = 0.8,
                       cv_folds = 0, mu = 10, tau = 2, solve_its = 3,
                       ncores = 1, warm_start = TRUE, min_nz = 1e-5){
  D <- dim(R_tce)[1]
  inspre_res <- inspre::inspre(
    X = R_tce, W = W, rho = rho, lambda = lambda,
    lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
    gamma = gamma, its = its, delta_target = delta_target, symmetrize = FALSE,
    verbose = verbose, train_prop = train_prop, cv_folds = cv_folds, mu = mu,
    tau = tau, solve_its = solve_its, ncores = ncores, warm_start = warm_start, min_nz = min_nz)
  inspre_res$R_hat <- array(0L, dim = dim(inspre_res$V))
  for(i in 1:length(inspre_res$lambda)){
    inspre_res$R_hat[ , , i] <-
      diag(D) - inspre_res$V[ , , i] / diag(inspre_res$V[ , , i])
  }
  dimnames(inspre_res$R_hat) <- list(rownames(R_tce), colnames(R_tce), inspre_res$lambda)
  return(inspre_res)
}
