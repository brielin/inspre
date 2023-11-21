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


inspre_primal <- function(X, W, U, V, lambda, gamma) {
  x_dist <- .5 * sum((W*(X - U))^2, na.rm = TRUE)
  l1_term <- lambda * sum(abs(off_diagonal(V)))
  det_term <- 0 # -1 * gamma * determinant(V)$modulus
  return(c(x_dist, l1_term, det_term))
}


inspre_lagrangian <- function(X, W, U, V, theta_UV, theta_VU, rhoUV, rhoVU, lambda, gamma) {
  UVmI <- U %*% V - diag(nrow(U))
  VUmI <- V %*% U - diag(nrow(U))
  sumUV <- sum(theta_UV * UVmI)
  sumVU <- sum(theta_VU * VUmI)
  sum_rhoUV <- 0.5 * rhoUV * sum(UVmI^2)
  sum_rhoVU <- 0.5 * rhoVU * sum(VUmI^2)
  return(c(inspre_primal(X, W, U, V, lambda, gamma), sumUV, sumVU, sum_rhoUV, sum_rhoVU))
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
filter_tce <- function(R_tce, SE_tce, max_R = 2, max_SE = 1, max_nan_perc = 0.5) {
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

  infs <- is.infinite(weights)
  if(!is.null(max_med_ratio)){
    med_weight <- stats::median(weights[!infs])
    max_weight <- med_weight * max_med_ratio
    weights[(weights > max_weight)&(!infs)] <- max_weight
  }
  max_weight <- max(weights[!infs])
  weights[infs] <- max_weight

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
#' @param constraint One of "UV" or "VU". Constraint to use.
#' @param DAG Bool. True to resitrict solutions to approximate DAGs.
inspre_worker <- function(X, W = NULL, rho = 1.0, lambda = 0.01,
                          its = 100, delta_target = 1e-4, warm_start = NULL,
                          symmetrize = FALSE, verbose = FALSE, gamma = 0.0,
                          mu = 5, tau = 1.5, solve_its = 3, ncores = 1,
                          constraint = "UV", DAG = FALSE) {
  D <- nrow(X)
  if(is.null(W)){
    W <- matrix(1.0, nrow = D, ncol = D)
  }
  if (is.null(warm_start)) {
    # TODO: consider initializing with psuedoinverse
    V <- diag(D)
    U <- diag(D)
    U[is.na(U)] <- 0
    theta <- matrix(0.0, D, D)
  } else {
    V <- warm_start$V
    U <- warm_start$U
    theta <- warm_start$theta
  }
  if(constraint == "UV"){
    L <- sum(inspre_lagrangian(X, W, U, V, theta, matrix(0.0, D, D), rho, 0, lambda, gamma))
  } else if(constraint == "VU"){
    L <- sum(inspre_lagrangian(X, W, U, V, matrix(0.0, D, D), theta, 0, rho, lambda, gamma))
  } else {
    stop("Not implemented")
  }
  min_delta_met <- 2
  delta_met <- 0
  max_wrong_way <- 3
  wrong_way <- 0

  start_time <- Sys.time()
  WX <- W * X
  WX[is.na(WX)] <- 0

  for (iter in 1:its) {
    if(constraint == "UV"){
      U_next <- fit_U_UV_const(WX, W, V, U, theta, rho, 0, solve_its, DAG, ncores)
      V_next <- fit_V_UV_const(V, U_next, theta, rho, lambda, solve_its, DAG, ncores)
      UV_next <- U_next %*% V_next

      constraint_resid <- UV_next - diag(D)
      non_constraint_resid <- V_next %*% U_next - diag(D)
      dual_resid <- rho * U_next %*% (V - V_next) %*% t(V_next)

      theta <- theta + rho * constraint_resid
      lagrangian <- inspre_lagrangian(X, W, U_next, V_next, theta, matrix(0L, D, D), rho, 0, lambda, gamma)

      constraint_norm <- max(sqrt(sum(UV_next**2)), sqrt(D))
      dual_norm <- max(sqrt(2*lagrangian[1]), sqrt(sum((theta %*% t(V))**2)))
    } else if(constraint == "VU"){
      U_next <- fit_U_VU_const(WX, W, V, U, theta, rho, 0, solve_its, DAG, ncores)
      V_next <- fit_V_VU_const(V, U_next, theta, rho, lambda, solve_its, DAG, ncores)
      VU_next <- V_next %*% U_next

      constraint_resid <- VU_next - diag(D)
      non_constraint_resid <- U_next %*% V_next - diag(D)
      dual_resid <- rho * t(V) %*% (V - V_next) %*% U_next

      theta <- theta + rho * constraint_resid
      lagrangian <- inspre_lagrangian(X, W, U_next, V_next, matrix(0L, D, D), theta, 0, rho, lambda, gamma)

      constraint_norm <- max(sqrt(sum(VU_next**2)), sqrt(D))
      dual_norm <- max(sqrt(2*lagrangian[1]), sqrt(sum((t(V) %*% theta)**2)))
    }
    L_next <- sum(lagrangian)
    if (symmetrize) {
      V_next <- 0.5 * (V_next + t(V_next)) *
        (abs(V_next) > 1e-10 & abs(t(V_next)) > 1e-10)
    }

    # tau_thresh <- sqrt(sqrt(mean(constraint_resid^2))/sqrt(mean(dual_resid^2)))
    # cat(tau_thresh, 1/tau_thresh, "\n")


    max_V <- max(abs(V_next))
    n_bad <- sum(abs(V_next) > 1.5)
    bad_sum <- sum(abs(V_next[abs(V_next) > 1.5]))
    L1V <- sum(abs(V_next))

    G <- diag(D) - V_next / diag(V_next)
    max_G <- max(abs(G))
    n_bad_G <- sum(abs(G) > 1.5)
    bad_sum_G <- sum(abs(G[abs(G) > 1.5]))

    max_U <- max(abs(U_next))
    n_bad_U <- sum(abs(U_next) > 1.5)
    bad_sum_U <- sum(abs(U_next[abs(U_next) > 1.5]))
    L2U <- sum(U_next**2)
    DU <- sum(diag(U_next))

    constraint_resid_norm <- sqrt(mean(constraint_resid^2))/constraint_norm
    non_constraint_resid_norm <- sqrt(mean(non_constraint_resid^2))
    dual_resid_norm <- sqrt(mean(dual_resid^2))/dual_norm
    if(constraint_resid_norm > mu*dual_resid_norm){
      rho <- tau*rho
    } else if(dual_resid_norm > mu*constraint_resid_norm){
      rho <- rho/tau
    }
    delta_sign <- sign(L_next - L)
    if(delta_sign == 1){
      wrong_way <- wrong_way + 1
    } else {
      wrong_way <- 0
    }

    if(wrong_way >= max_wrong_way){
      cat("L is inceasing, breaking and returning last U, V.\n")
      break
    }


    rmsd_u <- sqrt(mean((U - U_next)^2))
    rmsd_v <- sqrt(mean((V - V_next)^2))
    U <- U_next
    V <- V_next
    L_delta <- abs(L_next - L) / L
    L <- L_next
    if (L_delta > delta_target) {
      delta_met <- 0
    } else {
      delta_met <- delta_met + 1
    }

    if (verbose) {
      # cat(iter, L, rho, constraint_resid_norm, non_constraint_resid_norm, dual_resid_norm, rmsd_u, rmsd_v, max_V, n_bad, bad_sum, max_U, n_bad_U, bad_sum_U,
      # L_delta, Sys.time() - start_time, "\n")
      cat(iter, L, lagrangian[1], lagrangian[1] + lagrangian[2], delta_sign, rho, constraint_resid_norm, non_constraint_resid_norm, dual_resid_norm, bad_sum, bad_sum_U, bad_sum_G, L1V, L2U, DU, L_delta, Sys.time() - start_time,"\n")
    }

    # Converged if some iterations in a row have delta < target.
    if (delta_met >= min_delta_met) {
      break
    }
  }
  end_time <- Sys.time()
  return(list(
    "V" = V, "U" = U, "theta" = theta, "L" = L, "F_term" = lagrangian[1],
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
#' @param warm_start Boolean, TRUE to start next fit with result of previous
#'   fit. Default FALSE.
#' @param constraint One of "UV" or "VU". Constraint to use.
#' @param DAG Bool. True to resitrict solutions to approximate DAGs.
fit_inspre_sequence <- function(X, lambda, W = NULL, rho = 1.0,
                                its = 100, delta_target = 1e-4,
                                symmetrize = FALSE, verbose = 1,
                                gamma = 0.0, train_prop = 1.0,
                                mu = 5, tau = 2.5, solve_its = 3, ncores = 1,
                                warm_start = FALSE, constraint = "UV", DAG = FALSE) {
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
      solve_its = solve_its, ncores = ncores, constraint = constraint, DAG = DAG)
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
#' @param warm_start Boolean, TRUE to start next fit with result of previous
#'   fit. Default FALSE.
#' @param constraint One of "UV" or "VU". Constraint to use.
#' @param DAG Bool. True to resitrict solutions to approximate DAGs.
#' @export
inspre <- function(X, W = NULL, rho = 10.0, lambda = NULL,
                   lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                   gamma = NULL, its = 100, delta_target = 1e-4,
                   symmetrize = FALSE, verbose = 1, train_prop = 0.8,
                   cv_folds = 0, mu = 10, tau = 2, solve_its = 3, ncores = 1,
                   warm_start = TRUE, min_nz = 1e-5, constraint = "UV", DAG = FALSE) {
  D <- ncol(X)

  if(any(is.na(X))){
    if(is.null(W)){
      W <- matrix(1L, nrow = D, ncol = D)
    }
    W[is.na(X)] <- 0
  }

  if (is.null(lambda)) {
    lambda_max <- max(abs(off_diagonal(X)), na.rm = TRUE)
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
        ncores = ncores, warm_start = warm_start, constraint = constraint, DAG = DAG)$L

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
    ncores = ncores, warm_start = warm_start, constraint = constraint, DAG = DAG)

  if (cv_folds > 0) {
    xi_mat <- xi_mat <- array(0, dim = c(D, D, length(lambda)))
    error_matrix <- array(dim = c(cv_folds, length(lambda)))
    for (i in 1:cv_folds){
      if (verbose){
        cat(sprintf("Cross-validation iteration %d.\n", i))
      }
      cv_res <- fit_inspre_sequence(
        X = X, W = W, rho = rho, lambda = lambda, delta_target = delta_target,
        symmetrize = symmetrize, verbose = verbose, gamma = gamma,
        train_prop = train_prop, its = its, mu = mu, tau = tau,
        solve_its = solve_its, ncores = ncores, warm_start = warm_start, constraint = constraint, DAG = DAG)
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
#' @param constraint One of "UV" or "VU". Constraint to use.
#' @param DAG Bool. True to resitrict solutions to approximate DAGs.
#' @export
fit_inspre_from_R <- function(R_tce, W = NULL, rho = 100.0, lambda = NULL,
                              lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                              gamma = NULL, its = 100, delta_target = 1e-4,
                              verbose = 1, train_prop = 0.8,
                              cv_folds = 0, mu = 10, tau = 1.5, solve_its = 10,
                              ncores = 1, warm_start = FALSE, min_nz = 1e-5, constraint = "UV", DAG = FALSE){
  D <- dim(R_tce)[1]
  inspre_res <- inspre::inspre(
    X = R_tce, W = W, rho = rho, lambda = lambda,
    lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
    gamma = gamma, its = its, delta_target = delta_target, symmetrize = FALSE,
    verbose = verbose, train_prop = train_prop, cv_folds = cv_folds, mu = mu,
    tau = tau, solve_its = solve_its, ncores = ncores, warm_start = warm_start, min_nz = min_nz, constraint = constraint, DAG = DAG)
  inspre_res$R_hat <- array(0L, dim = dim(inspre_res$V))
  for(i in 1:length(inspre_res$lambda)){
    inspre_res$R_hat[ , , i] <-
      diag(D) - inspre_res$V[ , , i] / diag(inspre_res$V[ , , i])
  }
  dimnames(inspre_res$R_hat) <- list(rownames(R_tce), colnames(R_tce), inspre_res$lambda)
  return(inspre_res)
}


multiple_iv_reg <- function(target, .X, .targets){
  inst_obs <- .targets == target
  control_obs <- .targets == "control"
  n_target <- sum(inst_obs)
  n_control <- sum(control_obs)
  n_total <- n_target + n_control
  X_target <- .X[inst_obs, ]
  X_control <- .X[control_obs, ]
  this_feature <- which(colnames(.X) == target)
  X_exp_target <- X_target[, this_feature]
  X_exp_control <- X_control[, this_feature]
  beta_inst_obs <- sum(X_exp_target)/n_target
  sums_target <- colSums(X_target)
  beta_hat <- sums_target/sums_target[this_feature]
  resid_sum <- colSums((X_target - outer(X_exp_target, beta_hat))**2) +
    colSums((X_control - outer(X_exp_control, beta_hat))**2)
  var_resid <- resid_sum/(n_total*n_target*beta_inst_obs**2)
  se_hat <- sqrt(var_resid)
  names(beta_hat) <- paste(names(beta_hat), "beta_hat", sep="_")
  names(se_hat) <- paste(names(se_hat), "se_hat", sep="_")
  return(as.data.frame(c(list(target = target, beta_obs = beta_inst_obs), as.list(c(beta_hat, se_hat)))))
}


multiple_iv_reg_dumb_row <- function(target, .X, .targets){
  this_feature <- which(colnames(.X) == target)
  # X <- cbind(.X[, this_feature], 1)
  # Y <- .X[, -this_feature]
  # Z <- cbind(.targets == target, 1)
  rows <- .targets %in% c(target, "control")
  X <- cbind(.X[rows, this_feature], 1)
  Y <- .X[rows, -this_feature]
  Z <- cbind(.targets == target, 1)[rows, ]

  X_hat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  beta_hat <- solve(crossprod(X_hat), crossprod(X_hat, Y))
  u_hat <- Y - X %*% beta_hat
  C <- colSums(u_hat**2)/nrow(X)
  C <- C*solve(crossprod(X_hat))[1,1]
  beta_se <- sqrt(C)

  beta_hat <- append(beta_hat[1,], 1, after=this_feature-1)
  beta_se <-append(beta_se, 0, after=this_feature-1)

  names(beta_hat) <- paste(paste0("V", 1:D), "beta_hat", sep="_")
  names(beta_se) <- paste(paste0("V", 1:D), "se_hat", sep="_")
  return(as.data.frame(c(list(target = target, beta_obs = 0), as.list(c(beta_hat, beta_se)))))
}

multiple_iv_reg_dumb <- function(target, .X, .targets){
  D = ncol(.X)
  this_feature <- which(colnames(.X) == target)
  X <- cbind(.X[, this_feature], 1)
  Y <- .X[, -this_feature]
  Z <- cbind(.targets == target, 1)

  beta_obs <- solve(crossprod(Z), crossprod(Z, X))
  X_hat <- Z %*% beta_obs
  beta_hat <- solve(crossprod(X_hat), crossprod(X_hat, Y))
  u_hat <- Y - X %*% beta_hat
  C <- colSums(u_hat**2)/nrow(X)
  C <- C*solve(crossprod(X_hat))[1,1]
  beta_se <- sqrt(C)

  beta_hat <- append(beta_hat[1,], 1, after=this_feature-1)
  beta_se <-append(beta_se, 0, after=this_feature-1)

  names(beta_hat) <- paste(paste0("V", 1:D), "beta_hat", sep="_")
  names(beta_se) <- paste(paste0("V", 1:D), "se_hat", sep="_")
  return(as.data.frame(c(list(target = target, beta_obs = beta_obs[1,1]), as.list(c(beta_hat, beta_se)))))
}


calc_inst_effect <- function(target, .X, .targets){
  inst_obs <- .targets == target
  control_obs <- .targets == "control"
  n_target <- sum(inst_obs)
  n_control <- sum(control_obs)
  n_total <- n_target + n_control
  this_feature <- which(colnames(.X) == target)
  X_target <- .X[inst_obs, ]
  X_control <- .X[control_obs, ]
  X_exp_target <- X_target[, this_feature]
  X_exp_control <- X_control[, this_feature]
  Z <- c(rep(1, n_target), rep(0, n_control))
  inst_effect <- cor(Z, c(X_exp_target, X_exp_control))
  inst_se <- sqrt((1 - inst_effect**2)/(n_total - 2))
  return(as.data.frame(list(target = target, inst_effect = inst_effect, inst_se = inst_se, n = n_total)))
}


#' Calculates model predictions given new data, with cross-validation errors.
#'
#' This function uses the mean((X - XG - ZB)**2) error formulation.
#'
#' @param res Inspre result. Result of running fit_inspre for a single or
#'   sequence of hyperparameter values.
#' @param X observations x features data matrix. Colnames represent feature
#'   names. Rownames also correspond to feature names indicating an intervention
#'   was applied to that feature in that sample, with "control" indicating
#'   no intervention.
#' @param beta Sequence of floats. Obsevered scale estimates of instrument
#'   effects.
#' @export
predict_inspre <- function(res, .X, .beta, .targets){
  # Ensure that the beta and X cols have the same entries in the same order.
  vars_to_use <- intersect(colnames(.X), names(.beta))
  .beta <- .beta[vars_to_use]
  .X <- .X[, vars_to_use]

  ZB <- outer(.targets, colnames(.X), "==")
  ZB <- t(t(ZB) * .beta)

  D <- nrow(res$R_hat)
  ImGi <- purrr::map(1:length(res$lambda), ~ solve(diag(D) - res$R_hat[,,.x]))
  ImGi <- array(do.call(cbind, ImGi), dim=c(dim(ImGi[[1]]), length(ImGi)))

  eps_hat_G <- vector("numeric", length(res$lambda))
  eps_hat_I <- vector("numeric", length(res$lambda))
  for(i in 1:length(res$lambda)){
    X_hat_G <- .X %*% res$R_hat[,,i] + ZB
    X_hat_I <- ZB %*% ImGi[,,i]
    eps_hat_G[i] <- mean((.X - X_hat_G)**2)
    eps_hat_I[i] <- mean((.X - X_hat_I)**2)
  }
  return(list("eps_hat_G" = eps_hat_G, "eps_hat_I" = eps_hat_I))
}


#' Fits inverse sparse regression model.
#'
#' See also inspre::inspre() for more details.
#'
#' @param X observations x features data matrix. Colnames represent feature
#'   names. Rownames also correspond to feature names indicating an intervention
#'   was applied to that feature in that sample, with "control" indicating
#'   no intervention.
#' @param targets sequence of strings of length total number of observations
#'   (rows in X). Entries are either "control" to indicate no intervention
#'   or the name of a column in`X` to indicate the intervened on variable.
#' @param weighted Boolean. TRUE to calculate weights from SEs and use them,
#'   FALSE for unweighted. Default TRUE.
#' @param max_med_ratio Float or NULL. Ignored if `weight=FALSE`. If
#'   `weight=TRUE`, this is the ratio of the maximum weight to median weight.
#'   `NULL` for no restriction on the ratio. This can be useful to set
#'   if you have some entries with very small standard error, to prevent the
#'   algorithm from focusing exclusively on the entries with very small SE.
#' @param filter Bool. True to filter the produced TCE matrix with `fitler_tce`.
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
#'   training in cross-validation. NOT USED.
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
#' @param constraint One of "UV" or "VU". Constraint to use.
#' @param DAG Bool. True to restrict solutions to approximate DAGs. Useful to
#'   set to TRUE if you are having convergence issues with `DAG=FALSE` as the
#'   more restricted model can be easier to fit.
#' @export
fit_inspre_from_X <- function(X, targets, weighted = TRUE, max_med_ratio = NULL, filter = TRUE,
                              rho = 100.0, lambda = NULL,
                              lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                              gamma = NULL, its = 100, delta_target = 1e-4,
                              verbose = 1, train_prop = NULL,
                              cv_folds = 0, mu = 10, tau = 1.5, solve_its = 10,
                              ncores = 1, warm_start = FALSE, min_nz = 1e-5, constraint = "UV", DAG = FALSE){
  D <- ncol(X)
  target_names <- unique(targets)
  keep_features <- intersect(targets, colnames(X))
  keep_targets <- c(keep_features, "control")
  keep_obs <- targets %in% keep_targets
  targets <- targets[keep_obs]
  X <- X[keep_obs, keep_features]

  if(verbose){cat("Fitting model with full dataset")}
  inst_effects <- purrr::map(keep_features, ~ multiple_iv_reg(.x, X, targets))  %>% purrr::list_rbind()
  beta_obs <- inst_effects$beta_obs
  names(beta_obs) <- inst_effects$target
  R_hat <- data.matrix(dplyr::select(inst_effects, ends_with('beta_hat')) %>% dplyr::rename_with(~sub('_beta_hat', '', .x)))
  SE_hat <- data.matrix(dplyr::select(inst_effects, ends_with('se_hat')) %>% dplyr::rename_with(~sub('_se_hat', '', .x)))
  rownames(R_hat) <- inst_effects$target
  rownames(SE_hat) <- inst_effects$target
  if(filter){
    filtered <- filter_tce(R_hat, SE_hat)
    R_hat <- filtered$R
    SE_hat <- filtered$SE
    keep_rows <- rownames(R_hat)
    keep_cols <- colnames(R_hat)
  }
  if(weighted){
    W <- inspre::make_weights(SE_hat, max_med_ratio = max_med_ratio)
  } else{
    W <- NULL
  }

  full_res <- fit_inspre_from_R(R_hat, W = W, rho = rho, lambda = lambda,
                                lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
                                gamma = gamma, its = its, delta_target = delta_target,
                                verbose = verbose, train_prop = 1,
                                cv_folds = 0, mu = mu, tau = tau, solve_its = solve_its,
                                ncores = ncores, warm_start = warm_start, min_nz = min_nz, constraint = constraint, DAG = DAG)

  if(cv_folds > 0){
    folds <- caret::createFolds(1:sum(keep_obs), k=cv_folds)
    eps_hat_G <- array(0.0, length(full_res$lambda))
    eps_hat_I <- array(0.0, length(full_res$lambda))
    xi_mat <- array(0, dim = c(D, D, length(full_res$lambda)))
    for(i in seq_along(folds)){
      fold <- folds[[i]]
      if (verbose){
        cat(sprintf("Cross-validation iteration %d.\n", i))
      }
      X_train <- X[-fold, ]
      targets_train <- targets[-fold]
      inst_effects <- purrr::map(keep_features, ~ multiple_iv_reg(.x, X_train, targets_train))  %>% purrr::list_rbind()
      beta_obs <- inst_effects$beta_obs
      names(beta_obs) <- inst_effects$target
      R_hat_cv <- data.matrix(dplyr::select(inst_effects, ends_with('beta_hat')) %>% dplyr::rename_with(~sub('_beta_hat', '', .x)))
      SE_hat_cv <- data.matrix(dplyr::select(inst_effects, ends_with('se_hat')) %>% dplyr::rename_with(~sub('_se_hat', '', .x)))
      rownames(R_hat_cv) <- inst_effects$target
      rownames(SE_hat_cv) <- inst_effects$target
      if(filter){
        R_hat_cv <- R_hat_cv[keep_rows, keep_cols]
        SE_hat_cv <- SE_hat_cv[keep_rows, keep_cols]
        filtered <- filter_tce(R_hat_cv, SE_hat_cv, max_nan_perc = 1)
        R_hat_cv <- filtered$R
        SE_hat_cv <- filtered$SE
      }
      if(weighted){
        W_cv <- inspre::make_weights(SE_hat_cv, max_med_ratio = max_med_ratio)
      } else{
        W_cv <- NULL
      }
      cv_res <- fit_inspre_from_R(R_hat_cv, W = W_cv, rho = rho, lambda = lambda,
                                  lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
                                  gamma = gamma, its = its, delta_target = delta_target,
                                  verbose = verbose, train_prop = 1,
                                  cv_folds = 0, mu = mu, tau = tau, solve_its = solve_its,
                                  ncores = ncores, warm_start = warm_start, min_nz = min_nz, constraint = constraint, DAG = DAG)
      V_nz <- abs(cv_res$V) > min_nz
      xi_mat <- xi_mat + V_nz
      X_test <- X[fold, ]
      targets_test <- targets[fold]
      eps <- predict_inspre(cv_res, X_test, beta_obs, targets_test)
      eps_hat_G <- eps_hat_G + eps$eps_hat_G
      eps_hat_I <- eps_hat_I + eps$eps_hat_I
    }
    xi_mat <- xi_mat/cv_folds
    D_hat <- 2 * xi_mat * (1-xi_mat)
    D_hat <- apply(D_hat, 3, mean)
    D_hat_se <- apply(
      xi_mat, 3, function(x){ stats::sd(x)/sqrt(length(x)) })
    full_res$D_hat <- D_hat
    full_res$D_hat_se <- D_hat_se
    full_res$xi_mat <- xi_mat

    eps_hat_G <- eps_hat_G/cv_folds
    eps_hat_I <- eps_hat_I/cv_folds
    full_res$eps_hat_G <- eps_hat_G
    full_res$eps_hat_I <- eps_hat_I
  }

  # Hack until I have time to change all these variable names.
  full_res$G_hat <- full_res$R_hat
  full_res$R_hat <- R_hat
  full_res$SE_hat <- SE_hat
  full_res$W <- W
  return(full_res)
}


#' Creates an igraph from CDE matrix.
#'
#' @param R_cde Matrix of causal effects.
#' @param min_edge_value Minimum edge strength for pruning.
#' @param max_edge_value Set edges above this number to this.
#' @export
make_igraph <- function(R_cde, min_edge_value = 0.01, max_edge_value = NULL, scale_factor = NULL){
  adj_matrix <- R_cde
  adj_matrix[abs(adj_matrix) < min_edge_value] = 0
  if(!is.null(max_edge_value)){
    adj_matrix[abs(adj_matrix) > max_edge_value] = max_edge_value
  }
  if(!is.null(scale_factor)){
    adj_matrix <- adj_matrix/scale_factor
  }
  zeros <- adj_matrix == 0
  adj_matrix <- -log(abs(adj_matrix))
  adj_matrix[zeros] = 0
  return(igraph::graph_from_adjacency_matrix(
    adj_matrix, mode = "directed", weighted = TRUE))
}
