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

ilasso_loss <- function(X, V, lambda) {
  iloss <- sum((diag(nrow(X)) - X %*% V)^2)
  l1 <- lambda * sum(abs(off_diagonal(V)))
  nnz <- sum(V > 0)
  return(list("loss" = c(iloss, l1), "nnz" = nnz))
}

glasso_loss <- function(X, V, lambda) {
  I_hat <- X %*% V
  trloss <- sum(diag(I_hat))
  detV <- determinant(V)
  neglogdet <- -detV$sign * detV$modulus
  l1 <- lambda * sum(abs(off_diagonal(V)))
  nnz <- sum(V > 0)
  return(list("loss" = c(trloss, neglogdet, l1), "nnz" = nnz))
}

inspre_primal <- function(X, U, V, lambda, alpha) {
  x_dist <- .5 * sum((X - U)^2, na.rm = TRUE)
  l1_term <- lambda * sum(abs(off_diagonal(V)))
  det_term <- -1 * alpha * determinant(V)$modulus
  return(c(x_dist, l1_term, det_term))
}

inspre_lagrangian <- function(X, U, V, theta, rho, lambda, alpha) {
  VU_minus_I <- V %*% U - diag(nrow(U))
  VU_term1 <- sum(theta * VU_minus_I)
  VU_term2 <- .5 * rho * sum(VU_minus_I^2)
  return(c(inspre_primal(X, U, V, lambda, alpha), VU_term1, VU_term2))
}

#' Helper function to make weights for inspre.
#'
#' @param SE DxD matrix of standard errors.
#' @param max_min_ratio Float > 1. Ratio of maximum weight to minimum non-zero
#'   weight. Improves conditioning when some SEs are very small.
#' @export
make_weights <- function(SE, max_min_ratio = NULL) {
  weights <- 1 / SE^2
  weights[is.na(SE)] <- 0

  if (is.null(max_min_ratio)) {
    infs <- is.infinite(weights)
    weights[infs] <- 0
    max_weight <- max(weights)
    weights[infs] <- max_weight
  } else {
    max_weight <- min(weights[weights > 0]) * max_min_ratio
    weights[weights > max_weight] <- max_weight
  }

  weights <- weights / mean(weights)
  return(weights)
}

#' Worker function to fit inspre for a single value of lambda.
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
#' @param detpen Determinant penalty strength.
#' @param mu rho modification parameter for ADMM. Rho will be
#'   increased/decreased when the dual constrant and primal constraint are off
#'   by a factor of > mu.
#' @param tau rho modification parameter for ADMM. When called for, rho will be
#'   increased/decreased by the factor tau.
inspre_worker <- function(X, W = NULL, rho = 1.0, lambda = 0.01,
                          its = 1000, delta_target = 1e-4, warm_start = NULL,
                          symmetrize = FALSE, verbose = FALSE, detpen = 0.0,
                          mu = 5, tau = 2) {
  D <- nrow(X)
  theta <- matrix(0.0, D, D)
  if (is.null(warm_start)) {
    V <- diag(D)
    U <- diag(D)
    U[is.na(U)] <- 0
  } else {
    V <- warm_start$V
    U <- warm_start$U
  }
  L <- sum(inspre_lagrangian(X, U, V, theta, rho, lambda, alpha = detpen))

  min_delta_met <- 2
  delta_met <- 0
  for (iter in 1:its) {
    if (is.null(W)) {
      U_next <- qr.solve(
        diag(D) + rho * t(V) %*% V,
        X + rho * t(V) - t(V) %*% theta
      )
    } else {
      U_next <- matrix(0, D, D)
      WX <- W * X
      WX[is.na(WX)] <- 0
      rhs <- WX + rho * t(V) - t(V) %*% theta
      for (d in 1:D) {
        U_next[, d] <- solve(diag(W[, d]) + rho * t(V) %*% V, rhs[, d])
      }
    }

    glm_X <- sqrt(rho) * t(U_next)
    glm_Y <- ((rho + detpen) * diag(D) - t(theta)) / sqrt(rho)

    V_next <- matrix(0, D, D)
    # TODO(brielin): It's probably a lot faster to move this into the Cpp code.
    for (d in 1:D) {
      penalty_factor <- rep(1, D)
      penalty_factor[d] <- 0
      V_row <- V[d, ]
      lasso_one_iteration(
        glm_X, glm_Y[, d], V_row, lambda = penalty_factor * lambda)
      V_next[d, ] <- V_row
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
    lagrangian <- inspre_lagrangian(X, U, V, theta, rho, lambda, alpha = detpen)
    L_next <- sum(lagrangian)
    L_delta <- abs(L_next - L) / L
    L <- L_next

    if (L_delta > delta_target) {
      delta_met <- 0
    } else {
      delta_met <- delta_met + 1
    }

    if (verbose) {
      cat(
        iter, L, rho, constraint_resid_norm, dual_resid_norm, rmsd_u, rmsd_v, L_delta, "\n"
      )
    }

    # Converged if some iterations in a row have delta < target.
    if (delta_met >= min_delta_met) {
      break
    }
  }
  return(list(
    "V" = V, "U" = U, "L" = L, "F_term" = lagrangian[1],
    "l1_term" = lagrangian[2], "det_term" = lagrangian[3], "rho" = rho, "L_delta" = L_delta,
    "constraint_resid" = constraint_resid_norm,
    "dual_resid" = dual_resid_norm, "iter" = iter
  ))
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
#' @param detpen Determinant penalty strength.
#' @param train_prop Proportion of the dataset to use to train the model. Set
#'   to 1.0 to use all data.
fit_inspre_sequence <- function(X, lambda, W = NULL, rho = 1.0,
                                    its = 1000, delta_target = 1e-4,
                                    symmetrize = FALSE, verbose = 1,
                                    detpen = 0.0, train_prop = 1.0, warm_start = NULL) {
  # Break the matrix into training and test sets, equally and at random.
  D <- nrow(X)
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
  F_term_path <-
  train_error <- vector("numeric", length = length(lambda))
  test_error <- vector("numeric", length = length(lambda))
  for (i in 1:length(lambda)) {
    lambda_i <- lambda[i]
    if(length(detpen) > 1){
      detpen_i = detpen[i]
    } else {
      detpen_i = detpen
    }
    inspre_res <- inspre_worker(X = X, W = train_W,
                                lambda = lambda_i,
                                rho = rho, warm_start = warm_start, its = its,
                                delta_target = delta_target, symmetrize = symmetrize,
                                verbose = (verbose == 2), detpen = detpen_i)
    if (verbose) {
      cat(sprintf("Converged to L = %f in %d iterations for lambda = %f.\n L_delta = %e, constraint_resid = %e, dual_resid = %e.\n",
                  inspre_res$L, inspre_res$iter, lambda_i, inspre_res$L_delta, inspre_res$constraint_resid, inspre_res$dual_resid))
    }
    warm_start <- inspre_res
    V_all[, , i] <- inspre_res$V
    U_all[, , i] <- inspre_res$U
    rho_used[i] <- inspre_res$rho
    L_path[i] <- inspre_res$L
    if(is.null(train_W)){
      train_error[i] <- sqrt(mean(off_diagonal(X - inspre_res$U)^2, na.rm = TRUE))
    } else{
      train_error[i] <- sqrt(mean(off_diagonal(train_W * (X - inspre_res$U))^2, na.rm = TRUE))
    }
    test_error[i] <- sqrt(mean(off_diagonal(test_W * (X - inspre_res$U))^2, na.rm = TRUE))
  }
  return(list("V" = V_all, "U" = U_all, "lambda" = lambda, "rho" = rho_used, "L" = L_path,
              train_error = train_error, test_error = test_error))
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
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param symmetrize True to force the output to be symmetric. If the input
#'   is symmetric, the output isn't always perfectly symmetric due to numerical
#'   issues.
#' @param verbose 0, 1 or 2. 2 to print convergence progress for each lambda,
#'   1 to print convergence result for each lambda, 0 for no output.
#' @param detpen Determinant penalty strength.
#' @param train_prop Float between 0 and 1. Proportion of data to use for
#'   training in cross-validation.
#' @param cv_folds Integer. Number of cross-validation folds to perform.
#' @export
inspre <- function(X, W = NULL, rho = 1.0, lambda = NULL,
                   lambda_min_ratio = 1e-2, nlambda = 20,
                   its = 1000, delta_target = 1e-4,symmetrize = FALSE, verbose = 1,
                   detpen = 0.0, train_prop = 0.8, cv_folds = 0) {
  D <- ncol(X)
  if (is.null(lambda)) {
    lambda_max <- max(abs(off_diagonal(X)), na.rm = TRUE)
    lambda_min <- lambda_min_ratio * lambda_max
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  }

  if (verbose){
    cat("Fitting model with full dataset.\n")
  }
  full_res <- fit_inspre_sequence(
    X = X, W = W, rho = rho, lambda = lambda, delta_target = delta_target,
    symmetrize = symmetrize, verbose = verbose, detpen = detpen, its = its,
    train_prop = 1.0)

  if (cv_folds > 0) {
    xi_mat <- array(0, dim = c(D, D, length(lambda)))
    error_matrix <- array(dim = c(cv_folds, length(lambda)))
    for (i in 1:cv_folds){
      if (verbose){
        cat(sprintf("Cross-validation iteration %d.\n", i))
      }
      cv_res <- fit_inspre_sequence(
        X = X, W = W, rho = rho, lambda = lambda, delta_target = delta_target,
        symmetrize = symmetrize, verbose = verbose, detpen = detpen,
        train_prop = train_prop, its = its)
      error_matrix[i, ] = cv_res$test_error
      V_nz <- abs(cv_res$V) > 1e-8
      xi_mat <- xi_mat + V_nz
    }
    xi_mat <- xi_mat/cv_folds
    xi_mat <- 2 * xi_mat * (1-xi_mat)
    D_hat <- apply(xi_mat, 3, mean)
    D_hat_se <- apply(xi_mat, 3, function(x){ stats::sd(x)/sqrt(length(x)) })

    test_error <- apply(error_matrix, 2, mean)
    test_error_se <- apply(error_matrix, 2, function(x){ stats::sd(x)/sqrt(length(x))})
    full_res$test_error_se <- test_error_se
    full_res$test_error <- test_error
    full_res$D_hat <- D_hat
    full_res$D_hat_se <- D_hat_se
  }
  return(full_res)
}
