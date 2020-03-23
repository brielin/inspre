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

inspre_primal <- function(X, U, V, lambda, gamma) {
  VU_minus_I <- V %*% U - diag(nrow(U))
  x_dist <- .5 * sum((X - U)^2, na.rm = TRUE)
  l1 <- lambda * sum(abs(off_diagonal(V)))
  u_shrinkage <- gamma * sum((U - diag(nrow(U)))^2)
  nnz <- sum(V > 0)
  return(list("loss" = c(x_dist, u_shrinkage, l1), "nnz" = nnz))
}

inspre_lagrangian <- function(X, U, V, theta, rho, lambda, gamma) {
  VU_minus_I <- V %*% U - diag(nrow(U))
  VU_term1 <- sum(theta * VU_minus_I)
  VU_term2 <- .5 * rho * sum(VU_minus_I^2)
  return(c(inspre_primal(X, U, V, lambda, gamma)$loss, VU_term1, VU_term2))
}

#' Helper function to make weights for inspre.
#'
#' @param SE DxD matrix of standard errors.
#' @param max_weight Highest possible weight, useful when some SEs are
#'   near 0.
#' @export
make_weights <- function(SE, max_weight = NULL) {
  weights <- 1 / SE^2
  weights[is.na(SE)] <- 0

  if (is.null(max_weight)) {
    infs <- is.infinite(weights)
    weights[infs] <- 0
    max_weight <- max(weights)
    weights[infs] <- max_weight
  }

  weights[weights > max_weight] <- max_weight
  weights <- weights / mean(weights)
  return(weights)
}

#' Worker function to fit inspre for a single value of lambda.
#'
#' TODO(brielin): Remove gamma penalty.
#' TODO(brielin): Infer detpen automatically.
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param W DxD Matrix of weights.
#' @param rho Float, learning rate for ADMM.
#' @param lambda Float. L1 regularization strength on inverse of X.
#' @param gamma Float. L2 regularization on approximation to X that is inverted.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param warm_start List with elements "U" and "V" containing starting point.
#' @param symmetrize True to force the output to be symmetric. If the input
#'   is symmetric, the output isn't always perfectly symmetric due to numerical
#'   issues.
#' @param verbose Bool. True to print solution progress.
#' @param detpen Determinant penalty strength.
inspre_worker <- function(X, W = NULL, rho = 0.1, lambda = 0.01, gamma = 0.0,
                          its = 1000, delta_target = 1e-4, warm_start = NULL,
                          symmetrize = FALSE, verbose = FALSE, detpen = 0.1) {
  D <- nrow(X)
  theta <- matrix(0.0, D, D)
  if (is.null(warm_start)) {
    V <- diag(D)
    U <- X
    U[is.na(U)] <- 0
  } else {
    V <- warm_start$V
    U <- warm_start$U
  }

  L <- sum(inspre_lagrangian(X, U, V, theta, rho, lambda, gamma))
  converged <- TRUE
  sign_changes <- 0
  max_sign_changes <- 2
  last_sign <- 1
  min_delta_met <- 1
  delta_met <- 0
  for (iter in 1:its) {
    if (is.null(W)) {
      U_next <- solve(
        (1 + gamma) * diag(D) + rho * t(V) %*% V,
        X + rho * t(V) - t(V) %*% theta + gamma * diag(D)
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
    # if(symmetrize){
    #   U_next <- (0.5 * (U_next + t(U_next)))
    # }

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

    VU_minus_I <- V_next %*% U_next - diag(D)
    theta <- theta + rho * VU_minus_I

    rmsd_u <- sqrt(mean((U - U_next)^2))
    rmsd_v <- sqrt(mean((V - V_next)^2))
    U <- U_next
    V <- V_next

    L_next <- sum(inspre_lagrangian(X, U, V, theta, rho, lambda, gamma))
    L_delta <- abs(L_next - L) / L
    L_sign <- sign(L_next - L)
    L <- L_next

    if (L_delta > delta_target) {
      delta_met <- 0
    } else {
      delta_met <- delta_met + 1
    }
    if ((L_sign != last_sign) & (L_sign == 1)) {
      sign_changes <- sign_changes + 1
    }

    if (verbose) {
      vsvs <- svd(V, nu = 0, nv = 0)$d
      usvs <- svd(U, nu = 0, nv = 0)$d
      cat(
        iter, L, max(vsvs) / min(vsvs), max(usvs) / min(usvs),
        sqrt(mean(VU_minus_I^2)), rmsd_u, rmsd_v, L_delta, sign_changes, "\n"
      )
    }

    # Converged if some iterations in a row have delta < target.
    if (delta_met >= min_delta_met) {
      break
    }
    # Diverged if L switches from decreasing to increasing too many times.
    if (sign_changes >= max_sign_changes) {
      converged <- FALSE
      break
    }
    last_sign <- L_sign
  }
  if (isTRUE(converged) & last_sign == 1) {
    converged <- FALSE
  }
  return(list(
    "V" = V, "U" = U, "converged" = converged, "L" = L, "L_delta" = L_delta,
    "iter" = iter
  ))
}

#' Finds the inverse of X using Inverse Sparse Regression.
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param W DxD Matrix of weights.
#' @param rho Float, sequence of floats or NULL. Learning rate for ADMM.
#'   If NULL, a logarithmicallly spaced set of values between rho_min and
#'   rho_max will be used.
#' @param lambda Float, sequence of floats of NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmicallly spaced set of values between
#'   the maximimum absolute off diagonal element of X and lambda_min_ratio
#'   times this value will be used.
#' @param lambda_min_ratio Float, ratio of maximum lambda to minimum lambda.
#' @param nlambda Integer. Number of lambda values to try.
#' @param alpha Float. To be removed in future version.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param rho_min Float. Minimum rho to try.
#' @param rho_max Float. Maximum rho to try.
#' @param symmetrize True to force the output to be symmetric. If the input
#'   is symmetric, the output isn't always perfectly symmetric due to numerical
#'   issues.
#' @param verbose 0, 1 or 2. 2 to print convergence progress for each lambda,
#'   1 to print convergence result for each lambda, 0 for no output.
#' @param detpen Determinant penalty strength.
#' @export
inspre <- function(X, W = NULL, rho = NULL, lambda = NULL,
                   lambda_min_ratio = 1e-2, nlambda = 20, alpha = 1.0,
                   its = 1000, delta_target = 1e-4, rho_min = 0.1,
                   rho_max = 100, symmetrize = FALSE, verbose = 1,
                   detpen = 0.0) {
  D <- ncol(X)
  if (is.null(lambda)) {
    lambda_max <- max(abs(off_diagonal(X)), na.rm = TRUE)
    lambda_min <- lambda_min_ratio * lambda_max
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
    lambda <- alpha * lambda
    gamma <- (1 - alpha) * lambda
  }

  if (is.null(rho)) {
    rho <- exp(seq(log(rho_min), log(rho_max), length.out = 10))
  }

  warm_start <- NULL
  V_all <- array(dim = c(D, D, length(lambda)))
  j <- 1
  rho_used <- vector("numeric", length = length(lambda))
  L_path <- vector("numeric", length = length(lambda))
  for (i in 1:length(lambda)) {
    lambda_i <- lambda[i]
    gamma_i <- gamma[i]
    restart <- FALSE
    while (j <= length(rho)) {
      rho_j <- rho[j]
      inspre_res <- inspre_worker(X, W,
        lambda = lambda_i, gamma = gamma_i,
        rho = rho_j, warm_start = warm_start, its = its,
        delta_target, symmetrize = symmetrize,
        verbose = (verbose == 2), detpen = detpen
      )
      if (isFALSE(inspre_res$converged)) {
        if (verbose) {
          cat(sprintf(
            "Diverged in %d iterations at rho=%f for lambda=%f.
            Trying rho=%f.\n",
            inspre_res$iter, rho_j, lambda_i, rho[j + 1]
          ))
        }
        j <- j + 1
        rho_j <- rho[j]
        if (j == length(rho) & isFALSE(restart)) {
          restart <- TRUE
          j <- 1
        }
      } else {
        if (verbose) {
          cat(sprintf(
            "Converged to L=%f in %d iterations for lambda=%f using rho=%f.\n",
            inspre_res$L, inspre_res$iter, lambda_i, rho_j
          ))
        }
        if (j > 1) {
          j <- 1
        }
        break
      }
    }
    if (j == length(rho) + 1) {
      stop(sprintf("Error: did not converge for lambda=%f.\n", lambda_i))
    }
    warm_start <- inspre_res
    V_all[, , i] <- inspre_res$V
    rho_used[i] <- rho_j
    L_path[i] <- inspre_res$L
  }
  return(list("V" = V_all, "lambda" = lambda, "rho" = rho_used, "L" = L_path))
}

#' Fits simpler inverse LASSO model for one lambda value.
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param lambda Float, sequence of floats of NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmicallly spaced set of values between
#'   the maximimum absolute off diagonal element of X and lambda_min_ratio
#'   times this value will be used.
#' @param lambda_min_ratio Float, ratio of maximum lambda to minimum lambda.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param warm_start DxD matrix, initial solution.
#' @param verbose Bool. True to print solution progress.
ilasso_worker <- function(X, lambda, lambda_min_ratio = 1e-3, its = 1000,
                          delta_target = 1e-4, warm_start = NULL,
                          verbose = FALSE) {
  D <- nrow(X)
  resp <- diag(D)
  if (is.null(warm_start)) {
    V <- matrix(0.0, nrow = D, ncol = D)
  } else {
    V <- warm_start
  }

  # TODO(brielin): It's probably a lot faster to move this into the Cpp code.
  L_delta <- Inf
  L <- 0
  for (d in 1:D) {
    penalty_factor <- rep(1, D)
    penalty_factor[d] <- 0
    V_row <- V[d, ]
    for (iter in 1:its) {
      row_old <- V_row + 0.0
      L_next <- lasso_one_iteration(
        X, resp[, d], V_row, lambda = penalty_factor * lambda)
      v_delta <- sqrt(mean((V_row - row_old)**2))
      if (iter > 1) {
        L_delta <- (L - L_next) / L
      }
      L <- L_next
      if (verbose) {
        cat(d, iter, lambda, L, v_delta, L_delta, "\n")
      }
      if (L_delta < delta_target) {
        break
      }
    }
    V[d, ] <- V_row
  }
  return(V)
}

#' Fits simpler inverse LASSO model for one lambda value.
#'
#' @param X DxD Matrix to find approximate sparse inverse of.
#' @param lambda Float, sequence of floats or NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmically spaced sequence of lambda
#'
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param warm_start DxD matrix, initial solution.
#' @param verbose Bool. True to print solution progress.
#' @export
ilasso <- function(X, lambda = NULL, lambda_min_ratio = 1e-2, its = 1000,
                   delta_target = 1e-4, verbose = FALSE) {
  D <- ncol(X)
  if (is.null(lambda)) {
    lambda_max <- max(abs(off_diagonal(X)))
    lambda_min <- lambda_min_ratio * lambda_max
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = 20))
  }

  V_hat <- NULL
  V_all <- array(dim = c(D, D, length(lambda)))
  for (i in 1:length(lambda)) {
    lambda_i <- lambda[i]
    V_hat <- ilasso_worker(
      X, lambda_i, its = its, delta_target = delta_target,
      warm_start = V_hat, verbose = verbose)
    V_all[, , i] <- V_hat
  }
  return(list("V" = V_all, "lambda" = lambda))
}
