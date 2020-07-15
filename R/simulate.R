#' Calculates metrics for model evaluation.
#'
#' @param X DxD matrix of predicted parameters.
#' @param X_true DxD matrix of true parameters
#' @param eps float. Absolute values below eps will be considered 0.
calc_metrics <- function(X, X_true, eps = 1e-10) {
  D <- ncol(X)
  X <- X[!diag(D)]
  X_true <- X_true[!diag(D)]
  X[abs(X) < eps] <- 0
  X_true[abs(X_true) < eps] <- 0

  rmse <- sqrt(mean((X - X_true)^2))
  mae <- mean(abs(X - X_true))

  sign_X <- sign(X)
  sign_Xt <- sign(X_true)
  TN <- sum(!(abs(sign_X) | abs(sign_Xt)))
  TS <- sum(sign_X == sign_Xt) - TN
  FN <- sum((1 - abs(sign_X)) & abs(sign_Xt))
  FS <- sum(sign_X != sign_Xt) - FN

  acc <- mean(sign_X == sign_Xt)
  N_pos <- sum(X_true > 0)
  N_neg <- sum(X_true < 0)
  N_zero <- sum(X_true == 0)
  weights <- sign_Xt
  weights[sign_Xt > 0] <- 1 / N_pos
  weights[sign_Xt < 0] <- 1 / N_neg
  weights[sign_Xt == 0] <- 1 / N_zero
  weight_acc <- sum((sign_X == sign_Xt) * weights) / sum(weights)
  precision = TS / (TS + FS)
  recall = TS / (TS + FN)

  return(list("precision" = precision, "recall" = recall,
              "F1" = 2*precision*recall/(precision + recall),
              "rmse" = rmse, "mae" = mae, "acc" = acc,
              "weight_acc" = weight_acc))
}

#' Calculates STARS cross-validation.
#'
#' @param X N x D data matrix.
#' @param method Function that returns a list with two elements. "theta", a
#'   D x D x nlambda matrix and "lambda" a list of lambda values used.
#' @param train_prop Float, proportion of rows of X to use at each iteration.
#' @param cv_folds Integer, number of cross-validation folds.
stars_cv <- function(X, method, train_prop = 0.8, cv_folds = 10){
  N = dim(X)[1]
  D = dim(X)[2]
  S_full <- cor_w_se(X)$S_hat
  method_res <- method(S_full)
  theta_hat <- method_res$theta
  lambda <- method_res$lambda
  xi_mat <- array(0, dim = c(D, D, length(lambda)))
  for (fold in 1:cv_folds){
    train <- stats::runif(N) < train_prop
    X_train = X[train, ]
    cor_res <- cor_w_se(X_train)
    S_train <- cor_res$S_hat
    theta_cv <- method(S_train)$theta
    theta_nz <- abs(theta_cv) > 1e-8
    xi_mat <- xi_mat + theta_nz
  }
  xi_mat <- xi_mat/cv_folds
  xi_mat <- 2 * xi_mat * (1-xi_mat)
  D_hat <- apply(xi_mat, 3, mean)
  return(list(
    "theta" = theta_hat,
    "lambda" = lambda,
    "D_hat" = D_hat))
}

#' Estimates correlation matrix of X in the presence of missing data.
#'
#' Uses simple approximation sqrt((1-r^2)/(n-2)) which is technically
#' incorrect but works well in practice.
#'
#' @param X NxD matrix with optionally missing entries.
cor_w_se <- function(X) {
  S_hat <- stats::cor(X, use = "pairwise.complete.obs")
  N <- apply(X, 2, function(x) {
    apply(X, 2, function(y) {
      sum(!is.na(x) & !is.na(y))
    })
  })
  SE_S <- sqrt((1 - S_hat**2) / (N - 2))
  return(list("S_hat" = S_hat, "SE_S" = SE_S, "N" = N))
}
