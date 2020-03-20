#' Calculates metrics for model evaluation.
#'
#' @param X DxD matrix of predicted parameters.
#' @param X_true DxD matrix of true parameters
#' @param eps float. Absolute values below eps will be considered 0.
calc_metrics <- function(X, X_true, eps = 1e-10){
  D <- ncol(X)
  X <- X[!diag(D)]
  X_true <- X_true[!diag(D)]
  X[abs(X) < eps] <- 0
  X_true[abs(X_true) < eps] <- 0

  rmse <- sqrt(mean((X - X_true)^2))
  mae <- median(abs(X-X_true))

  sign_X <- sign(X)
  sign_Xt <- sign(X_true)
  TN <- sum(!(abs(sign_X) | abs(sign_Xt)))
  TS <- sum(sign_X == sign_Xt) - TN
  FN <- sum((1-abs(sign_X)) & abs(sign_Xt))
  FS <- sum(sign_X != sign_Xt) - FN

  acc <- mean(sign_X == sign_Xt)
  N_pos <- sum(X_true > 0)
  N_neg <- sum(X_true < 0)
  N_zero <- sum(X_true == 0)
  weights <- sign_Xt
  weights[sign_Xt > 0] <- 1/N_pos
  weights[sign_Xt < 0] <- 1/N_neg
  weights[sign_Xt == 0] <- 1/N_zero
  weight_acc <- sum((sign_X == sign_Xt)*weights)/sum(weights)

  return(list("precision" = TS/(TS+FS), "recall" = TS/(TS+FN), "rmse" = rmse, "mae" = mae, "acc" = acc, "weight_acc" = weight_acc))
}

#' Generates correlation/precision matrices for gaussian graphical model.
#'
#' @param D Integer. Dimensionality.
#' @param p_net Float. Connection density.
#' @param sd_net Float. Standard deviation of initial covariance.
#' @param normalize Bool. True to return data and parameters with mean 0,
#'   var 1.
generate_dataset <- function(D, N, p_net, sd_net = 1.0, normalize = TRUE,
                             min_missing = 0.0, max_missing = 0.0){
  B <- matrix(rnorm(D*D, sd=sd_net), nrow=D)
  theta <- 0.5*(B + t(B))
  zeros <- stats::runif(D*(D-1)/2) > p_net
  theta[lower.tri(theta)][zeros] <- 0
  theta <- t(theta)
  theta[lower.tri(theta)][zeros] <- 0
  min_ev <- eigen(theta, only.values = TRUE)$values[D]
  theta <- theta + (1-min_ev)*diag(D)
  sigma <- solve(theta)

  X <- mvtnorm::rmvnorm(N, sigma = sigma)
  if(max_missing > 0){
    missing <- stats::runif(D, min = min_missing, max = max_missing)
    drop <- t(matrix(stats::runif(N*D) < missing, nrow = D))
    X[drop] <- NA
  }

  if(normalize){
    sigma <- stats::cov2cor(sigma)
    theta <- solve(sigma)
    X <- scale(X)
  }

  return(list("X" = X, "sigma" = sigma, "theta" = theta))
}

#' Estimates correlation matrix of X in the presence of missing data.
#'
#' Uses simple approximation sqrt((1-r^2)/(n-2)) which is technically
#' incorrect but works well in practice.
cor_w_se <- function(X){
  S_hat <- cor(X, use = "pairwise.complete.obs")
  N <- apply(X, 2, function(x){apply(X, 2, function(y){sum(!is.na(x) & !is.na(y))})})
  SE_S <- sqrt((1-S_hat**2)/(N-2))
  return(list("S_hat" = S_hat, "SE_S" = SE_S, "N" = N))
}
