require(igraph)
require(mc2d)


scale_free <- function(D, alpha=0.1, beta=0.4, gamma=0.5, a_in = 1, a_out = 1, n = 1){
  A = matrix(c(0, 0, 1, 0), nrow = 2)
  while(nrow(A) < D){
    d_in <- colSums(A)
    d_out <- rowSums(A)
    choice <- stats::runif(1)
    D_curr <- nrow(A)
    if(choice < alpha){ # Edge from new to existing
      w <- sample(D_curr, 1, prob = d_in + a_in)
      a_next <- rep(0, D_curr)
      a_next[w] <- 1
      A <- cbind(rbind(A, a_next, deparse.level = 0), rep(0, D_curr + 1))
    } else if(choice < alpha + beta){ # Edge between existing
      v <- sample(D_curr, 1, prob = d_out + a_out)
      w <- sample(D_curr, 1, prob = d_in + a_in)
      if(v != w) A[v, w] = 1
    } else{ # Edge from existing to new
      v <- sample(D_curr, 1, prob = d_out + a_out)
      a_next <- rep(0, D_curr)
      a_next[v] <- 1
      A <- rbind(cbind(A, a_next, deparse.level = 0), rep(0, D_curr + 1))
    }
    if(n > 1){
      for(i in 2:n){
        v <- sample(D_curr, 1, prob = d_out + a_out)
        w <- sample(D_curr, 1, prob = d_in + a_in)
        if(v != w) A[v, w] = 1
      }
    }
  }
  return(A)
}


random <- function(D, p = 3/D){
  A <- matrix(as.integer(runif(10*10) < 0.3), nrow=10)
  diag(A) <- 0
  return(A)
}


#' Simulates random graph, optionally a DAG
#'
#' Edges are sampled from a PERT distribution with random sign. By
#' default these are between 0.05 and 0.8 with a mode of 0.2 and
#' intended to represent "root variance-explained".
#'
#' @param D Integer. Number of nodes to simulate.
#' @param graph String. One of 'scale-free' or 'random'. Type of network
#'   to simulate.
#' @param v_mode Float. Mode of the pert distribution.
#' @param v_min Float. Min of the pert distribution.
#' @param v_max Float. Max of the pert distribuion.
#' @param DAG Bool. TRUE to ensure the returned graph is a DAG.
#' @export
generate_network <- function(D, graph = 'scale-free', v_mode = 0.2,
                             v_max = 0.8, v_min = 0.05, DAG = FALSE){
  if(graph == 'scale-free'){
    A <- scale_free(D, n = 3*(DAG + 1))
  } else if(graph == 'random'){
    A <- random(D, p = 3/D * (DAG+1))
  } else{
    stop("NotImplementedException.")
  }
  if(DAG){
    perm <- sample.int(D)
    A <- A[perm, perm]
    A <- A * upper.tri(A)
  }
  G = A
  G[A != 0] <- sample(c(1, -1), sum(A), replace = T, prob = c(0.6, 0.4)) *
    mc2d::rpert(sum(A), min = v_min, mode = v_mode, max = v_max)
  rownames(G) <- paste0("V", 1:D)
  colnames(G) <- paste0("V", 1:D)
  return(G)
}


#' Turns a fully observed network into a partially observed network.
#'
#' This turns a fully observed network matrix into a partially observed one
#' where TCEs of missing nodes are integrated into pseudo-DCEs for the observed
#' nodes.
#'
#' @param G D x D matrix representing direct causal graph.
#' @param p Float 0 to 1, proportion of nodes to hide.
#' @export
censor_network <- function(G, p){
  R <- get_tce(get_observed(G))
  D <- nrow(R)
  keep_cols <- sort(sample.int(D, size = round((1-p)*D)))
  R_cens <- R[keep_cols, keep_cols]
  G_cens <- fit_exact(R_cens)$R_hat
  return(G_cens)
}


generate_data_inspre <- function(G, N_cont, N_int, size=N_int/10, min_int=N_int/10, int_r2=0.1,
                                 int_dir='negative', noise='gaussian'){
  D <- nrow(G)
  r2 <- mc2d::rpert(D, min = 0.05, mean = int_r2, max = 0.25)
  int_sizes <- rnbinom(D, mu=N_int - min_int, size=size) + min_int

  beta <- sqrt(((int_sizes + N_cont)/int_sizes)*r2)
  if(int_dir == 'negative'){
    beta <- -beta
  } else if(int_dir == 'both'){
    beta <- (mc2d::rbern(D)*2 - 1)*beta
  } else{
    stop('NotImplementedException')
  }

  Ncs <- cumsum(c(N_cont, int_sizes))
  N <- Ncs[length(Ncs)]
  XB <- matrix(0, nrow=sum(N), ncol=D)

  for(d in 1:D){
    start = Ncs[d]
    end = Ncs[d+1]
    XB[(start+1):end, d] = 1
  }
  XB <- t(t(XB) * beta)

  net_vars <- colSums(G**2)
  inst_vars <- apply(XB, 2, var)
  eps_vars <- 1 - net_vars - inst_vars
  if(any(eps_vars < 0)){
    stop("Noise Variance cannot be negative.")
  }
  if(noise == 'gaussian'){
    eps <- t(matrix(rnorm(D*N, sd=sqrt(eps_vars)), nrow=D, ncol=N))
  } else{
    stop('NotImplementedError')
  }
  Y <- (XB + eps) %*% solve(diag(D) - G)
  # This is required in order for the inst regression to have intercept term 0.
  mu_cont <- colMeans(Y[1:N_cont, ])
  Y <- t((t(Y) - mu_cont))

  colnames(Y) <- paste0("V", 1:D)
  targets <- c(rep("control", N_cont), paste0("V", rep(1:D, times=int_sizes)))
  return(list(Y = Y, targets = targets))
}


#' Simulates an intervention dataset.
#'
#' Simulates a graph-intervention dataset with corresponding graph.
#'
#' @param D Integer. Number of nodes to simulate.
#' @param N_cont Integer. Number of control samples to simulate.
#' @param N_int Integer. Mean of number of intervention samples to simulate per node.
#' @param size Float. Size parameter for NB distribution for int samples.
#' @param int_r2 Float. Mean variance in node explained by intervention.
#' @param int_dir string. 'positive', 'negative', or 'both'. Direction of
#'  effect of intervention on node.
#' @param graph String. One of 'scale-free' or 'random'. Type of network
#'   to simulate.
#' @param v_mode Float. Mode of the pert distribution.
#' @param v_min Float. Min of the pert distribution.
#' @param v_max Float. Max of the pert distribuion.
#' @param DAG Bool. TRUE to ensure the returned graph is a DAG.
#' @param confoudning Float >= 0. Proportion of nodes to "hide" to simulate confounding.
#' @param noise String. Noise model to simulate, currently just "gaussian".
#' @param model Data generating model. One of "inspre" or "dotears".
generate_dataset <- function(D, N_cont, N_int, size=N_int/10, min_int=N_int/10, int_r2=0.1, int_dir='negative',
                             graph = 'scale-free', v_mode = 0.2, v_max = 0.8, v_min = 0.05,
                             DAG = FALSE, confounding = 0, noise = 'gaussian', model = 'inspre'){
  G <- generate_network(D, graph, v_mode, v_max, v_min, DAG)
  if(confounding > 0){
    G <- censor_network(G, confounding)
  }
  data = generate_data_inspre(G, N_cont, N_int, size, min_int, int_r2, int_dir, noise)
  return(list(Y=data$Y, G=G, targets=data$targets))
}


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


dataset_to_score <- function(dataset){
  targets <- c(list(integer(0)), as.list(1:ncol(dataset$Y)))
  target.index <- dataset$targets == "control"
  for(i in seq_along(colnames(dataset$Y))){
    name_i <- colnames(dataset$Y)[[i]]
    target.index = target.index + (i+1)*(dataset$targets == name_i)
  }
  return(new("GaussL0penIntScore", data = dataset$Y, targets = targets,
             target.index = target.index, lambda = 0.5*log(nrow(dataset$Y)),
             intercept = FALSE, use.cpp = TRUE))
}
