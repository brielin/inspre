## usethis namespace: start
#' @useDynLib inspre, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
## usethis namespace: end
NULL


#' Simple parses for hdf5 metadata into a data.frame.
#'
#' Note that this was written with GPWS data in mind, thus
#' it can be used to parse the 'obs' and 'var' matrices
#' that correspond to the columns and rows of the data object
#' stored under 'X' in those files. YMMV for other hdf5 files.
#'
#' @param hfile hdf5r object.
#' @param entry String. Name of the entry in the hfile to parse.
#' @export
parse_hdf5_df <- function(hfile, entry = "obs"){
  colnames <- names(hfile[[entry]])

  cat_vals = NULL
  if('__categories' %in% colnames){
    cat_names = names(hfile[[entry]][['__categories']])
    get_cat <- function(cat_name){
      return(hfile[[entry]][['__categories']][[cat_name]][])
    }
    cat_vals = lapply(cat_names, get_cat)
    names(cat_vals) = cat_names
  }

  drop <- colnames %in% c('__categories')
  colnames <- colnames[!drop]
  get_col <- function(name){
    if(name %in% names(cat_vals)){
      indices = hfile[[entry]][[name]][]
      map = cat_vals[[name]]
      values = map[indices+1]  # Yay, R! :'[
    }
    else{
      values = hfile[[entry]][[name]][]
    }
    return(values)
  }

  df <- as.data.frame(lapply(colnames, get_col))
  names(df) <- colnames
  return(df)
}


#' Implements simple multiple-outcome instrumental variables regression.
#'
#' @param inst_id String. Id of instrument to calculate outcomes from.
#' @param target String. Target feature of the instrument.
#' @param X H5D. Example: `X=hfile[['X']]` where hfile is an hdf5r object
#'   containing a data matrix stored under label 'X'.
#' @param X_control Matrix. features by control observations. Usually
#'   pulled from X but stored in memory so you don't have to pull the
#'   same control observations off of disk for every instrument.
#' @param X_ids Sequence of strings with length equal to number of columns in
#'   X. Entries correspond to the instrument applied to that entry.
#' @param X_vars Sequence of strings with length equal to the number of rows
#'   in X. Name of the feature measured in that row.
#' @param vars_to_use Sequence of strings. Entries in X_vars to keep. Default
#'   NULL to keep all.
#' @param obs_to_use
#' @export
multiple_iv_reg_h5X <- function(inst_id, target, X, X_control, X_ids, X_vars, vars_to_use = NULL, obs_to_use = NULL){
  inst_obs <- X_ids == inst_id
  if(!is.null(obs_to_use)){
    inst_obs <- inst_obs & obs_to_use
  }
  n_target <- sum(inst_obs)
  n_control <- ncol(X_control)
  n_total <- n_target + n_control
  if(is.null(vars_to_use)){
    keep_vars <- rep(TRUE, length(X_vars))
  } else {
    keep_vars <- X_vars %in% vars_to_use
  }
  X_target <- X[keep_vars, inst_obs]
  this_feature <- which(X_vars[keep_vars] == target)
  beta_inst_obs <- sum(X_target[this_feature, ])/n_target
  X_exp_target <- X_target[this_feature, ]
  X_exp_control <- X_control[which(X_vars == target), ]  # X_control has all features.
  sums_target <- rowSums(X_target)
  beta_hat <- sums_target/sums_target[this_feature]
  resid_sum <- rowSums((X_target - outer(beta_hat, X_exp_target))**2) +
    rowSums((X_control[keep_vars,] - outer(beta_hat, X_exp_control))**2)
  var_resid <- resid_sum/(n_total*n_target*beta_inst_obs**2)
  se_hat <- sqrt(var_resid)
  names(beta_hat) <- paste(X_vars[keep_vars], "beta_hat", sep="_")
  names(se_hat) <- paste(X_vars[keep_vars], "se_hat", sep="_")
  return(as.data.frame(c(list(inst_id = inst_id, target = target, beta_obs = beta_inst_obs), as.list(c(beta_hat, se_hat)))))
}


#' Calculates correlation and standard error of instrument with target feature.
#'
#' @param inst_id String. Id of instrumental to calculate outcomes from.
#' @param target String. Target feature of the instrument.
#' @param X H5D. Example: `X=hfile[['X']]` where hfile is an hdf5r object
#'   containing a data matrix stored under label 'X'.
#' @param X_control Matrix. features by control observations. Usually
#'   pulled from X but stored in memory so you don't have to pull the
#'   same control observations off of disk for every instrument.
#' @param X_ids Sequence of strings with length equal to number of columns in
#'   X. Entries correspond to the instrument applied to that entry.
#' @param X_vars Sequence of strings with length equal to the number of rows
#'   in X. Name of the feature measured in that row.
#' @export
calc_inst_effect_h5X <- function(inst_id, target, X, X_control, X_ids, X_vars){
  inst_obs <- X_ids == inst_id
  n_target <- sum(inst_obs)
  n_control <- ncol(X_control)
  n_total <- n_target + n_control
  this_feature <- which(X_vars == target)

  X_exp_target <- X[this_feature, inst_obs]
  X_exp_control <- X_control[this_feature, ]
  Z <- c(rep(1, n_target), rep(0, n_control))
  inst_cor <- cor(Z, c(X_exp_target, X_exp_control))
  inst_beta <- cov(Z, c(X_exp_target, X_exp_control))/var(Z)
  cor_se <- sqrt((1 - inst_cor**2)/(n_total - 2))
  return(as.data.frame(list(inst_id = inst_id, target = target, inst_beta = inst_beta, inst_cor = inst_cor, cor_se = cor_se, n = n_total)))
}


#' Calculates model predictions given new data, with cross-validation errors.
#'
#' This function uses the mean((X - XG - ZB)**2) error formulation.
#' TODO: At the moment this pulls all the data into memory which is probably
#'   fine since its normally only 10-20% of the total data size. Consider
#'   modifying to pull one target at a time to match fit functionality.
#'
#' @param res Inspre result. Result of running fit_inspre for a single or
#'   sequence of hyperparameter values.
#' @param X H5D. Example: `X=hfile[['X']]` where hfile is an hdf5r object
#'   containing a data matrix stored under label 'X'.
#' @param beta Sequence of floats. Obsevered scale estimates of instrument
#'   effects.
#' @param X_targets Sequence of strings with length equal to number of columns in
#'   X. Entries correspond to the target of the instrument applied to that entry.
#' @param X_vars Sequence of strings with length equal to the number of rows
#'   in X. Name of the feature measured in that row.
#' @param vars_to_use Sequence of strings. Entries in X_vars to keep. Default
#'   NULL to keep all.
#' @param obs_to_use Sequence of bools. Indicator of columns of X to use in calculations.
#'   Useful for cross validation. NULL to keep all.
predict_inspre_G_h5X <- function(res, X, X_ntc, beta, X_targets, X_vars,
                                 vars_to_use = NULL, obs_to_use = NULL){
  if(is.null(vars_to_use)){
    keep_vars <- rep(TRUE, length(X_vars))
  } else {
    keep_vars <- X_vars %in% vars_to_use
  }
  if(is.null(obs_to_use)){
    obs_to_use <- rep(TRUE, length(X_targets))
  }
  # Drops control samples and anything not in our target set.
  obs_to_use <- obs_to_use & (X_targets %in% X_vars[keep_vars])
  ZB <- outer(X_vars[keep_vars], X_targets[obs_to_use], "==")
  ZB <- ZB * beta
  ZB <- cbind(ZB, matrix(0, nrow=nrow(ZB), ncol=ncol(X_ntc)))

  data <- cbind(X[keep_vars, obs_to_use], X_ntc[keep_vars,])
  X_hat <- array(0L, dim = c(dim(data), length(res$lambda)))
  eps_hat <- vector("numeric", length(res$lambda))
  for(i in 1:length(res$lambda)){
    X_hat[,,i] <- t(res$R_hat[,,i]) %*% data + ZB
    eps_hat[i] <- mean((data - X_hat[,,i])**2)
  }
  return(eps_hat)
}


#' Calculates model predictions given new data, with cross-validation errors.
#'
#' This function uses the mean((ZB(I-G)^-1)**2) error formulation.
#' TODO: At the moment this pulls all the data into memory which is probably
#'   fine since its normally only 10-20% of the total data size. Consider
#'   modifying to pull one target at a time to match fit functionality.
#'
#' @param res Inspre result. Result of running fit_inspre for a single or
#'   sequence of hyperparameter values.
#' @param X H5D. Example: `X=hfile[['X']]` where hfile is an hdf5r object
#'   containing a data matrix stored under label 'X'.
#' @param beta Sequence of floats. Obsevered scale estimates of instrument
#'   effects. Entry name corresponds to targeted gene.
#' @param X_targets Sequence of strings with length equal to number of columns in
#'   X. Entries correspond to the instrument applied to that entry.
#' @param X_vars Sequence of strings with length equal to the number of rows
#'   in X. Name of the feature measured in that row.
#' @param vars_to_use Sequence of strings. Entries in X_vars to keep. Default
#'   NULL to keep all.
#' @param obs_to_use Sequence of bools. Indicator of columns of X to use in calculations.
#'   Useful for cross validation. NULL to keep all.
predict_inspre_inv_h5X <- function(res, X, beta, X_targets, X_vars, vars_to_use = NULL, obs_to_use = NULL){
  if(is.null(vars_to_use)){
    keep_vars <- rep(TRUE, length(X_vars))
  } else {
    keep_vars <- X_vars %in% vars_to_use
  }
  if(is.null(obs_to_use)){
    obs_to_use <- rep(TRUE, length(X_targets))
  }
  D <- nrow(res$R_hat)
  # This is probably not working for CV iters
  # Use approx instead?
  # ImGi <- map(1:length(res$lambda), ~ solve(diag(D) - res$R_hat[,,.x]))
  ImGi <- map(1:length(res$lambda), ~ res$U[,,.x] * diag(res$V[,,.x])) # = D[V] %*% U
  ImGi <- array(do.call(cbind, ImGi), dim=c(dim(ImGi[[1]]), length(ImGi)))

  target_list <- unlist(intersect(setdiff(unique(X_targets[obs_to_use]), c('non-targeting')), names(beta)))
  calc_error <- function(target, ImGi_target, beta_target){
    Xd <- X[keep_vars, (X_targets == target) & obs_to_use]
    Xd_hat <- ImGi_target * beta_target
    eps <- Xd - Xd_hat
    return(sum(eps**2))
  }

  eps_hat <- vector("numeric", length(res$lambda))
  for(i in 1:length(res$lambda)){
    eps_hat_i <- sum(unlist(pmap(list(target_list, array_branch(ImGi[,,i], 1), beta[target_list]), calc_error)))
    eps_hat[i] <- eps_hat_i/sum(obs_to_use)
  }
  return(eps_hat)
}

#' Fits inverse sparse regression model.
#'
#' See also inspre::inspre() for more details.
#'
#' @param X H5D. Example: `X=hfile[['X']]` where hfile is an hdf5r object
#'   containing a data matrix stored under label 'X'.
#' @param X_control Matrix. features by control observations. Usually
#'   pulled from X but stored in memory so you don't have to pull the
#'   same control observations off of disk for every instrument.
#' @param X_ids Sequence of strings with length equal to number of columns in
#'   X. Entries correspond to the instrument applied to that entry.
#' @param X_vars Sequence of strings with length equal to the number of rows
#'   in X. Name of the feature measured in that row.
#' @param targets List. Entries correspond to variables in X (`X_vars`), names
#'   correspond to the instrument targeting that variable (`X_ids`). These targets
#'   will be used to calculate causal effect sizes.
#' @param vars_to_use Sequence of strings. Entries in X_vars to keep. Default
#'   NULL to keep all.
#' @param obs_to_use Sequence of bools. Indicator of columns of X to use in calculations.
#'   Useful for cross validation.
#' @param max_med_ratio Float or NULL. Passed through to `make_weights`, NULL
#'   for no weights.
#' @param filter Bool. True to filter the produced TCE matrix with `fitler_tce`.
#' @param rho Float. Initial learning rate for ADMM.
#' @param lambda Float, sequence of floats of NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmicly spaced set of values between
#'   the maximum absolute off diagonal element of X and lambda_min_ratio
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
#' @param DAG Bool. True to resitrict solutions to approximate DAGs.
#' @export
fit_inspre_from_h5X <- function(X, X_control, X_ids, X_vars, targets,
                                max_med_ratio = NULL, filter = TRUE,
                                rho = 100.0, lambda = NULL,
                                lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                                gamma = NULL, its = 100, delta_target = 1e-4,
                                verbose = 2, cv_folds = 0, mu = 10, tau = 1.5, solve_its = 10,
                                ncores = 1, min_nz = 0.01, warm_start = FALSE,
                                constraint = "UV", DAG = FALSE){
  X_targets = map(X_ids, ~ targets[[.x]])
  if(verbose){cat("Fitting model with full dataset\n")}
  inst_effects <- map2(names(targets), targets, ~ multiple_iv_reg_h5X(
    .x, .y, X, X_ntc, X_ids, X_vars, targets)) %>% list_rbind()
  beta_obs <- inst_effects$beta_obs
  names(beta_obs) <- inst_effects$inst_id
  R_hat <- data.matrix(select(inst_effects, ends_with('beta_hat')) %>%
                         rename_with(~sub('_beta_hat', '', .x)))
  SE_hat <- data.matrix(select(inst_effects, ends_with('se_hat')) %>%
                          rename_with(~sub('_se_hat', '', .x)))
  rownames(R_hat) <- names(beta_obs)
  rownames(SE_hat) <- names(beta_obs)

  if(filter){
    filtered <- filter_tce(R_hat, SE_hat)
    R_hat <- filtered$R
    SE_hat <- filtered$SE
    keep_rows <- rownames(R_hat)
    keep_cols <- colnames(R_hat)
  }
  D <- ncol(R_hat)

  weights <- NULL
  if(!is.null(max_med_ratio)){
    weights <- inspre::make_weights(SE_hat, max_med_ratio = max_med_ratio)
  }
  full_res <- fit_inspre_from_R(R_hat, W = weights, rho = rho, lambda = lambda,
                                lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
                                gamma = gamma, its = its, delta_target = delta_target,
                                verbose = verbose, train_prop = 1,
                                cv_folds = 0, mu = mu, tau = tau, solve_its = solve_its,
                                ncores = ncores, warm_start = warm_start, min_nz = min_nz,
                                constraint = constraint, DAG = DAG)

  if(cv_folds > 0){
    # For simplicity later, we make two sets here and use the case indices from
    # the first (all) and control indicies from the second (control).
    folds_all <- caret::createFolds(X_ids, k=cv_folds)
    folds_control <- caret::createFolds(1:ncol(X_control), k=cv_folds)
    eps_hat_G <- array(0.0, length(full_res$lambda))
    eps_hat_I <- array(0.0, length(full_res$lambda))
    eps_hat_Rtr <- array(0.0, length(full_res$lambda)) # W from tr
    eps_hat_Rte <- array(0.0, length(full_res$lambda))# W from te
    xi_mat <- array(0, dim = c(D, D, length(full_res$lambda)))
    for(i in 1:cv_folds){
      if (verbose){
        cat(sprintf("Cross-validation iteration %d.\n", i))
      }
      fold_all <- folds_all[[i]]
      fold_control <- folds_control[[i]]
      train_obs_all <- !(1:length(X_ids) %in% fold_all)
      train_obs_control <-!(1:ncol(X_control) %in% fold_control)
      inst_effects <- map2(names(targets), targets, ~ multiple_iv_reg_h5X(
        .x, .y, X, X_ntc[, train_obs_control], X_ids, X_vars, targets, train_obs_all)) %>% list_rbind()
      beta_obs <- inst_effects$beta_obs
      names(beta_obs) <- inst_effects$inst_id
      R_hat_cv <- data.matrix(dplyr::select(inst_effects, ends_with('beta_hat')) %>%
                             dplyr::rename_with(~sub('_beta_hat', '', .x)))
      SE_hat_cv <- data.matrix(dplyr::select(inst_effects, ends_with('se_hat')) %>%
                              dplyr::rename_with(~sub('_se_hat', '', .x)))
      rownames(R_hat_cv) <- names(beta_obs)
      rownames(SE_hat_cv) <- names(beta_obs)

      if(filter){
        # Still want to set large values to NA but need to use the same columns/rows.
        R_hat_cv <- R_hat_cv[keep_rows, keep_cols]
        SE_hat_cv <- SE_hat_cv[keep_rows, keep_cols]
        filtered <- filter_tce(R_hat_cv, SE_hat_cv, max_nan_perc = 1)
        R_hat_cv <- filtered$R
        SE_hat_cv <- filtered$SE
      }
      weights <- NULL
      if(!is.null(max_med_ratio)){
        weights <- inspre::make_weights(SE_hat_cv, max_med_ratio = max_med_ratio)
      }
      cv_res <- fit_inspre_from_R(R_hat_cv, W = weights, rho = rho, lambda = lambda,
                                  lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
                                  gamma = gamma, its = its, delta_target = delta_target,
                                  verbose = verbose, train_prop = 1,
                                  cv_folds = 0, mu = mu, tau = tau, solve_its = solve_its,
                                  ncores = ncores, warm_start = warm_start, min_nz = min_nz,
                                  constraint = constraint, DAG = DAG)
      V_nz <- abs(cv_res$V) > min_nz
      xi_mat <- xi_mat + V_nz

      test_obs_all <- 1:length(X_ids) %in% fold_all
      test_obs_control <- 1:ncol(X_control) %in% fold_control

      eps_G <- predict_inspre_G_h5X(cv_res, X, X_ntc[,test_obs_control], beta_obs,
                                    X_targets, X_vars, targets, test_obs_all)
      beta_test <- beta_obs
      names(beta_test) <- inst_effects$target
      eps_I <- predict_inspre_inv_h5X(cv_res, X, beta_test, X_targets, X_vars, targets, test_obs_all)


      # Testing new approach.
      test_obs_all <- (1:length(X_ids) %in% fold_all)
      test_obs_control <-(1:ncol(X_control) %in% fold_control)
      inst_effects_te <- map2(names(targets), targets, ~ multiple_iv_reg_h5X(
        .x, .y, X, X_ntc[, test_obs_control], X_ids, X_vars, targets, test_obs_all)) %>% list_rbind()
      beta_obs_te <- inst_effects_te$beta_obs
      names(beta_obs_te) <- inst_effects_te$inst_id
      R_hat_te <- data.matrix(dplyr::select(inst_effects_te, ends_with('beta_hat')) %>%
                                dplyr::rename_with(~sub('_beta_hat', '', .x)))
      SE_hat_te <- data.matrix(dplyr::select(inst_effects_te, ends_with('se_hat')) %>%
                                 dplyr::rename_with(~sub('_se_hat', '', .x)))
      rownames(R_hat_te) <- names(beta_obs)
      rownames(SE_hat_te) <- names(beta_obs)
      if(filter){
        # Still want to set large values to NA but need to use the same columns/rows.
        R_hat_te <- R_hat_te[keep_rows, keep_cols]
        SE_hat_te <- SE_hat_te[keep_rows, keep_cols]
        filtered <- filter_tce(R_hat_te, SE_hat_te, max_nan_perc = 1)
        R_hat_te <- filtered$R
        SE_hat_te <- filtered$SE
      }
      weights_te <- NULL
      if(!is.null(max_med_ratio)){
        weights_te <- inspre::make_weights(SE_hat_te, max_med_ratio = max_med_ratio)
      }
      eps_Rte <- unlist(purrr::map(1:length(full_res$lambda), ~ mean((weights_te * (cv_res$U[,,.x] - R_hat_te))**2, na.rm=T)))
      eps_Rtr <- unlist(purrr::map(1:length(full_res$lambda), ~ mean((weights * (cv_res$U[,,.x] - R_hat_te))**2, na.rm=T)))
      ##

      eps_hat_G <- eps_hat_G + eps_G
      eps_hat_I <- eps_hat_I + eps_I
      eps_hat_Rte <- eps_hat_Rte + eps_Rte
      eps_hat_Rtr <- eps_hat_Rtr + eps_Rtr
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
    eps_hat_Rte <- eps_hat_Rte/cv_folds
    eps_hat_Rtr <- eps_hat_Rtr/cv_folds
    full_res$eps_hat_G <- eps_hat_G
    full_res$eps_hat_I <- eps_hat_I
    full_res$eps_hat_Rte <- eps_hat_Rte
    full_res$eps_hat_Rtr <- eps_hat_Rtr
  }
  # Hack until I have time to change all these variable names.
  # TODO: add these to other methods.
  full_res$G_hat <- full_res$R_hat
  full_res$R_hat <- R_hat
  full_res$SE_hat <- SE_hat
  return(full_res)
}
