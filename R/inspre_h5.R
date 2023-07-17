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
#' @param obs_to_use Sequence of bools. Indicator of columns of X to use in calculations.
#'   Useful for cross validation.
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
  return(as.data.frame(c(list(inst_id = inst_id, beta_obs = beta_inst_obs), as.list(c(beta_hat, se_hat)))))
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
calc_inst_effect_h5X <- function(inst_id, target, X, X_control, X_ids, X_vars){
  inst_obs <- X_ids == inst_id
  n_target <- sum(inst_obs)
  n_control <- ncol(X_control)
  n_total <- n_target + n_control
  this_feature <- which(X_vars == target)

  X_exp_target <- X[this_feature, inst_obs]
  X_exp_control <- X_control[this_feature, ]
  Z <- c(rep(1, n_target), rep(0, n_control))
  inst_effect <- cor(Z, c(X_exp_target, X_exp_control))
  inst_se <- sqrt((1 - inst_effect**2)/(n_total - 2))
  return(as.data.frame(list(inst_id = inst_id, target = target, inst_effect = inst_effect, inst_se = inst_se, n = n_total)))
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
predict_inspre_G_h5X <- function(res, X, beta, X_targets, X_vars, vars_to_use = NULL, obs_to_use = NULL){
  if(is.null(vars_to_use)){
    keep_vars <- rep(TRUE, length(X_vars))
  } else {
    keep_vars <- X_vars %in% vars_to_use
  }
  if(is.null(obs_to_use)){
    obs_to_use <- rep(TRUE, length(X_targets))
  }
  ZB <- outer(X_vars[keep_vars], X_targets[obs_to_use], "==")
  ZB <- ZB * beta

  data <- X[keep_vars, obs_to_use]
  X_hat <- array(0L, dim = c(dim(data), length(res$lambda)))
  eps_hat <- vector("numeric", length(res$lambda))
  for(i in 1:length(res$lambda)){
    X_hat[,,i] <- t(res$R_hat[,,i]) %*% data + ZB
    eps_hat[i] <- mean((data - X_hat[,,i])**2)
  }
  return(list(X_hat = X_hat, eps_hat = eps_hat))
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
  ImGi <- map(1:length(res$lambda), ~ solve(diag(D) - res$R_hat[,,.x]))
  ImGi <- array(do.call(cbind, ImGi), dim=c(dim(ImGi[[1]]), length(ImGi)))

  target_list <- intersect(setdiff(unique(X_targets[obs_to_use]), c('non-targeting')), names(beta))
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
  return(list(eps_hat = eps_hat))
}

