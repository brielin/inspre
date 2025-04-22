## usethis namespace: start
#' @importFrom foreach %do%
## usethis namespace: end
NULL


#' Reads bigbrain formatted summary statistics.
#'
#' Reads the bigbrain formatter summary statistics file. Reads one line at
#' a time and saves each as Exposure, Outcome, Instrument, Effect size. Also
#' saves instrument information.
#'
#' @param filename String. Path to file you want to parse.
#' @param beta_col String. Name of the column containing beta.
#' @param p_col String. Name of the column containing p values.
#' @param id_col String. Name of the column containing variant ids (usually
#'  rsids).
#' @param se_col String. Name of the column containing standard errors of beta
#'  estimates.
#' @param p_thresh Float. Threshold to use for instrument inclusion.
#' @param chr_col String. Name of the column containing the SNP chromosome.
#' @param pos_col String. Name of the column containing the SNP position.
#' @param ref_col String. Name of the column containing the reference allele.
#' @param alt_col String. Name of the column containing the alternative allele.
#' @param feature_col String. Name of the column containing the feature.
#' @param delim String. File field delimiter.
#' @param filter_col String. Column to filter on. NULL for no additional filtering.
#' @param filter_min Float. Minimum value of `filter_col` for inclusion.
#' @param filter_max Float. Maximum value of `filter_col` for inclusion.
#'
#' @export
read_bigbrain <- function(filename, beta_col, p_col, id_col = "variant_id",
                          se_col = NULL, p_thresh = 5e-08, chr_col = "chr",
                          pos_col = "pos", ref_col = "ref", alt_col = "alt",
                          feature_col = "feature", delim = "\t",
                          filter_col = NULL, filter_max = NULL,
                          filter_min = NULL){
  last_snp <- NULL
  inst_df <- NULL
  effect_df <- NULL
  result_df <- NULL

  f <- function(x, pos){
    if(is.null(last_snp)){
      print("Initializing")
      # Initialize all matrices
      inst_df <<- tibble::tibble(inst=x[[id_col]], chr=x[[chr_col]], pos=x[[pos_col]],
                        ref=x[[ref_col]], alt=x[[alt_col]], n_eff=1)
      effect_df <<- tibble::tibble(inst=x[[id_col]], feature=x[[feature_col]],
                          beta=x[[beta_col]], se=x[[se_col]], p=x[[p_col]])
      last_snp <<- x[[id_col]]
      result_df <<- tibble::tibble(exp_out = character(), Exposure=character(), Outcome=character(),
                                   inst=list(),
                                   beta_exp=list(), se_exp=list(),
                                   beta_out=list(), se_out=list())

    }
    else if(x[[id_col]] != last_snp){
      cat(paste0("\nProcessing SNP ", x[[id_col]], ".\n"))
      effect_df <- effect_df %>% dplyr::distinct()
      inst_exposures <- dplyr::filter(effect_df, p < p_thresh)[[feature_col]]
      print(paste("It is an instrument for", length(inst_exposures), "exposures."))
      exp_data <- foreach::foreach(exp=inst_exposures, .combine=dplyr::bind_rows) %do% {
        data <- effect_df %>% dplyr::select(-p) %>%
          dplyr::filter(feature==exp) %>%
          dplyr::rename(Exposure=feature, beta_exp=beta, se_exp=se) %>%
          dplyr::cross_join(dplyr::select(effect_df, -inst, -p) %>%
                             dplyr::rename(Outcome=feature, beta_out=beta, se_out=se)) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(exp_out=paste(Exposure, Outcome, sep="_"), inst=list(inst),
                        beta_exp=list(beta_exp), se_exp=list(se_exp), beta_out=list(beta_out), se_out=list(se_out)) %>% dplyr::ungroup()
      }
      if(!is.null(exp_data)){
        print(paste("Joining with existing data."))
        start.time <- Sys.time()
        idx <- 1:nrow(exp_data)
        match_idx <- match(exp_data$exp_out, result_df$exp_out)
        idx[is.na(match_idx)] <- NA
        if(length(na.omit(match_idx)) > 0){
          print(length(na.omit(match_idx)))
          i=na.omit(idx)
          m=na.omit(match_idx)
          result_df[m, ] <<- result_df[m, ] %>% dplyr::mutate(
            inst = purrr::map2(result_df[m, "inst"][[1]], exp_data[i, "inst"][[1]], \(x,y) append(x, y)),
            beta_exp = purrr::map2(result_df[m, "beta_exp"][[1]], exp_data[i, "beta_exp"][[1]], \(x,y) append(x, y)),
            se_exp = purrr::map2(result_df[m, "se_exp"][[1]], exp_data[i, "se_exp"][[1]], \(x,y) append(x, y)),
            beta_out = purrr::map2(result_df[m, "beta_out"][[1]], exp_data[i, "beta_out"][[1]], \(x,y) append(x, y)),
            se_out = purrr::map2(result_df[m, "se_out"][[1]], exp_data[i, "se_out"][[1]], \(x,y) append(x, y)))
        }
        no_match <- exp_data[is.na(match_idx),]
        result_df <<- dplyr::bind_rows(result_df, no_match)
        print(Sys.time() - start.time)
      } else{
        print("Moving on.")
      }
      # Reset effect_df for new SNP
      effect_df <<- tibble::tibble(inst=x[[id_col]], feature=x[[feature_col]],
                          beta=x[[beta_col]], se=x[[se_col]], p=x[p_col])
      # Add new SNP info to inst_df
      inst_df <<- rbind(inst_df, tibble::tibble(inst=x[[id_col]], chr=x[[chr_col]],
                                        pos=x[[pos_col]], ref=x[[ref_col]],
                                        alt=x[[alt_col]], n_eff=1))
      last_snp <<- x[[id_col]]
    }
    else{
      # Continue adding SNP effects to effect_df.
      # Optionally skip SNPs based on filer_col.
      if(!is.null(filter_col)){
        if(is.null(filter_min) & is.null(filter_max)){
          stop("Filter col specified but no values given.")
        } else if(!is.null(filter_min)) {
          if(x[[filter_col]] < filter_min) next
        } else if(!is.null(filter_max)) {
          if(x[[filter_col]] > filter_max) next
        }

      }
      if(x[[se_col]] > 0){
        effect_df <<- rbind(effect_df, tibble::tibble(
          inst=x[[id_col]], feature=x[[feature_col]], beta=x[[beta_col]],
          se=x[[se_col]], p=x[[p_col]]))
      }
    }
  }

  callback_fn <- readr::SideEffectChunkCallback$new(f)
  readr::read_delim_chunked(filename, callback_fn, chunk_size=1, delim=delim)
  return(list(result_df=result_df, inst_df=inst_df))
}


#' Join bigbrain datasets that have been estimated on different chunks
#'
#' @param ... datasets to join
#'
#' @export
join_bigbrain <- function(...){
  input_res <- list(...)

  inst_df <<- tibble::tibble(inst=character(), chr=character(), pos=numeric(),
                             ref=character(), alt=character(), n_eff=numeric())
  result_df <<- tibble::tibble(exp_out = character(), Exposure=character(), Outcome=character(),
                               inst=list(),
                               beta_exp=list(), se_exp=list(),
                               beta_out=list(), se_out=list())
  aggregate_res <- function(res){
    inst_df <<- rbind(inst_df, res$inst_df)

    print(paste("Joining with existing data."))
    start.time <- Sys.time()
    idx <- 1:nrow(res$result_df)
    match_idx <- match(res$result_df$exp_out, result_df$exp_out)
    idx[is.na(match_idx)] <- NA
    if(length(na.omit(match_idx)) > 0){
      print(length(na.omit(match_idx)))
      i=na.omit(idx)
      m=na.omit(match_idx)
      result_df[m, ] <<- result_df[m, ] %>% dplyr::mutate(
        inst = purrr::map2(result_df[m, "inst"][[1]], res$result_df[i, "inst"][[1]], \(x,y) append(x, y)),
        beta_exp = purrr::map2(result_df[m, "beta_exp"][[1]], res$result_df[i, "beta_exp"][[1]], \(x,y) append(x, y)),
        se_exp = purrr::map2(result_df[m, "se_exp"][[1]], res$result_df[i, "se_exp"][[1]], \(x,y) append(x, y)),
        beta_out = purrr::map2(result_df[m, "beta_out"][[1]], res$result_df[i, "beta_out"][[1]], \(x,y) append(x, y)),
        se_out = purrr::map2(result_df[m, "se_out"][[1]], res$result_df[i, "se_out"][[1]], \(x,y) append(x, y)))
    }
    no_match <- res$result_df[is.na(match_idx),]
    result_df <<- dplyr::bind_rows(result_df, no_match)
    print(Sys.time() - start.time)
  }
  purrr::map(input_res, aggregate_res)
  return(list(inst_df=inst_df, result_df=result_df))
}
