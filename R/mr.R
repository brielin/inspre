#' Uses MR to fit R matrix
#'
#' @param input_data A dataframe with columns `Exposure`, `Outcome`, `inst`,
#'  `beta_exp`, `beta_out`, `se_exp` and `se_out`. Each column is a list-column
#'  with values for the corresponding instrument for that pair of phenotypes.
#' @param mr_method MR Method to use. Currently supported: inverse variance
#'  weighting (IVW).
#' @param mr_function function taking 4 inputs `b_exp`, `b_out`, `se_exp`
#'  and `se_out`. Use to specify a custom MR function not implemented above.
#'
#' @export
fit_R_MR <- function(input_data, mr_method, mr_function=NULL){
  mr_method_func <- switch(
    mr_method,
    ps = function(beta_exp, beta_out, se_exp, se_out, weight){
      suppressWarnings(mr.raps::mr.raps(
        b_exp = beta_exp, b_out = beta_out, se_exp = se_exp, se_out = se_out))
      },
    aps = function(beta_exp, beta_out, se_exp, se_out, weight){
      suppressWarnings(mr.raps::mr.raps(
        b_exp = beta_exp, b_out = beta_out, se_exp = se_exp, se_out = se_out,
        over.dispersion = TRUE))
      },
    raps = function(beta_exp, beta_out, se_exp, se_out, weight){
      suppressWarnings(mr.raps::mr.raps(
        b_exp = beta_exp, b_out = beta_out, se_exp = se_exp, se_out = se_out,
        over.dispersion = TRUE, loss.function = "huber"))
      },
    egger_p = function(beta_exp, beta_out, se_exp, se_out, ...) {
      input <- MendelianRandomization::mr_input(bx = beta_exp, bxse = se_exp, by = beta_out, byse = se_out)
      egger_res <- MendelianRandomization::mr_egger(input, robust = FALSE, penalized = TRUE)
      return(list("beta.hat" = egger_res$Estimate, "beta.se" = egger_res$StdError.Est, "beta.p.value" = egger_res$Pvalue.Est))
      },
    egger = egger,
    mbe = function(beta_exp, beta_out, se_exp, se_out, ...){
      input <- MendelianRandomization::mr_input(bx = beta_exp, bxse = se_exp, by = beta_out, byse = se_out)
      mbe_res <- MendelianRandomization::mr_mbe(input, seed = NA, stderror = "simple")
      return(list("beta.hat" = mbe_res$Estimate, "beta.se" = mbe_res$StdError, "beta.p.value" = mbe_res$Pvalue))
      },
    median = function(beta_exp, beta_out, se_exp, se_out, ...){
      input <- MendelianRandomization::mr_input(bx = beta_exp, bxse = se_exp, by = beta_out, byse = se_out)
      median_res <- MendelianRandomization::mr_median(input)
      return(list("beta.hat" = median_res$Estimate, "beta.se" = median_res$StdError, "beta.p.value" = median_res$Pvalue))
      },
    ivw = function(beta_exp, beta_out, se_exp, se_out, ...){
      input <- MendelianRandomization::mr_input(bx = beta_exp, bxse = se_exp, by = beta_out, byse = se_out)
      ivw_res <- MendelianRandomization::mr_ivw(input)
      return(list("beta.hat" = ivw_res$Estimate, "beta.se" = ivw_res$StdError, "beta.p.value" = ivw_res$Pvalue))
      },
    mr_presso = function(beta_exp, beta_out, se_exp, se_out, ...){
      input <- data.frame("b_exp" = beta_exp, "b_out" = beta_out, "se_exp" = se_exp, "se_out" = se_out)
      mr_presso_res <- MRPRESSO::mr_presso(data = input, BetaOutcome = "b_out", BetaExposure = "b_exp", SdOutcome = "se_out", SdExposure = "se_exp",
                                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000)$`Main MR results`
      result_index = 2
      if(is.na(mr_presso_res$`Causal Estimate`[[result_index]])){
        result_index = 1
        }
      return(list("beta.hat" = mr_presso_res$`Causal Estimate`[[result_index]], "beta.se" = mr_presso_res$Sd[[result_index]], "beta.p.value" = mr_presso_res$`P-value`[[result_index]]))
      },
    mr_mix = function(beta_exp, beta_out, se_exp, se_out, ...){
      # TODO(brielin): This seems to flip the result?? Double check this.
      mrmix_res <- MRMix::MRMix(beta_exp, -beta_out, se_exp, se_out)
      return(list("beta.hat" = mrmix_res$theta, "beta.se" = mrmix_res$SE_theta, "beta.p.value" = mrmix_res$pvalue_theta))
      },
    custom = mr_function
  )

  mr_res <- input_data %>% dplyr::select(beta_exp, beta_out, se_exp, se_out) %>%
    purrr::pmap(mr_method_func)
  input_data %>% dplyr::select(-inst, -beta_exp, -beta_out, -se_exp, -se_out) %>%
    dplyr::mutate(R_hat = mr_res$beta.hat, SE_hat = mr_res$beta.se, p_val = mr_res$beta.p.value)
}
