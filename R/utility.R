calc_beta_res <- function(column, data){
  betareg_res <- betareg::betareg(as.formula(paste("eg_cent_mod ~ ", column)), data=data)
  return(betareg_res)
}


get_coeff <- function(betareg_res){
  summary_res <- summary(betareg_res)
  column <- rownames(summary_res$coefficients$mean)[2]
  return(as.data.frame(as.list(summary_res$coefficients$mean[2,]), row.names = column))
}


mae <- function(.x) mean(abs(.x))
sem <- function(.x) sd(abs(.x))/sqrt(length(.x))


# Based on https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/
plot_qq <- function(dat, dist=function(.x) .x, ci=0.95, kind='two-sided', decreasing=FALSE){
  n_dat <- length(dat)
  qqdata <- tibble(observed = sort(dat, decreasing), expected=dist(ppoints(n_dat)),
                   clower = dist(qbeta(p=(1-ci)/2, shape1=seq(n_dat), shape2=rev(seq(n_dat)))),
                   cupper = dist(qbeta(p=(1+ci)/2, shape1=seq(n_dat), shape2=rev(seq(n_dat)))))
  plot <- qqdata %>% ggplot(aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    geom_line() + labs(x = "Expected", y = "Observed")
  if(kind == 'two-sided'){
    plot <- plot + geom_segment(data = . %>% filter(expected == max(expected)),
                                aes(x = -expected, xend = expected, y = -expected, yend = expected),
                                linewidth = 1, alpha = 0.5, color = "grey30", lineend = "round")
  } else if(kind == 'one-sided'){
    plot <- plot + geom_segment(data = . %>% filter(expected == max(expected)),
                                aes(x = 0, xend = expected, y = 0, yend = expected),
                                linewidth = 1, alpha = 0.5, color = "grey30", lineend = "round")
  }
  return(plot)
}


two_side_log <- function(p) ifelse(p > 0.5, log10(2*(1-p)), -log10(2*p))
