#' @title  Simplified version of gen_dat() with no W covariate.
#' @description Generate fake data for simulation studies
#' @inheritParams gen_dat
#' @param control.sd.Y1 Make correlation of random intercept and random slope
#'   such that the variance of the Y1s is 1.0, Default: TRUE
#' @param n.bar average people per site
#' @param J Number of sites
#' @param p prop treated, Default: 0.5
#' @param tau.11.star Total amount of cross site treatment variation
#' @param ICC The ICC, Default: 0.7
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param variable.n Allow n to vary around n.bar, Default: TRUE
#' @param control.sd.Y1 Make correlation of random intercept and random slope
#' @return Dataframe of individual data from a MLM DGP.
#' @rdname gen_dat_no_cov
#' @export
gen_dat_no_cov <- function(n.bar = 10, J = 30, p = 0.5, tau.11.star = 0.3, ICC = 0.7, gamma.00 = 0, gamma.10 = 0.2, verbose = FALSE, 
  variable.n = TRUE, control.sd.Y1 = TRUE, ... ) {
  tau.00 <- ICC
  tau.11 <- tau.11.star
  if (control.sd.Y1) {
    tau.01 <- -tau.11 / 2
  } else {
    tau.01 <- 0
  }
  sigma2.e <- 1 - ICC
  
  if (verbose) {
    scat( "tau.11* = %.2f\tICC = %.2f\n", tau.11.star, ICC)
    scat( "tau.00* = %.2f\n",  tau.00)
    scat( "tau.11* = %.2f\n",  tau.11)
    scat( "sigma2.e* = %.2f\n", sigma2.e)
  }
  gen_dat_model(n.bar = n.bar, J = J, p = p, gamma.00 = gamma.00, gamma.10 = gamma.10, gamma.01 = 0, gamma.11 = 0, tau.00 = tau.00, tau.01 = tau.01, tau.11 = tau.11,
    sigma2.e = sigma2.e, verbose = verbose, variable.n = variable.n, ...)
}