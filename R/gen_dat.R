#' @title gen_dat
#' @description Generate data for a multisite trial with a given collection of
#'   features.
#' @param n.bar average site size, Default: 10
#' @param J number sites, Default: 30
#' @param p prop treated, Default: 0.5
#' @param tau.11.star Total amount of cross site treatment variation, both
#'   explained by covariate and not, Default: 0.3
#' @param rho2.0W Explanatory power of W for control outcomes, Default: 0.1
#' @param rho2.1W Explanatory power of W for average treatment impact, Default:
#'   0.5
#' @param ICC The ICC, Default: 0.7
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param zero.corr TRUE means treatment impact and mean site outcome are not
#'   correlated.  TRUE means they are negatively correlated to make the variance
#'   of the treatment group 1, Default: FALSE
#' @param ... Further parameters passed to gen_dat_model()
#' @return Dataframe of data!
#' @rdname gen_dat
#' @export
gen_dat <- function(n.bar = 10, J = 30, p = 0.5, tau.11.star = 0.3, rho2.0W = 0.1, rho2.1W = 0.5, ICC = 0.7, gamma.00 = 0, gamma.10 = 0.2, verbose = FALSE, zero.corr = FALSE, ... ) {
  sigma2.W <- 1
  gamma.01 <- sqrt(rho2.0W * ICC / sigma2.W)
  gamma.11 <- sqrt(rho2.1W * tau.11.star)
  tau.00 <- (1 - rho2.0W) * ICC
  if (zero.corr) {
    tau.01 <- 0
  } else {
    tau.01 <- -tau.11.star / 2 - gamma.01 * gamma.11 * sigma2.W
  }
  tau.11 <- (1 - rho2.1W) * tau.11.star
  sigma2.e <- 1 - ICC
  if (verbose) {
    scat( "tau.11* <- %.2f\tICC <- %.2f\trho2.Ws <- %.2f, %.2f\n", tau.11.star, ICC, rho2.0W, rho2.1W)
    scat( "tau.00* <- %.2f\n", gamma.01 ^ 2 * sigma2.W + tau.00)
    scat( "tau.11* <- %.2f\n", gamma.11 ^ 2 * sigma2.W + tau.11)
    scat( "sigma2.e* <- %.2f\n", sigma2.e)
  }
  gen_dat_model(n.bar = n.bar, J = J, p = p, gamma.00, gamma.01, gamma.10, gamma.11, tau.00, tau.01, tau.11, sigma2.e, verbose = verbose, ...)
}