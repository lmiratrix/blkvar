#' @title Calculate power for cross-site variation
#' @description Calculate power for detecting cross site variation with a ratio
#'   test.  This uses work and results in Appendix A of Bloom and Spybrook
#'
#' @param J Number of sites
#' @param n.bar average people per site
#' @param tau Cross site VARIANCE of site-level ATEs.
#' @param p.tx proportion treated, Default: 0.5
#' @param ICC Control-side ICC, Default: 0.5
#' @param K Number of individual covariates used to explain control-side variation, Default: 0
#' @param R2.X Explanatory power of these individual level covariates, Default: 0
#' @param alpha Significance level, Default: 0.05
#' @return Power
#'
#' @examples
#' # Example in the text
#' calc_power( 80, 60, 0.02, p.tx=0.6, ICC=0.2, K = 0, R2.X=0.25 )
#'
#' @rdname calc_power
#' @export

calc_power <- function(J, n.bar, tau, p.tx = 0.5, ICC = 0.5, K = 0, R2.X = 0, alpha = 0.05) {
  omega <- calc_omega(n.bar, tau, p.tx, ICC, R2.X)
  df.num <- J - 1
  df.denom <- J * (n.bar - 2) - K
  F.crit <- qf(1 - alpha, df.num, df.denom)
  1 - pf( F.crit / (omega + 1), df.num, df.denom)
}