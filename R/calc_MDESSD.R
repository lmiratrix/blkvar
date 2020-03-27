#' @title Calculate MDESSD for cross-site variation
#' @description Calculate MDESSD (minimum detectable effect size for cross site standard deviation of effects
#' assuming a ratio test.  This uses work and results found in Appendix A of Bloom and Spybrook
#'
#' @inheritParams calc_power
#' @return MDESSD for given scenario.
#' @seealso calc_power
calc_MDESSD <- function(J, n.bar, p.tx = 0.5, ICC = 0.5, K = 0, R2.X = 0, alpha=0.05) {
  df.num <- J - 1
  df.denom <- J * (n.bar - 2) - K
  F.crit <- qf(1 - alpha, df.num, df.denom)
  F.crit.20 <- qf(0.20, df.num, df.denom)
  rat <- (1 - ICC) * (1 - R2.X) / (n.bar * p.tx * (1 - p.tx))
  sqrt(rat * (F.crit / F.crit.20 - 1))
}