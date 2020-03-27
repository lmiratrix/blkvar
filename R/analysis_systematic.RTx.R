#' Systematic test with the random ideosyncratic variation.  This tests for
#' systematic variation, ignoring explicitly modeled random treatment variation.
#'
#' @rdname analysis_idio
#' @param df Dataframe to analyze.
#' @importFrom arm se.coef
#' @export
analysis_systematic.RTx <- function(df) {
  M0 <- lme4::lmer( Yobs ~ 1 + Z * X + (Z|sid), data = df)
  tstat <- lme4::fixef(M0) / arm::se.coef(M0)$fixef
  2 * pnorm( - abs(tstat[["Z:X"]]))
}