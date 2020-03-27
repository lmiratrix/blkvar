#' Test for cross site variation by testing for a covariate predictive of variation.
#'
#' Comment: Should this model allow for random tx variation or no?  (Currently it does not.)
#' @param df Dataframe to analyze.
#' @importFrom arm se.coef
#' @export
analysis_systematic <- function(df) {
  M0 <- lme4::lmer(Yobs ~ 1 + Z * X + (1|sid), data = df)
  tstat <- lme4::fixef(M0) / arm::se.coef(M0)$fixef
  2 * pnorm( - abs(tstat[["Z:X"]]))
}