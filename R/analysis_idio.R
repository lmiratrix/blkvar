#' Test for cross site variation using a liklihood ratio test vs. a model with no random slope (but random intercept).
#' This is the RIRC version of such a test.
#' @param df Dataframe to analyze.
#' @importFrom lme4 lmer
#' @importFrom stats deviance pchisq
#' @export

analysis_idio <- function(df) {
  M0 <- lme4::lmer(Yobs ~ 1 + Z + (Z|sid), data = df, REML = FALSE )
  M0.null <- lme4::lmer(Yobs ~ 1 + Z + (1|sid), data = df, REML = FALSE)
  # I _think_ this is what is suggested to handle the boundary by Snijders and Bosker
  td <- deviance( M0.null) - deviance(M0)
  0.5 * pchisq(td, 2, lower.tail = FALSE) + 0.5 * pchisq(td, 1, lower.tail = FALSE)
  # tst <- lrtest( M0, M0.null )
  # tst[[5]][[2]]
}