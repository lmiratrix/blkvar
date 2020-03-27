#' Test for treatment variation across site.
#'
#' This is a combination test: it uses
#' a liklihood ratio test on model with both systematic and ideosyncratic
#' variation
#' @param df Dataframe to analyze.
#' @rdname analysis_idio
#' @export
analysis_combination <- function(df) {
  M0 <- lme4::lmer( Yobs ~ 1 + Z + X + Z:X + (Z|sid), data = df, REML = FALSE)
  M0.null <- lme4::lmer( Yobs ~ 1 + Z + X + (1|sid), data = df, REML = FALSE)
  # M0 <- lmer( Yobs ~ 1 + Z + X + Z:X + (Z|sid), data=df, REML <- TRUE )
  # M0.null <- lmer( Yobs ~ 1 + Z + X + (1|sid), data=df, REML <- TRUE )
  tst <- lmtest::lrtest(M0, M0.null)
  td <- deviance( M0.null) - deviance(M0)
  0.5 * pchisq(td, 3, lower.tail = FALSE) + 0.5 * pchisq(td, 2, lower.tail = FALSE)
  # tst[[5]][[2]]
}