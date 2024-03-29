

#' @title Multilevel model based tests for cross site variation
#'
#' @description These minimal methods test for cross-site variation in
#'   a multi-site randomized trial using multilevel models. They
#'   assume a passed dataframe with specifically named columns.
#'
#'   In particular, Yobs is the outcome column, Z is the treatment
#'   column, W is a site-level covariate column, and sid is the site
#'   ID column (categorical).
#'
#' @name analysis_functions
#' @param df Dataframe to analyze.  Columns of \code{Z}, \code{sid},
#'   \code{Yobs}, and \code{W} assumed to exist, with those names.
#' @return P-value from the test.
NULL


#' @describeIn analysis_functions
#'
#' The idiosyncratic version tests for cross site variation using a likelihood
#' ratio test vs. a model with no random slope (but random intercept). This is
#' the RIRC version of such a test.
#'
#' @export
#' @importFrom lme4 lmer
#' @importFrom stats deviance pchisq
analysis_idio <- function(df) {
    M0 <- lme4::lmer(Yobs ~ 1 + Z + (1+Z|sid), data = df, REML = FALSE )
    M0.null <- lme4::lmer(Yobs ~ 1 + Z + (1|sid), data = df, REML = FALSE)

    # I _think_ this is what is suggested to handle the boundary by Snijders and Bosker
    td <- deviance( M0.null) - deviance(M0)
    0.5 * pchisq(td, 2, lower.tail = FALSE) + 0.5 * pchisq(td, 1, lower.tail = FALSE)
    # tst <- lrtest( M0, M0.null )
    # tst[[5]][[2]]
}



#' @describeIn analysis_functions
#'
#'   This is a combination test: it uses a liklihood ratio test on a
#'   model with both systematic and idiosyncratic variation (i.e., it
#'   has an interaction of the covariate and outcome, and also
#'   includes a random effect for treatment), comparing to a model
#'   which has neither.
#'
#' @export
analysis_combination <- function(df) {
  M0 <- lme4::lmer( Yobs ~ 1 + Z + W + Z:W + (Z|sid), data = df, REML = FALSE)
  M0.null <- lme4::lmer( Yobs ~ 1 + Z + W + (1|sid), data = df, REML = FALSE)
  # M0 <- lmer( Yobs ~ 1 + Z + W + Z:W + (Z|sid), data=df, REML <- TRUE )
  # M0.null <- lmer( Yobs ~ 1 + Z + W + (1|sid), data=df, REML <- TRUE )
  tst <- lmtest::lrtest(M0, M0.null)
  td <- deviance( M0.null) - deviance(M0)
  0.5 * pchisq(td, 3, lower.tail = FALSE) + 0.5 * pchisq(td, 2, lower.tail = FALSE)
  # tst[[5]][[2]]
}


#' @describeIn analysis_functions
#'
#'   Test for cross site variation by testing for a covariate
#'   predictive of variation. This version does not allow for random
#'   tx variation (random slope).
#'
#' @importFrom arm se.coef
#' @export
analysis_systematic <- function(df) {
    M0 <- lme4::lmer(Yobs ~ 1 + Z * W + (1|sid), data = df)
    tstat <- lme4::fixef(M0) / arm::se.coef(M0)$fixef
    2 * pnorm( - abs(tstat[["Z:W"]]))
}



#' @describeIn  analysis_functions
#'
#'   Systematic test with the random idiosyncratic variation.  This
#'   tests for systematic variation, ignoring any explicitly modeled
#'   random (idiosyncratic) treatment variation.
#'
#' @importFrom arm se.coef
#' @export

analysis_systematic_RTx <- function(df) {
    M0 <- lme4::lmer( Yobs ~ 1 + Z * W + (1+Z|sid), data = df)
    tstat <- lme4::fixef(M0) / arm::se.coef(M0)$fixef
    2 * pnorm( - abs(tstat[["Z:W"]]))
}




