#'
#' This file contains functions implementing the various methods for detecting
#' and estimating treatment effect variation across sites in a multi-site trial
#' using random-slope, random-intercept approaches.
#'

#' We have:
#'   * Ideosyncratic test
#'   * Testing for systematic covariate (two versions, one with random coeff and one not)
#'   * Hybrid test doing a likelihood ratio of the hybrid model vs. no-variation model
#'


#' Test for cross site variation using a liklihood ratio test vs. a model with no random slope (but random intercept).
#' This is the RIRC version of such a test.
#'
#' @importFrom lme4 lmer
#' @importFrom stats deviance pchisq
#' @export
analysis.idio = function( df ) {
    require( lme4 )
    M0 = lme4::lmer( Yobs ~ 1 + Z + (Z|sid), data=df, REML = FALSE )
    M0.null = lme4::lmer( Yobs ~ 1 + Z + (1|sid), data=df, REML= FALSE )

    # I _think_ this is what is suggested to handle the boundary by Snijders and Bosker
    td = deviance( M0.null ) - deviance( M0 )
    0.5 * pchisq(td, 2, lower.tail = FALSE ) + 0.5*pchisq(td, 1, lower.tail = FALSE )
    #tst = lrtest( M0, M0.null )
    #tst[[5]][[2]]
}


#' Test for cross site variation by testing for a covariate predictive of variation.
#'
#' Comment: Should this model allow for random tx variation or no?  (Currently it does not.)
#'
#' @export
analysis.systematic = function( df ) {
    require( lme4 )

    M0 = lme4::lmer( Yobs ~ 1 + Z * X + (1|sid), data=df )
    tstat = fixef( M0 ) / se.coef(M0)$fixef
    2 * pnorm( - abs( tstat[["Z:X"]] ) )
}

#' Systematic test with the random ideosyncratic component
#'
#' @export
analysis.systematic.RTx = function( df ) {
    require( lme4 )

    M0 = lme4::lmer( Yobs ~ 1 + Z * X + (Z|sid), data=df )
    tstat = fixef( M0 ) / se.coef(M0)$fixef
    2 * pnorm( - abs( tstat[["Z:X"]] ) )
}


#' Test for treatment variation across site.
#' This is a combination test: it uses a liklihood ratio test on model with both systematic and ideosyncratic variation
#'
#' @export
analysis.combination = function( df ) {
    require( lme4 )

    M0 = lme4::lmer( Yobs ~ 1 + Z + X + Z:X + (Z|sid), data=df, REML = FALSE )
    M0.null = lme4::lmer( Yobs ~ 1 + Z + X + (1|sid), data=df, REML = FALSE )
    #M0 = lmer( Yobs ~ 1 + Z + X + Z:X + (Z|sid), data=df, REML = TRUE )
    #M0.null = lmer( Yobs ~ 1 + Z + X + (1|sid), data=df, REML = TRUE )

    tst = lmtest::lrtest( M0, M0.null )
    tst
    td = deviance( M0.null ) - deviance( M0 )
    0.5 * pchisq(td, 3, lower.tail = FALSE ) + 0.5*pchisq(td, 2, lower.tail = FALSE )
    #tst[[5]][[2]]
}



#' Fit Random-intercept, random-slope model
#'
#' This is a wrapper for a simpler lmer call.
#'
#' @importFrom lme4 lmer VarCorr
#' @export
estimate.ATE.RIRC = function( Yobs, Z, sid, data=NULL, REML = FALSE ) {
    M0.full = lmer( Yobs ~ 1 + Z + (1+Z|sid), data=data, REML = REML )

    # get ATE and SE
    a = summary( M0.full )
    a = a$coefficients[2,]

    # Cross site variation
    tau.hat = sqrt( VarCorr( M0.full )$sid[2,2] )
    tau.hat

    return( list(ATE = a[[1]], SE.ATE = a[[2]],
                 tau.hat = tau.hat, SE.tau=NA))

}



if ( FALSE ) {
    library( blkvar )
    source( "R/multisite_data_generators.R")
    dat = catherine.gen.dat( 0.2, 0.2, 30, 50 )
    head( dat )
    analysis.idio.FIRC.pool( Yobs, Z, sid, data=dat )
    fit.FIRC.pool( Yobs, Z, sid, data=dat )

    estimate.ATE.RIRC( Yobs, Z, sid, data=dat )

    fit.FIRC( dat$Yobs, dat$Z, dat$sid )
}

