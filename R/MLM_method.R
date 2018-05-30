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

#' Systematic test with the random ideosyncratic variation.  This tests for
#' systematic variation, ignoring explicitly modeled random treatment variation.
#'
#' @rdname analysis.idio
#' @export
analysis.systematic.RTx = function( df ) {
    require( lme4 )

    M0 = lme4::lmer( Yobs ~ 1 + Z * X + (Z|sid), data=df )
    tstat = fixef( M0 ) / se.coef(M0)$fixef
    2 * pnorm( - abs( tstat[["Z:X"]] ) )
}


#' Test for treatment variation across site.
#'
#' This is a combination test: it uses
#' a liklihood ratio test on model with both systematic and ideosyncratic
#' variation
#' @rdname analysis.idio
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


estimate.ATE.RIRC <- function( Yobs, Z, sid, data=NULL, REML = FALSE, include.testing=TRUE) {

    # get our variables
    if ( is.null( data ) ) {
        data = data.frame( Yobs = Yobs, Z = Z, sid= factor(sid) )
        #Yobs = eval( substitute( Yobs ), data )
        #Z = eval( substitute( Z ), data )
        #sid = eval( substitute( sid ), data )
    } else {
        sid.name = as.character( quote( sid ) )
        data[ sid.name ] = factor( eval( substitute( sid ), data ) )
    }

    #fit multilevel model and extract tau
    method = ifelse( REML, "REML", "ML" )

    re.mod <- nlme::lme(Yobs ~ 1 + Z,
                        data = data,
                        random = ~ 1 + Z | sid,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))

    if ( include.testing ) {

        # Test for cross site variation (???)
        re.mod.null <- nlme::lme(Yobs ~ 1 + Z,
                                 data = data,
                                 random = ~ 1 | sid,
                                 weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                                 method = method,
                                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))

        #M0.null = lm( Yobs ~ 0 + sid + Z, data=data )
        td = as.numeric( deviance( re.mod.null ) - deviance( re.mod ) )
        p.variation = 0.5 * pchisq(td, 2, lower.tail = FALSE ) + 0.5 * pchisq(td, 1, lower.tail = FALSE )
    } else {
        p.variation = NA
        td = NA
    }

    # get ATE and SE
    ATE = nlme::fixef( re.mod )[[1]]
    SE.ATE = sqrt( vcov( re.mod )[1,1] )

    #extract tau
    vc <- nlme::VarCorr(re.mod)
    suppressWarnings(storage.mode(vc) <- "numeric")
    tau.hat <- vc["Z","StdDev"]


    return( list(ATE = ATE, SE.ATE = SE.ATE,
                 tau.hat = tau.hat, SE.tau=NA,
                 p.variation = p.variation,
                 deviance = td ))
}


localsource = function( filename ) {
    source( file.path( dirname( rstudioapi::getActiveDocumentContext()$path ), filename ) )
}


if ( FALSE ) {
    localsource( "MLM_method.R" )
    localsource( "multisite_data_generators.R" )

    dat = catherine.gen.dat( 0.2, 5, 30, 50 )
    head( dat )
    describe.data( dat, Y0 = "y0", Y1="y1" )


    estimate.ATE.RIRC( Yobs, Z, sid, data=dat )

    dat = catherine.gen.dat( 0.2, 0, 30, 50 )
    describe.data( dat, Y0 = "y0", Y1="y1" )

    estimate.ATE.RIRC( Yobs, Z, sid, data=dat )
}





#' Fit Random-intercept, random-slope model
#'
#' This is a wrapper for a simpler lmer call.
#'
#' @importFrom lme4 lmer VarCorr
#' @export
estimate.ATE.RIRC.pool = function( Yobs, Z, sid, data=NULL, REML = FALSE, include.testing=TRUE ) {
    M0.full = lme4::lmer( Yobs ~ 1 + Z + (1+Z|sid), data=data, REML = REML )


    if ( include.testing ) {
        M0.null = lme4::lmer( Yobs ~ 1 + Z + (1|sid), data=data, REML= FALSE )

        # I _think_ this is what is suggested to handle the boundary by Snijders and Bosker
        td = deviance( M0.null ) - deviance( M0.full )
        pv = 0.5 * pchisq(td, 2, lower.tail = FALSE ) + 0.5*pchisq(td, 1, lower.tail = FALSE )
    } else {
        pv = NA
        td = NA
    }

    # get ATE and SE
    a = summary( M0.full )
    a = a$coefficients[2,]

    # Cross site variation
    tau.hat = sqrt( VarCorr( M0.full )$sid[2,2] )
    tau.hat

    res = list(ATE = a[[1]], SE.ATE = a[[2]],
                 tau.hat = tau.hat, SE.tau=NA )
    if ( include.testing ) {
        res$p.variation = pv
        res$deviance = td
    }

    res

}


##
## Testing code
##

if ( FALSE ) {
    library( blkvar )
    source( "R/multisite_data_generators.R")
    dat = catherine.gen.dat( 0.2, 0.2, 30, 50 )
    head( dat )
    analysis.idio.RIRC.pool( Yobs, Z, sid, data=dat )
    fit.RIRC.pool( Yobs, Z, sid, data=dat )

    estimate.ATE.RIRC( Yobs, Z, sid, data=dat )

    fit.RIRC( dat$Yobs, dat$Z, dat$sid )
}

