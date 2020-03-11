#
# This file contains functions implementing the various methods for detecting
#  and estimating treatment effect variation across sites in a multi-site trial
# using random-slope, random-intercept approaches.
#

# We have:
#   * Ideosyncratic test
#   * Testing for systematic covariate (two versions, one with random coeff and one not)
#   * Hybrid test doing a likelihood ratio of the hybrid model vs. no-variation model
#




#' Test for cross site variation using a liklihood ratio test vs. a model with no random slope (but random intercept).
#' This is the RIRC version of such a test.
#'
#' @importFrom lme4 lmer
#' @importFrom stats deviance pchisq
#' @export
analysis_idio = function( df ) {
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
analysis_systematic = function( df ) {
    require( lme4 )

    M0 = lme4::lmer( Yobs ~ 1 + Z * X + (1|sid), data=df )
    tstat = lme4::fixef( M0 ) / se.coef(M0)$fixef
    2 * pnorm( - abs( tstat[["Z:X"]] ) )
}

#' Systematic test with the random ideosyncratic variation.  This tests for
#' systematic variation, ignoring explicitly modeled random treatment variation.
#'
#' @rdname analysis_idio
#' @export
analysis_systematic.RTx = function( df ) {
    require( lme4 )

    M0 = lme4::lmer( Yobs ~ 1 + Z * X + (Z|sid), data=df )
    tstat = lme4::fixef( M0 ) / se.coef(M0)$fixef
    2 * pnorm( - abs( tstat[["Z:X"]] ) )
}


#' Test for treatment variation across site.
#'
#' This is a combination test: it uses
#' a liklihood ratio test on model with both systematic and ideosyncratic
#' variation
#' @rdname analysis_idio
#' @export
analysis_combination = function( df ) {
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


#' Estimate the ATE using Random-Intercept, Random-Coefficient (RIRC) Models
#'
#' @inheritParams estimate_ATE_FIRC
#' @rdname estimate_ATE_RIRC
#'
#' @export
estimate_ATE_RIRC <- function( Yobs, Z, B, data=NULL, REML = FALSE, include.testing=TRUE, pool = FALSE,
                               control.formula = NULL ) {

    stopifnot( !( include.testing && REML ) )
    if ( !is.null( control.formula ) ) {
        stopifnot( !is.null( data ) )
        stopifnot( !missing( "Yobs" ) )
    }

    if( !is.null(data) ){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
        } else {
            d2 = data
            d2$Yobs = eval( substitute( Yobs ), data )
            d2$Z = eval( substitute( Z ), data )
            d2$B = eval( substitute( B ), data )
            data = d2
            rm( d2 )
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
    }
    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )

    #fit multilevel model and extract tau
    method = ifelse( REML, "REML", "ML" )

    formula = make_base_formula( control.formula=control.formula, data=data )
    if ( pool ) {
        re.mod <- nlme::lme(formula,
                            data = data,
                            random = ~ 1 + Z | B,
                            na.action=na.exclude,
                            method = method,
                            control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    } else {
        re.mod <- nlme::lme(formula,
                        data = data,
                        random = ~ 1 + Z | B,
                        weights = nlme::varIdent(form = ~ 1 | Z),
                        na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    }

    if ( include.testing ) {

        # Test for cross site variation (???)
        if ( pool ) {
            re.mod.null <- nlme::lme(formula,
                                     data = data,
                                     random = ~ 1 | B, na.action=na.exclude,
                                     method = method,
                                     control=nlme::lmeControl(opt="optim",returnObject=TRUE))

        } else {
            re.mod.null <- nlme::lme(formula,
                                 data = data,
                                 random = ~ 1 | B,
                                 weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                                 method = method,
                                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))
        }
        #M0.null = lm( Yobs ~ 0 + B + Z, data=data )
        td = as.numeric( deviance( re.mod.null ) - deviance( re.mod ) )
        p.variation = 0.5 * pchisq(td, 2, lower.tail = FALSE ) + 0.5 * pchisq(td, 1, lower.tail = FALSE )
    } else {
        p.variation = NA
        td = NA
    }

    # get ATE and SE
    ATE = nlme::fixef( re.mod )[[2]]
    SE.ATE = sqrt( vcov( re.mod )[2,2] )

    #extract tau
    vc <- nlme::VarCorr(re.mod)
    suppressWarnings(storage.mode(vc) <- "numeric")
    tau.hat <- vc["Z","StdDev"]


    return( list(ATE = ATE, SE.ATE = SE.ATE,
                 tau.hat = tau.hat, SE.tau=NA,
                 p.variation = p.variation,
                 deviance = td ))
}



#' Estimate the ATE using Random-Intercept, Constant-Coefficient (RICC) Model.
#'
#' This model has a single treatment coefficient, and a random intercept for the
#' site control average. So it is analogous to a fixed effect model, but with a
#' random effect.
#'
#' There is no test for cross site variation for this method, since we assume none.
#'
#' @inheritParams estimate_ATE_FIRC
#' @rdname estimate_ATE_RIRC
#'
#' @export
estimate_ATE_RICC <- function( Yobs, Z, B, data=NULL, REML = FALSE,
                               control.formula = NULL ) {
    if ( !is.null( control.formula ) ) {
        stopifnot( !is.null( data ) )
        stopifnot( !missing( "Yobs" ) )
    }

    if( !is.null(data) ){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
        } else {
            d2 = data
            d2$Yobs = eval( substitute( Yobs ), data )
            d2$Z = eval( substitute( Z ), data )
            d2$B = eval( substitute( B ), data )
            data = d2
            rm( d2 )
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
    }
    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )

    #fit multilevel model and extract tau
    method = ifelse( REML, "REML", "ML" )
    formula = make_base_formula( control.formula=control.formula, data=data )

    re.mod <- nlme::lme(formula,
                        data = data,
                        random = ~ 1 | B,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))


    # get ATE and SE
    ATE = nlme::fixef( re.mod )[[2]]
    SE.ATE = sqrt( vcov( re.mod )[2,2] )

    return( list(ATE = ATE, SE.ATE = SE.ATE,
                 tau.hat = NA, SE.tau=NA,
                 p.variation = NA,
                 deviance = NA ))
}




localsource = function( filename ) {
    source( file.path( dirname( rstudioapi::getActiveDocumentContext()$path ), filename ) )
}


if ( FALSE ) {
    localsource( "MLM_method.R" )
    localsource( "multisite_data_generators.R" )

    dat = catherine_gen_dat( 0.2, 5, 30, 50 )
    head( dat )
    describe_data( dat, Y0 = "Y0", Y1="Y1" )


    estimate_ATE_RIRC( Yobs, Z, B, data=dat )

    dat = catherine_gen_dat( 0.2, 0, 30, 50 )
    describe_data( dat, Y0 = "Y0", Y1="Y1" )

    estimate_ATE_RIRC( Yobs, Z, B, data=dat )
}





#' Fit Random-intercept, random-slope model
#'
#' This is a wrapper for a simpler lmer call.
#'
#' @importFrom lme4 lmer VarCorr
#' @export
estimate_ATE_RIRC_pool = function( Yobs, Z, B, data=NULL, REML = FALSE, include.testing=TRUE,
                                   control.formula = NULL ) {
    if ( !is.null( control.formula ) ) {
        stopifnot( !is.null( data ) )
        stopifnot( !missing( "Yobs" ) )
    }

    if( !is.null(data) ){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
        } else {
            if ( !is.null( siteID ) ) {
                siteID = data[[siteID]]
                stopifnot( !is.null( siteID ) )
            }
            d2 = data
            if ( !is.null( siteID ) ) {
                d2$siteID = data[[siteID]]
                stopifnot( !is.null( d2$siteID ) )
            }
            d2$Yobs = eval( substitute( Yobs ), data )
            d2$Z = eval( substitute( Z ), data )
            d2$B = eval( substitute( B ), data )
            data = d2
            rm( d2 )
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
        }
    }
    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )

    M0.full = lme4::lmer( Yobs ~ 1 + Z + (1+Z|B), data=data, REML = REML )

    if ( include.testing ) {
        M0.null = lme4::lmer( Yobs ~ 1 + Z + (1|B), data=data, REML = REML )

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
    tau.hat = sqrt( VarCorr( M0.full )$B[2,2] )
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
    dat = catherine_gen_dat( 0.2, 0.2, 30, 50 )
    head( dat )
    #analysis_idio.RIRC.pool( Yobs, Z, B, data=dat )
    #fit.RIRC.pool( Yobs, Z, B, data=dat )

    estimate_ATE_RIRC( Yobs, Z, sid, data=dat )
    estimate_ATE_RIRC_pool( Yobs, Z, sid, data=dat )

    #    fit.RIRC( dat$Yobs, dat$Z, dat$sid )
}

