##
## FIRC model methods
##
## This file has implementation of the FIRC model (both with pooled variances
## across Tx and Co units and unpooled where Tx and Co get their own variance
## terms)
##
## The unpooled functions are taken from Catherine's implementation of Weiss et
## al methods.
##

library( dplyr )

## ---------- Both pooled and unpooled FIRC method ------------


#' Fit the FIRC model to estimate (1) ATE across sites and (2) cross site
#' treatment variation.
#'
#' Acknowledgement: Unpooled version taken and adapted from Catherine's weiss.tau() method.
#'
#' @param anova Use the anova() method to do the test for significance between
#'   the models.  FALSE means do the modified chi-squared test.
#' @param pool  TRUE means tx and co have same reBual variance. FALSE gives seperate estimates for each (recommended, default).
#' @param B Name of the block indicator.
#' @param siteID Character name of the ID variable of site (blocks are conBered nested in site).  If omitted, then blocks are considered sites (the default).
#' @export
estimate_ATE_FIRC <- function( Yobs, Z, B, siteID = NULL, control.formula = NULL,
                               data=NULL, REML = FALSE, include.testing=TRUE, anova=FALSE, pool = FALSE ) {

    stopifnot( !( include.testing && REML ) )

    if ( !is.null( control.formula ) ) {
        stopifnot( !is.null( data ) )
        stopifnot( !missing( "Yobs" ) )
    }

    # This code block takes the parameters of
    # Yobs, Z, B, siteID = NULL, data=NULL, ...
    # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
    if(!is.null(data)){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
            n.tx.lvls = length( unique( data$Z ) )
            stopifnot( n.tx.lvls == 2 )
            stopifnot( is.numeric( data$Yobs ) )
        } else {
            Yobs.n = as.character( substitute( Yobs ) )
            Z.n = as.character( substitute( Z ) )
            B.n = as.character( substitute( B ) )
            data = rename( data,
                            Yobs = Yobs.n,
                            Z = Z.n,
                            B = B.n )
            if ( !is.null( siteID ) ) {
                data$siteID = data[[siteID]]
                stopifnot( !is.null( data$siteID ) )
            }
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
        }
    }

    if ( is.null( siteID ) ) {
        data$siteID = data$B
    }


    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )


    #fit multilevel model and extract tau
    method = ifelse( REML, "REML", "ML" )

    # Make control variable function
    formula = make_FE_formula( "Yobs", "Z", "B", control.formula = control.formula, data = data)

    if ( pool ) {
        re.mod <- nlme::lme(formula,
                            data = data,
                            random = ~ 0 + Z | siteID,
                            method = method,
                            control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    } else {
        re.mod <- nlme::lme(formula,
                        data = data,
                        random = ~ 0 + Z | siteID,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    }

    if ( include.testing ) {
        # Test for cross site variation
        if ( pool ) {
            re.mod.null <- nlme::gls(formula,
                                     data=data,
                                     method = method,
                                     control=nlme::lmeControl(opt="optim",returnObject=TRUE))

        } else {
            re.mod.null <- nlme::gls(formula,
                                 data=data,
                                 weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                                 method = method,
                                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))
        }
        #M0.null = lm( Yobs ~ 0 + B + Z, data=data )
        if ( anova ) {
            stopifnot( REML == FALSE )
            # This tosses a warning that the linear model is being swapped into a MLM
            # "original model was of class "gls", updated model is of class "lme""
            myanova = suppressWarnings( lmtest::lrtest( re.mod.null, re.mod ) ) #anova(re.mod.null, re.mod)
            p.value.anova = myanova[2,5] #myanova[2,9]
            p.variation = ( p.value.anova / 2 )  # divide by 2 by same logic as chi-squared test.
            td = NA
        } else {
            td = as.numeric( deviance( re.mod.null ) - deviance( re.mod ) )
            p.variation = 0.5 * pchisq(td, 1, lower.tail = FALSE )
        }
    } else {
        p.variation = NA
        td = NA
    }

    # get ATE and SE
    ATE = nlme::fixef( re.mod )[[1]]
    SE.ATE = sqrt( vcov( re.mod )[1,1] )

    # extract tau (estimated cross site variation)
    vc <- nlme::VarCorr(re.mod)
    suppressWarnings(storage.mode(vc) <- "numeric")
    tau.hat <- vc["Z","StdDev"]

    # extract se of tau
    ## not source code: https://stackoverflow.com/questions/31694812/standard-error-of-variance-component-from-the-output-of-lmer/31704646#31704646
    ## soure code: https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-standard-errors-for-variance-components-from-mixed-models/

    # let logged sds first
    var <- re.mod$apVar

    if ( is.character( var ) ) {
        SE.tau = NA
    }
    else {
        par<-attr(var,"Pars")
        #transform to variances
        var.re <- exp(par)^2
        #use delta method from msm package
        SE.tau <- msm::deltamethod(~ exp(x1)^2,par,var)
    }

    return( list(ATE = ATE, SE.ATE = SE.ATE,
                 tau.hat = tau.hat, SE.tau=SE.tau,
                 p.variation = p.variation,
                 deviance = td ))
}







# For the debugging code to get the DGP files
localsource = function( filename ) {
    source( file.path( dirname( rstudioapi::getActiveDocumentContext()$path ), filename ) )
}


# Testing
if ( FALSE ) {

    localsource( "multisite_data_generators.R")

    dat = catherine_gen_dat( 0.2, 1.0, 30, 50 )
    head( dat )
    describe_data( dat )


    #debug( estimate_ATE_FIRC )
    estimate_ATE_FIRC( Yobs, Z, sid, data=dat )

    dat$X = dat$Y0 + rnorm( nrow(dat) )
    estimate_ATE_FIRC( Yobs, Z, sid, data=dat, control.formula = ~ X )

    dat = catherine_gen_dat( 0.2, 0, 30, 50 )
    head( dat )
    describe_data( dat )

    #debug( estimate_ATE_FIRC )
    estimate_ATE_FIRC( Yobs, Z, B, dat )
}








if ( FALSE ) {

    localsource( "multisite_data_generators.R")

    dat = catherine_gen_dat( 0.2, 0.0, 30, 50 )
    describe_data( dat )
    head( dat )

    estimate_ATE_FIRC( Yobs, Z, sid, data=dat, pool=TRUE )

    dat = catherine_gen_dat( 0.2, 0.5, 30, 50 )
    describe_data( dat )

    estimate_ATE_FIRC( Yobs, Z, B, data=dat, pool=TRUE )

}




# Testing
if ( FALSE ) {

    df = gen_dat_no_cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

    # head( df )
    # describe_data( df )

    analysis.FIRC( Yobs, Z, B, df )

}



# Testing
# This testing is based on the DGP for small sample simulations
if ( FALSE ) {

  df = gen_dat_no_cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

  # head( df )
  # describe_data( df )

  analysis.FIRC( df )

}


## ---------------- unpooled FIRC: accounting for covariates -----------------
##


#' Covariate adjusted test for cross-site variation
#'
#' This method fits a FIRC model but also includes a site-level covariate, X,
#' potentially predictive of treatment impact.
#'
#' This fits unpooled FIRC models.
#'
#' @param df  Dataframe to fit.  Needs Z, X, B, and Yobs columns.
#' @param X Site-level covariate ideally predictive of treatment variation.
#'
analysis_FIRC_cov <- function( Yobs, Z, B, X, siteID = NULL, data=NULL, REML = FALSE, anova=FALSE ) {


   # stopifnot( !( include.testing && REML ) )

    # This code block takes the parameters of
    # Yobs, Z, B, siteID = NULL, data=NULL, ...
    # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
    if(!is.null(data)){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]],
                               X = data[[4]] )
            n.tx.lvls = length( unique( data$Z ) )
            stopifnot( n.tx.lvls == 2 )
            stopifnot( is.numeric( data$Yobs ) )
        } else {
            if ( !is.null( siteID ) ) {
                siteID = data[[siteID]]
            }
            data = data.frame( Yobs = eval( substitute( Yobs ), data ),
                               Z = eval( substitute( Z ), data ),
                               B = eval( substitute( B ), data),
                               X = eval( substitute( X ), data) )
            data$siteID = siteID
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B,
                           X = X)
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
        }
    }

    if ( is.null( data$siteID ) ) {
        data$siteID = data$B
    }

    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )


    #fit multilevel model and extract tau
    method = ifelse( REML, "REML", "ML" )

  #fit multilevel model and extract pvalue
  re.mod <- nlme::lme(Yobs ~ 0 + Z + Z:X + B,
                      data = data,
                      random = ~ 0 + Z | B,
                      weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                      method="ML",
                      control=nlme::lmeControl(opt="optim",returnObject=TRUE))

  # Test for cross site variation
  re.mod.null <- nlme::gls(Yobs ~ 0 + Z + Z:X + B,
                           data = data,
                           weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                           method = "ML",
                           control=nlme::lmeControl(opt="optim", returnObject=TRUE))


  # Generate the different tests we can have for cross site variation in this model

  # Test 1: LR test (combination)
  if ( anova ) {
    stopifnot( REML == FALSE )
    myanova = anova(re.mod.null, re.mod)

    # TODO: Fix this!
    p.value.comb = myanova[2,9] / 2  # divide by 2 by same logic as chi-squared test.
  } else {
    td = abs(as.numeric( deviance( re.mod ) - deviance( re.mod.null )))
    p.value.comb = 0.5 * pchisq(td, 2, lower.tail = FALSE ) + 0.5 * pchisq(td, 1, lower.tail = FALSE )
  }

  # Test 2: Test the systematic coefficient
  SEs = sqrt( vcov( re.mod )["Z:X","Z:X"] )
  tstat = nlme::fixef( re.mod )[["Z:X"]] / SEs
  p.value.sys = 2 * pnorm( - abs( tstat ) )

  data.frame( method=c("FIRC combination", "FIRC systematic.RTx" ),
              pvalue = c( p.value.comb, p.value.sys ))

}



