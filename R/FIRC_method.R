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
#' @param siteID ID of site (blocks are conBered nested in site).  If omitted, then blocks are conBered sites (the default).
#' @export
estimate.ATE.FIRC <- function( Yobs, Z, B, siteID = NULL, data=NULL, REML = FALSE, include.testing=TRUE, anova=FALSE, pool = FALSE ) {

    stopifnot( !( include.testing && REML ) )

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
            if ( !is.null( siteID ) ) {
                siteID = data[[siteID]]
            }
            data = data.frame( Yobs = eval( substitute( Yobs ), data ),
                               Z = eval( substitute( Z ), data ),
                               B = eval( substitute( B ), data) )
            data$siteID = siteID
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
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

    if ( pool ) {
        re.mod <- nlme::lme(Yobs ~ 0 + Z + B,
                            data = data,
                            random = ~ 0 + Z | siteID,
                            method = method,
                            control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    } else {
        re.mod <- nlme::lme(Yobs ~ 0 + Z + B,
                        data = data,
                        random = ~ 0 + Z | siteID,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    }

    if ( include.testing ) {
        # Test for cross site variation (???)
        if ( pool ) {
            re.mod.null <- nlme::gls(Yobs ~ 0 + Z + B,
                                     data=data,
                                     method = method,
                                     control=nlme::lmeControl(opt="optim",returnObject=TRUE))

        } else {
            re.mod.null <- nlme::gls(Yobs ~ 0 + Z + B,
                                 data=data,
                                 weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                                 method = method,
                                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))
        }
        #M0.null = lm( Yobs ~ 0 + B + Z, data=data )
        if ( anova ) {
            stopifnot( REML == FALSE )
            myanova = anova(re.mod.null, re.mod)
            p.value.anova = myanova[2,9]
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

    #extract tau
    vc <- nlme::VarCorr(re.mod)
    suppressWarnings(storage.mode(vc) <- "numeric")
    tau.hat <- vc["Z","StdDev"]

    #extract se of tau
    ## not source code: https://stackoverflow.com/questions/31694812/standard-error-of-variance-component-from-the-output-of-lmer/31704646#31704646
    ## soure code: https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-standard-errors-for-variance-components-from-mixed-models/

    #let logged sds first
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

    dat = catherine.gen.dat( 0.2, 1.0, 30, 50 )
    head( dat )
    describe.data( dat )

    #debug( estimate.ATE.FIRC )
    estimate.ATE.FIRC( Yobs, Z, B, dat )


    dat = catherine.gen.dat( 0.2, 0, 30, 50 )
    head( dat )
    describe.data( dat )

    #debug( estimate.ATE.FIRC )
    estimate.ATE.FIRC( Yobs, Z, B, dat )
}








if ( FALSE ) {

    localsource( "multisite_data_generators.R")

    dat = catherine.gen.dat( 0.2, 0.0, 30, 50 )
    describe.data( dat )

    estimate.ATE.FIRC( Yobs, Z, B, data=dat, pool=TRUE )

    dat = catherine.gen.dat( 0.2, 0.5, 30, 50 )
    describe.data( dat )

    estimate.ATE.FIRC( Yobs, Z, B, data=dat, pool=TRUE )

}




# Testing
if ( FALSE ) {

    df = gen.dat.no.cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

    # head( df )
    # describe.data( df )

    analysis.FIRC( Yobs, Z, B, df )

}



# Testing
# This testing is based on the DGP for small sample simulations
if ( FALSE ) {

  df = gen.dat.no.cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

  # head( df )
  # describe.data( df )

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
#'
analysis.FIRC.cov <- function( df, REML = FALSE, anova=FALSE ) {

  #fit multilevel model and extract pvalue
  re.mod <- nlme::lme(Yobs ~ 0 + Z * X + B,
                      data = df,
                      random = ~ 0 + Z | B,
                      weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                      method="ML",
                      control=nlme::lmeControl(opt="optim",returnObject=TRUE))

  # Test for cross site variation
  re.mod.null <- nlme::gls(Yobs ~ 0 + Z * X + B,
                           data=df,
                           weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                           method = "ML",
                           control=nlme::lmeControl(opt="optim", returnObject=TRUE))


  if ( anova ) {
    stopifnot( REML == FALSE )
    myanova = anova(re.mod.null, re.mod)
    p.value.anova = myanova[2,9]
    return( p.value.anova / 2 )  # divide by 2 by same logic as chi-squared test.
  } else {
    td = abs(as.numeric( deviance( re.mod ) - deviance( re.mod.null )))
    p.value = 0.5 * pchisq(td, 1, lower.tail = FALSE )
    return( p.value )
  }

}



