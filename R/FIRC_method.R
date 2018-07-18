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

library( tidyverse )

## ---------- unpooled methods (from Catherine)  ------------


#' Fit the FIRC model to estimate (1) ATE across sites and (2) cross site
#' treatment variation.
#'
#' Acknowledgement: taken and adapted from Catherine's weiss.tau() method.
#'
#' @export
estimate.ATE.FIRC <- function( Yobs, Z, sid, data=NULL, REML = FALSE, include.testing=TRUE ) {

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

    re.mod <- nlme::lme(Yobs ~ 0 + Z + sid,
                        data = data,
                        random = ~ 0 + Z | sid,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))

    if ( include.testing ) {
        # Test for cross site variation (???)
        re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                                 data=data,
                                 weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                                 method = method,
                                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))

        #M0.null = lm( Yobs ~ 0 + sid + Z, data=data )
        td = as.numeric( deviance( re.mod.null ) - deviance( re.mod ) )
        p.variation = 0.5 * pchisq(td, 1, lower.tail = FALSE )
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
    estimate.ATE.FIRC( Yobs, Z, sid, dat )


    dat = catherine.gen.dat( 0.2, 0, 30, 50 )
    head( dat )
    describe.data( dat )

    #debug( estimate.ATE.FIRC )
    estimate.ATE.FIRC( Yobs, Z, sid, dat )
}



## ---------------- pooled methods -----------------

#' Fit the FIRC model, but using the same residual variance for both treatment
#' and control groups.  This can be done with a simple call to `lmer()`.
#'
#' For testing for cross site variation, this will adjust the ratio statistic's
#' p-value to account for the mixture of chi-squared distributions.

#' @return List of things including a pvalue for cross-site variation
#'
#' @rdname estimate.ATE.FIRC
#' @export
estimate.ATE.FIRC.pool = function(  Yobs, Z, sid, data=NULL, include.testing = TRUE ) {

    #    M0 = lmer( Yobs ~ 0 + Z + sid + (0+Z|sid), data=data, REML = FALSE )
    #    M0.null = lm( Yobs ~ 0 + Z + sid, data=data )

    # obtain the estimate of cross-site variation
    #    tau.hat.FIRC.pool = sqrt( VarCorr( M0 )$sid[1,1] )

    # Can we get deviance from the lm and the lmer models?
    #    td = deviance( M0.null ) - deviance( M0 )

    re.mod <- nlme::lme(Yobs ~ 0 + Z + sid,
                        data = data,
                        random = ~ 0 + Z | sid,
                        method = "ML",
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))

    if ( include.testing ) {
        # Test for cross site variation (???)
        re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                                 data=data,
                                 method = "ML",
                                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))


        td = as.numeric( deviance( re.mod.null ) - deviance( re.mod ) )
        p.variation = 0.5 * pchisq(td, 1, lower.tail = FALSE )
    } else {
        p.variation = NA
        td = NA
    }
    #pvalue = 0.5 * pchisq(-td, 1, lower.tail = FALSE )
    #tst = lrtest( M0, M0.null )
    #tst[[5]][[2]]

    vc <- nlme::VarCorr(re.mod)
    suppressWarnings(storage.mode(vc) <- "numeric")
    tau.hat.FIRC.pool <- vc["Z","StdDev"]

    list( ATE = fixef( re.mod )[[1]],
          tau.hat = tau.hat.FIRC.pool,
          p.variation = p.variation,
          deviance = td)
}






if ( FALSE ) {

    localsource( "multisite_data_generators.R")

    dat = catherine.gen.dat( 0.2, 0.0, 30, 50 )
    describe.data( dat )

    estimate.ATE.FIRC.pool( Yobs, Z, sid, data=dat )

    dat = catherine.gen.dat( 0.2, 0.5, 30, 50 )
    describe.data( dat )

    estimate.ATE.FIRC.pool( Yobs, Z, sid, data=dat )

}


## ---------------- unpooled FIRC: small sample sim -----------------
##
## This is an updated function from Masha's small sample simultions that
## extracts the correct p-value

#' Test for treatment heterogeniety with the FIRC model.
#'
#' @param anova Use the anova() method to do the test for significance between
#'   the models.  FALSE means do the modified chi-squared test.
analysis.FIRC <- function( data, REML = FALSE, anova=FALSE ) {

    # get variables
    #if ( is.null( df ) ) {
    #    data = data.frame( Yobs = Yobs, Z = Z, sid= factor(sid) )
    #} else {
    #    sid.name = as.character( quote( sid ) )
    #    data[ sid.name ] = factor( eval( substitute( sid ), df ) )
    #}

    #fit multilevel model and extract pvalue
    method = ifelse( REML, "REML", "ML" )

    re.mod <- nlme::lme(Yobs ~ 0 + Z + sid,
                        data = data,
                        random = ~ 0 + Z | sid,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))

    # Test for cross site variation
    re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                             data=data,
                             weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                             method = method,
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

# Testing
if ( FALSE ) {

    df = gen.dat.no.cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

    # head( df )
    # describe.data( df )

    analysis.FIRC( Yobs, Z, sid, df )

}



# Testing
# This testing is based on the DGP for small sample simulations
if ( FALSE ) {

  df = gen.dat.no.cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

  # head( df )
  # describe.data( df )

  analysis.FIRC( df )

}


## ---------------- unpooled FIRC: accoutning for covariates -----------------
##


analysis.FIRC.cov <- function( df ) {

  #fit multilevel model and extract pvalue
  re.mod <- nlme::lme(Yobs ~ 0 + Z * X + sid,
                      data = df,
                      random = ~ 0 + Z | sid,
                      weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                      method="ML",
                      control=nlme::lmeControl(opt="optim",returnObject=TRUE))

  # Test for cross site variation
  re.mod.null <- nlme::gls(Yobs ~ 0 + Z * X + sid,
                           data=df,
                           weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                           method = "ML",
                           control=nlme::lmeControl(opt="optim", returnObject=TRUE))


  myanova = anova(re.mod.null, re.mod)
  myanova[2,9]

}


