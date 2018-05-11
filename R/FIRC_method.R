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


#' Fit the FIRC model to estimate ATE across sites as well as cross site treatment variation.
#'
#' Acknowledgement: taken and adapted from Catherine's weiss.tau() method.
#'
#' @export
estimate.ATE.FIRC <- function( Yobs, Z, sid, data=NULL, REML = FALSE ) {

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

    # Test for cross site variation (???)
    re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                             data=data,
                        weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                        method = method,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))
    #M0.null = lm( Yobs ~ 0 + sid + Z, data=data )
    td = as.numeric( deviance( re.mod ) - deviance( re.mod.null ) )
    p.variation = 0.5 * pchisq(td, 1, lower.tail = FALSE )


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
                 p.variation = p.variation ))
}




# Testing
if ( FALSE ) {
    source( "multisite_data_generators.R")

    dat = catherine.gen.dat( 0.2, 0.2, 30, 50 )
    head( dat )

    debug( estimate.ATE.FIRC )
    estimate.ATE.FIRC( Yobs, Z, sid, dat )
}



## ---------------- pooled methods -----------------

#' Fit the FIRC model, but using the same residual variance for both treatment
#' and control groups.  This is a simple call to `lmer()`.
#'
#' For testing for cross site variation, this will adjust the ratio stat p-value to account for the mixture of
#' chi-squared distributions.

#' @return List of things including a pvalue for cross-site variation
#'
#' @rdname estimate.ATE.FIRC
#' @export
estimate.ATE.FIRC.pool = function(  Yobs, Z, sid, data=NULL ) {

    M0 = lmer( Yobs ~ 0 + sid + Z + (0+Z|sid), data=data, REML = FALSE )
    #M0.null = lm( Yobs ~ 0 + sid + Z, data=df )
    #td = anova( M0, M0.null )$Chisq[[2]]
    #p.v = 0.5 * pchisq(td, 1, lower.tail = FALSE )
    M0.null = lm( Yobs ~ 0 + sid + Z, data=data )

    # obtain the estimate of cross-site variation
    tau.hat.FIRC.pool = sqrt( VarCorr( M0 )$sid[1,1] )

    # I _think_ this is what is suggested to handle the boundary by Snijders and Bosker
    td = deviance( M0.null ) - deviance( M0 )
    pvalue = 0.5 * pchisq(td, 1, lower.tail = FALSE ) + 0.5
    #tst = lrtest( M0, M0.null )
    #tst[[5]][[2]]

    list( ATE = tau.hat.FIRC.pool, p.variation = pvalue )
}




if ( FALSE ) {
    library( blkvar )
    source( "R/multisite_data_generators.R")
    dat = catherine.gen.dat( 0.2, 0.2, 30, 50 )
    head( dat )

    analysis.idio.FIRC.pool( Yobs, Z, sid, data=dat )

    fit.FIRC.pool( Yobs, Z, sid, data=dat )

    fit.FIRC( Yobs, Z, sid, data=dat )

    fit.FIRC( dat$Yobs, dat$Z, dat$sid )
}
