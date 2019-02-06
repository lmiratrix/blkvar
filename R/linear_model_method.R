#
# Using linear regression to deal with multi-site or blocked experiments
#




#' Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


grab.SE = function( MOD, coef="Z" ) {
    vc = vcov( MOD )
    sqrt( vc[ coef, coef ] )
}


fixed.effect.estimators = function( dat ) {

    # simple linear model
    M0 = lm( Yobs ~ 0 + Z + blk, data=dat )
    SE.lm = summary( M0 )$coeff["Z",2]

    # Huber-White SEs
    vcov_sand = sandwich::vcovHC(M0, type = "HC1")
    SE.lm.sand <-sqrt( vcov_sand[1,1] )

    # Cluster robust SEs
    vcov_clust = sandwich::vcovCL( M0, dat$blk )
    SE.lm.clust = sqrt( vcov_clust[1,1] )

    FEmodels = data.frame( method=c("fixed effects", "fixed effects (sand SE)", "fixed effects (cluster SE)" ),
                           tau = c( coef(M0)[["Z"]], coef(M0)[["Z"]], coef(M0)[["Z"]] ),
                           SE = c( SE.lm, SE.lm.sand, SE.lm.clust ),
                           stringsAsFactors = FALSE )

    FEmodels
}


#' Individual Weighted regression models using the lm and precision weighting
#'
#' WARNING: This is the wrong kind of regression.  The weights are correct, but
#' the use is wrong.
#'
weighted.linear.estimators.naive = function( dat ) {
    Z.bar = mean( dat$Z )
    n = nrow( dat )
    J = length( unique( dat$blk ) )
    n.bar = n / J
    dat = dat %>% dplyr::group_by( blk ) %>%
        dplyr::mutate( p = mean( Z ),
                       nj = n(),
                       w.orig = ifelse( Z, 1/p, 1/(1-p) ),
                       weight = ifelse( Z, Z.bar / p, (1-Z.bar)/(1-p) ),
                       weight.site = weight * n.bar / nj ) %>%
        dplyr::ungroup()

    M0w = lm( Yobs ~ Z + blk, weights=dat$w.orig, data=dat )
    SE.w = summary( M0w )$coeff["Z",2]

    M0w2 = lm( Yobs ~ Z + blk, weights=dat$weight, data=dat )
    SE.w2 = summary( M0w2 )$coeff["Z",2]

    # Site weighted regression models
    M0w.site = lm( Yobs ~ Z + blk, weights=dat$weight.site, data=dat )
    tau.w.site = coef( M0w.site )[[2]]
    SE.w.site = summary( M0w.site )$coeff["Z",2]

    weightModels = data.frame( method=c("IPTW weighted regression (lm, naive)", "IPTW weighted regression (lm)", "IPTW weighted regression (lm, site)"),
                               tau = c( coef( M0w )[["Z"]], coef( M0w2 )[["Z"]], coef( M0w.site )[["Z"]] ),
                               SE = c( SE.w, SE.w2, SE.w.site ),
                               stringsAsFactors = FALSE )

    weightModels
}




#' Survey-weighted adjusted linear regression
#'
#' Use survey weight regression to reweight blocks to target unbiased ATE estimators.
#'
#' @param dat Dataframe to analyze, with Y, Z, and blk columns.
#' @return Dataframe of results for different estimators.
#' @importFrom survey svydesign svyglm
#' @importFrom stats gaussian
weighted.linear.estimators = function( dat ) {
    require( survey )
    dat$blk<-as.factor(dat$blk)

    Z.bar = mean( dat$Z )
    n = nrow( dat )
    J = length( unique( dat$blk ) )
    n.bar = n / J
    dat = dat %>% dplyr::group_by( blk ) %>%
        dplyr::mutate( p = mean( Z ),
                       nj = n(),
                       w.orig = ifelse( Z, 1/p, 1/(1-p) ),
                       weight = ifelse( Z, Z.bar / p, (1-Z.bar)/(1-p) ),
                       weight.site = weight * n.bar / nj ) %>%
        dplyr::ungroup()


    M0w = svyglm( Yobs ~ 0 + Z + blk,
                  design=svydesign(id=~1, weights=~w.orig, data=dat ),
                  family = gaussian() )
    SE.w = grab.SE( M0w )

    M0w2 = svyglm( Yobs ~ 0 + Z + blk,
                  design=svydesign(id=~1, weights=~weight, data=dat ),
                  family = gaussian() )
    SE.w2 = grab.SE( M0w2 )

    # Site weighted regression models
    M0w.site = svyglm( Yobs ~ 0 + Z + blk,
                   design=svydesign(id=~1, weights=~weight.site, data=dat ),
                   family = gaussian() )
    tau.w.site = coef( M0w.site )[["Z"]]
    SE.w.site = summary( M0w.site )$coeff["Z",2]

    weightModels = data.frame( method=c("IPTW weighted regression (naive)", "IPTW weighted regression", "IPTW weighted regression (site)"),
                               tau = c( coef( M0w )[["Z"]], coef( M0w2 )[["Z"]], coef( M0w.site )[["Z"]] ),
                               SE = c( SE.w, SE.w2, SE.w.site ),
                               stringsAsFactors = FALSE )

    weightModels
}


#' Interacted linear regression models
#'
#' These linear models have block by treatment interaction terms.  The final ATE
#' estimates are then weighted average of the block (site) specific ATE
#' estimates.
#'
#' SEs come from the overall variance-covariance matrix.
#'
#' @return Dataframe of the different versions of this estimator (person and site weighted)
interacted.linear.estimators = function( dat ) {
    require( multiwayvcov )
    require( tidyverse )

    J = length( unique( dat$blk ) )
    nj = table( dat$blk )
    n = nrow( dat )

    M0.int = lm( Yobs ~ 0 + Z * blk - Z, data = dat )
    VC = vcov( M0.int )

    tau.hats = coef(M0.int)[J + 1:J]

    wts = c( rep( 0, J ), rep( 1/J, J ) )
    tau.site = mean( tau.hats )
    SE.site = sqrt( t(wts) %*% VC %*% wts )

    wts.indiv = c( rep( 0, J ), nj/n )
    tau.indiv = weighted.mean( tau.hats, nj )
    SE.indiv = sqrt( t(wts.indiv) %*% VC %*% wts.indiv )

    interactModels = data.frame( method=c("fixed effects interact (site)", "fixed effects interact (indiv)"),
                                 tau = c( tau.site, tau.indiv ),
                                 SE = c( SE.site, SE.indiv ),
                                 stringsAsFactors = FALSE)
}


#' Estimate a series of linear models using different weighting schemes and
#' standard errors.
#'
#' @importFrom dplyr group_by ungroup mutate
#' @importFrom sandwich vcovHC vcovCL
#' @importFrom stats coef
#' @return Data frame of the various results.
#' @export
linear.model.estimators = function( Y, Z, B, data=NULL ) {
    if ( !is.null( data ) ) {
        dat = data
        names( dat )[1:3] = c( "Yobs", "Z", "blk" )
    }  else {
        dat = data.frame( Yobs = Y, Z = Z, blk = B )
    }

    dat$blk<-as.factor(dat$blk)

    FEmodels = fixed.effect.estimators( dat )

    weightModels = weighted.linear.estimators( dat )

    interactModels = interacted.linear.estimators( dat )

    # combine and return results
    bind_rows( FEmodels, weightModels, interactModels )
}



if ( FALSE ) {
    library( blkvar )
    library( tidyverse )
    dat = make.obs.data(p = 0.2)
    head( dat )
    fixed.effect.estimators( dat )

    weighted.linear.estimators.naive( dat )
    weighted.linear.estimators( dat )

    linear.model.estimators( dat$Yobs, dat$Z, dat$blk )
}



