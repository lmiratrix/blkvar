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


# See https://www.jepusto.com/handmade-clubsandwich/
clubsandwich.variance = function( w, tau.hat.b, tau.hat ) {
    W = sum( w )
    V = (1/W^2) * sum( ( w^2 * (tau.hat.b - tau.hat)^2 ) / (1 - w/W) )

    df.inv = sum( w^2 / (W - w)^2 ) - (2/W) *sum( w^3 / (W-w)^2 ) + (1/W^2)*sum( w^2/(W-w) )^2

    list( var.hat=V, df = 1/df.inv )
}


fixed.effect.estimators = function( Yobs, Z, B, siteID = NULL, data=NULL, block.stats = NULL ) {
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


    # make sites RA blocks if there are no sites.
    if ( is.null( data$siteID ) ) {
        data$siteID = data$B
    }

    # simple linear model
    M0 = lm( Yobs ~ 0 + Z + B, data=data )
    SE.lm = summary( M0 )$coeff["Z",2]

    # est ATE
    tau.hat = coef(M0)[["Z"]]

    # Huber-White SEs
    vcov_sand = sandwich::vcovHC(M0, type = "HC1")
    SE.lm.sand <-sqrt( vcov_sand[1,1] )

    # Cluster robust SEs (clustering at site level)
    vcov_clust = sandwich::vcovCL( M0, siteID )
    SE.lm.clust = sqrt( vcov_clust[1,1] )

    # Cluster robust SEs (clustering at site level using clubSandwich)
    # aggregate!
    if ( is.null( block.stats ) ) {
        block.stats = calc.summary.stats(Yobs, Z, B, data=data, siteID=siteID, add.neyman = FALSE )
    }
    block.stats = mutate( block.stats, tau.hat = Ybar1 - Ybar0,
                          prec = n * (n0/n) * (n1/n),
                          w = 1/prec )
    if ( !is.null( siteID ) ) {
        # aggregate blocks into sites and calculate site weights and tau.hats
        block.stats = block.stats %>% group_by( siteID ) %>%
            summarise( tau.hat = sum( w * tau.hat ) / sum( w ),
                       w = 1 / sum( prec ) )
    }
    cs.var = clubsandwich.variance( block.stats$w, block.stats$tau.hat, tau.hat )
    SE.lm.clust.club = sqrt( cs.var$var.hat )

    FEmodels = data.frame( method=c("FE", "FE (sand)", "FE (cluster)", "FE (club)" ),
                           tau = rep( tau.hat, 4 ),
                           SE = c( SE.lm, SE.lm.sand, SE.lm.clust, SE.lm.clust.club ),
                           stringsAsFactors = FALSE )

    FEmodels
}


#' Individual Weighted regression models using the lm and precision weighting
#'
#' WARNING: This is the wrong kind of regression.  The weights are correct, but
#' the use is wrong.
#'
#' NOTE: This method has not been updated for site
weighted.linear.estimators.naive = function( Yobs, Z, B, data=NULL ) {
    if ( missing( "Z" ) && is.null( data ) ) {
        data = Yobs
    }
    if ( is.null( data ) ) {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
    } else {
        if ( missing( "Z" ) ) {
            stopifnot( all( c( "Yobs", "Z", "B" ) %in% names(data) ) )
        } else {
            data = rename_( data,
                            Yobs = as.character( substitute( Yobs ) ),
                            Z = as.character( substitute( Z ) ),
                            B = as.character( substitute( B ) ) )
        }
    }
    dat = data

    Z.bar = mean( dat$Z )
    n = nrow( dat )
    J = length( unique( dat$B ) )
    n.bar = n / J
    dat = dat %>% dplyr::group_by( B ) %>%
        dplyr::mutate( p = mean( Z ),
                       nj = n(),
                       w.orig = ifelse( Z, 1/p, 1/(1-p) ),
                       weight = ifelse( Z, Z.bar / p, (1-Z.bar)/(1-p) ),
                       weight.site = weight * n.bar / nj ) %>%
        dplyr::ungroup()

    M0w = lm( Yobs ~ Z + B, weights=dat$w.orig, data=dat )
    SE.w = summary( M0w )$coeff["Z",2]

    M0w2 = lm( Yobs ~ Z + B, weights=dat$weight, data=dat )
    SE.w2 = summary( M0w2 )$coeff["Z",2]

    # Site weighted regression models
    M0w.site = lm( Yobs ~ Z + B, weights=dat$weight.site, data=dat )
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
#' @param data Dataframe to analyze, with Y, Z, and B columns.
#' @return Dataframe of results for different estimators.
#' @importFrom survey svydesign svyglm
#' @importFrom stats gaussian
weighted.linear.estimators = function( Yobs, Z, B, data=NULL, siteID = NULL ) {

    if ( missing( "Z" ) && is.null( data ) ) {
        data = Yobs
    }
    if ( is.null( data ) ) {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
    } else {
        if ( missing( "Z" ) ) {
            stopifnot( all( c( "Yobs", "Z", "B" ) %in% names(data) ) )
        } else {
            data = rename_( data,
                            Yobs = as.character( substitute( Yobs ) ),
                            Z = as.character( substitute( Z ) ),
                            B = as.character( substitute( B ) ) )
        }
    }
    dat = data

    require( survey )
    dat$B<-as.factor(dat$B)

    Z.bar = mean( dat$Z )
    n = nrow( dat )
    J = length( unique( dat$B ) )
    n.bar = n / J
    dat = dat %>% dplyr::group_by( B ) %>%
        dplyr::mutate( p = mean( Z ),
                       nj = n(),
                       w.orig = ifelse( Z, 1/p, 1/(1-p) ),
                       weight = ifelse( Z, Z.bar / p, (1-Z.bar)/(1-p) ),
                       weight.site = weight * n.bar / nj ) %>%
        dplyr::ungroup()

    if ( !is.null( siteID ) ) {
        dat$siteID = dat[[siteID]]

        n.site = length( unique( dat$siteID ) )
        n.bar = n / n.site

        # adjust weights to account for blocks nested within sites
        dat = dat %>% group_by( siteID ) %>%
            mutate( weight.site = weight * n.bar / n() )
    }

    M0w = svyglm( Yobs ~ 0 + Z + B,
                  design=svydesign(id=~1, weights=~w.orig, data=dat ),
                  family = gaussian() )
    SE.w = grab.SE( M0w )

    M0w2 = svyglm( Yobs ~ 0 + Z + B,
                  design=svydesign(id=~1, weights=~weight, data=dat ),
                  family = gaussian() )
    SE.w2 = grab.SE( M0w2 )

    # Site weighted regression models
    M0w.site = svyglm( Yobs ~ 0 + Z + B,
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
#' If siteID passed, it will weight the RA blocks within site and then average
#' these site estimates.
#'
#' SEs come from the overall variance-covariance matrix.
#'
#' @return Dataframe of the different versions of this estimator (person and
#'   site weighted)
interacted.linear.estimators = function( Yobs, Z, B, siteID = NULL, data=NULL ) {

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

    require( multiwayvcov )

    J = length( unique( data$B ) )
    nj = table( data$B )
    n = nrow( data )

    data$B = as.factor( data$B )
    M0.int = lm( Yobs ~ 0 + Z * B - Z, data=data )
    VC = vcov( M0.int )

    tau.hats = coef(M0.int)[J + 1:J]

    # site weighting
    if ( !is.null( data$siteID ) ) {
        # aggregate!
        wts = data %>% group_by( B, siteID ) %>%
            summarise( n = n() ) %>%
            group_by( siteID ) %>%
            mutate( wts = n / sum( n ) )

        # some checks to make sure we are matching RA blocks and sites to the right things
        stopifnot( nrow( wts ) == J )
        nms = gsub( "Z:B", "", names( tau.hats ) )
        stopifnot( all( nms == wts$B ) )

        wts = wts$wts / sum(wts$wts)
    } else {
        wts = rep( 1/J, J )
    }

    # the block SEs from our linear model
    SE.hat = diag( VC )[(J+1):(2*J)]

    tau.site = weighted.mean( tau.hats, wts )

    # Calculate SE for tau.site
    SE.site = sqrt( sum( wts^2 * SE.hat ) )

    wts.indiv = nj/n
    tau.indiv = weighted.mean( tau.hats, wts.indiv )
    SE.indiv = sqrt( sum( wts.indiv^2 * SE.hat ) )

    # This is the cautious way we don't need since we have 0s in the off diagonal
    # SE.site = sqrt( t(wts) %*% VC %*% wts )
    # faster way---this should work easily.
    #sqrt( t(wts.indiv) %*% VC %*% wts.indiv )

    interactModels = data.frame( method=c("FE interact (site)", "FE interact (indiv)"),
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
linear.model.estimators = function( Yobs, Z, B, siteID = NULL, data=NULL, block.stats = NULL ) {
    if ( missing( "Z" ) && is.null( data ) ) {
        data = Yobs
    }
    if ( is.null( data ) ) {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
    } else {
        if ( missing( "Z" ) ) {
            stopifnot( all( c( "Yobs", "Z", "B" ) %in% names(data) ) )
        } else {
            data = rename_( data,
                            Yobs = as.character( substitute( Yobs ) ),
                            Z = as.character( substitute( Z ) ),
                            B = as.character( substitute( B ) ) )
        }
    }
    dat = data
    dat$B<-as.factor(dat$B)

    FEmodels = fixed.effect.estimators( Yobs, Z, B, siteID = siteID, data=dat, block.stats = block.stats )

    weightModels = weighted.linear.estimators( Yobs, Z, B, data=dat )

    interactModels = interacted.linear.estimators( Yobs, Z, B, siteID = siteID, data=dat )

    # combine and return results
    bind_rows( FEmodels, weightModels, interactModels )
}



if ( FALSE ) {
    library( blkvar )
    library( dplyr )
    dat = make.obs.data(p = 0.2)
    head( dat )

    fixed.effect.estimators( Yobs, Z, blk, data=dat )

    #debug( weighted.linear.estimators.naive )
    weighted.linear.estimators.naive( Yobs, Z, blk, data=dat )

    dat2 = rename( dat, B= blk )
    head( dat2 )
    weighted.linear.estimators.naive( dat2 )

    weighted.linear.estimators( Yobs, Z, blk, data=dat )

    debug( linear.model.estimators)
    linear.model.estimators( dat$Yobs, dat$Z, dat$blk )
}



