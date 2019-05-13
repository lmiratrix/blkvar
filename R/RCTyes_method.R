


#' Implementation of the formula found in RCT-YES documentation.
#'
#' This can handle a superpopulation model, non-clustered but blocked.
#'
#' Taken from page 83 of Shochet RCT-YES paper (eq 6.25).
#'
#' The `superpop` variant is a modification of the original 'superpop.original',
#' pulling the weights from inside the squared term to outside. This method was
#' suggested in personal correspondance with Schochet.  If the weights are not
#' all 1, this can make a difference.
#'
#' @param sum_tab Table of summary statistics by block, from, e.g.,
#'   `block.data()`
#' @param siteID Vector of site IDs if there are randomization blocks nested in
#'    site that should be aggregated (will change results for site weighting only).
#' @param weight Individual weight (i.e., number of individuals in each block)
#'   or site weight (average site estimates (which will be considered block
#'   estimates if siteID is null)).
#' @param method finite, superpop, or superpop2 to give SEs that either capture
#'   uncertainty due to targeting a superpopulation quantity or not.
#' @return dataframe with calculated impacts and standard errors.
#' @export
estimate.ATE.design.based = function( sum_tab, siteID = NULL,
                    method = c( "finite", "superpop", "superpop.original" ),
                    weight = c( "individual", "site" ) ) {

    stopifnot( is.data.frame( sum_tab ) )
    stopifnot( all( c("Ybar0","Ybar1","n", "var1","var0", "n1","n0" ) %in% names( sum_tab ) ) )

    if ( !( "Ybar1" %in% names( sum_tab ) ) ) {
        sum_tab = convert.table.names( sum_tab )
    }

    method = match.arg( method )
    weight = match.arg( weight )
    h = nrow( sum_tab )

    # Get our block weights depending on target estimand
    if ( weight == "individual" ) {
        sum_tab$.weight = sum_tab$n
    } else {
        if ( !is.null( siteID ) ) {
            sum_tab = sum_tab %>% dplyr::group_by_( siteID ) %>%
                dplyr::mutate( .weight = n / sum( n ) ) %>% ungroup()
        } else {
            sum_tab$.weight = rep( 1, h )
        }
    }

    # calculate individual block treatment impact estimates
    sum_tab = mutate( sum_tab, tau.hat.b = Ybar1 - Ybar0 )

    # calculate overall ATE estimate by taking a weighted average of the
    # individual
    tau.hat = with( sum_tab, sum( tau.hat.b * .weight ) / sum( .weight ) )

    # Now do the SEs.
    if ( method == "finite" ) {
        # finite pop (Neyman)
        w.tot = sum( sum_tab$.weight )

        # calculate SEs for each block by itself
        sum_tab = mutate( sum_tab,
                          block.vars = (var1 / n1 + var0 / n0 ) )

        # and then take a weighted sum of these
        var = with( sum_tab, sum ( .weight^2 * block.vars ) / w.tot^2 )

        SE = sqrt( var )
    } else {  # superpopulation!

        # First aggregate to get sites, if needed
        if ( !is.null( siteID ) ) {
            sum_tab = sum_tab %>% group_by_( siteID ) %>%
                summarise( tau.hat.b = sum( tau.hat.b * .weight ) / sum( .weight ),
                           .weight = sum( .weight ) )
            h = nrow( sum_tab )
        }

        # Calculate average weight across sites
        wbar = mean( sum_tab$.weight )

        if ( method == "superpop.original" ) {
            # This is the formula 6.25
            asyVar = with( sum_tab, sum( (.weight * tau.hat.b - wbar * tau.hat)^2 ) / ((h-1)*h * wbar^2 ) )
        } else if ( method == "superpop" ) {
            # This is based on the email chain with Weiss, Pashley, etc.
            asyVar = with( sum_tab, sum( .weight^2 * (tau.hat.b - tau.hat)^2 ) / ((h-1)*h * wbar^2 ) )
        }

        SE = sqrt( asyVar )

    }

    data.frame( tau.hat = tau.hat, SE = SE, weight=weight, method=method,
                stringsAsFactors = FALSE )
}







# testing
if ( FALSE ) {
    library( tidyverse )

    dat = make.obs.data.linear( method="big")
    #dat = make.obs.data( method="small")
    head( dat )
    #write_csv( dat, path="some_fake_data.csv" )

    sdat = calc.summary.stats( dat )
    sdat
    #write_csv( sdat, path="summary_fake_data.csv" )

    estimate.design.based( sdat, weight="individual", method="finite" )
    estimate.design.based( sdat, weight="site", method="finite" )
    estimate.design.based( sdat, weight="individual", method="superpop" )
    estimate.design.based( sdat, weight="site", method="superpop" )

    compare_methods( dat$Yobs, dat$Z, dat$B )
}
