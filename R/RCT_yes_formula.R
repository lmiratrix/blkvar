


#' Implementation of the formula found in RCT-YES documentation.
#'
#' This can handle a superpopulation model, non-clustered but blocked.
#'
#' Taken from page 83 of Shochet RCT-YES paper (eq 6.25).
#'
#' The `superpop2` variant is a modification of the original, pulling the
#' weights from inside the squared term to outside.  If the weights are not all
#' 1, this can make a difference.
#'
#' @param sum_tab Table of summary statistics by block, from, e.g.,
#'   `block.data()`
#' @param weight Individual weight (i.e., number of individuals in each block)
#'   or site weight (average block estimates).
#' @param method finite, superpop, or superpop2 to give SEs that either capture
#'   uncertainty due to targeting a superpopulation quantity or not.
#' @export
calc.RCT.Yes.SE = function( sum_tab,
                    method = c( "finite", "superpop", "superpop.original" ),
                    weight = c( "individual", "site" ) ) {

    if ( !( "Y.bar.1" %in% names( sum_tab ) ) ) {
        sum_tab = convert.table.names( sum_tab )
    }

    method = match.arg( method )
    weight = match.arg( weight )
    h = nrow( sum_tab )

    # Get our block weights depending on target estimand
    if ( weight == "individual" ) {
        w = sum_tab$n
    } else {
        w = rep(1, h )
    }

    # calculate individual block treatment impact estimates
    tau.hat.b = with( sum_tab, Y.bar.1 - Y.bar.0 )

    # calculate overall ATE estimate by taking a weighted average of the
    # individual
    tau.hat = sum( tau.hat.b * w ) / sum( w )

    if ( method == "superpop.original" ) {
        wbar = mean( w )

        # This is the formula 6.25
        asyVar = sum( (w * tau.hat.b - wbar * tau.hat)^2 ) / ((h-1)*h * wbar^2 )
        SE = sqrt( asyVar)
    } else if ( method == "superpop" ) {
        # This is based on the email chain with Weiss, Pashley, etc.
        wbar = mean( w )

        asyVar = sum( w^2 * (tau.hat.b - tau.hat)^2 ) / ((h-1)*h * wbar^2 )
        SE = sqrt( asyVar)
    } else {
        w.tot = sum( w )
        # calculate SEs for each block by itself
        block.vars = with( sum_tab, (var.1 / n.1 + var.0 / n.0 ) )

        # and then take a weighted sum of these
        var =  sum ( w^2 * block.vars ) / w.tot^2

        SE = sqrt( var )
    }
    data.frame( tau.hat = tau.hat, SE = SE, weight=weight, method=method,
                stringsAsFactors = FALSE )
}







# testing
if ( FALSE ) {
    library( tidyverse )

    dat = make.obs.data( method="big")
    #dat = make.obs.data( method="small")
    head( dat )
    write_csv( dat, path="some_fake_data.csv" )

    sdat = calc.summary.stats( dat )
    sdat
    write_csv( sdat, path="summary_fake_data.csv" )

    calc.RCT.Yes.SE( sdat, weight="individual", method="finite" )
    calc.RCT.Yes.SE( sdat, weight="site", method="finite" )
    calc.RCT.Yes.SE( sdat, weight="individual", method="superpop" )
    calc.RCT.Yes.SE( sdat, weight="site", method="superpop" )

    compare.methods( dat$Yobs, dat$Z, dat$blk )
}
