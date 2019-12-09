



#' Summarise data by block.
#'
#' Given dataframe or list of variables, return table of stats for each
#' randomization block.
#'
#' @param dat Dataframe with defined Yobs, Z, and B variables.
#' @param siteID If not null, name of siteID that has randomization blocks
#'   nested inside.
#'
#' @param  formula  Formula of form Yobs ~ Z * B (ONLY)
#'
#' @return dataframe with summary statistics by block
#' @export
calc.summary.stats.formula = function( formula, data = NULL, siteID = NULL, add.neyman = FALSE ) {

    # Figure out the covariates we are using
    if(length(formula.tools::lhs.vars(formula)) != 1 | length(formula.tools::rhs.vars(formula)) != 2){
        stop("The formula argument must be of the form outcome ~ treatment*block_id.")
    }

    main.vars <- formula.tools::get.vars(formula, data=data)
    if(any(!(main.vars %in% colnames(data)))){
        browser()
        stop("Some variables in formula are not present in your data.")
    }

    out.name = formula.tools::lhs.vars(formula)[[1]]
    main.name = formula.tools::rhs.vars(formula)

    # Copy over the variables to our names
    data$Yobs = data[[out.name]]
    data$Z = data[[ main.name[[1]] ]]
    data$B = data[[ main.name[[2]] ]]
    if ( !is.null( siteID ) && (length( siteID ) == 1 ) ) {
        data$siteID = data[[ siteID ]]
    } else {
        data$siteID = siteID
    }

    calc.summary.stats( Yobs, Z, B, data=data, siteID = siteID, add.neyman = add.neyman )

}







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
estimate.ATE.design.based = function( formula,
                                      control.formula = NULL,
                                      siteID = NULL,
                                      data,
                                      method = c( "finite", "superpop", "superpop.original" ),
                                      weight = c( "individual", "site" ) ) {

    if ( is.null( control.formula ) ) {

        data.table<-calc.summary.stats.formula(formula, data=data, siteID=siteID, add.neyman = TRUE )
        estimate.ATE.design.based.from.stats( data.table, siteID = siteID, method=method, weight=weight )
    } else {
        estimate.ATE.design.based.adjusted( formula=formula,
                                            control.formula=control.formula,
                                            siteID = siteID,
                                            data=data,
                                            method = method,
                                            weight = weight )
    }
}


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
estimate.ATE.design.based.from.stats = function( sum_tab, siteID = NULL,
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
    library( dplyr )

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





# R implementation of methods described in Schochet paper
# Miratrix, 2019

library( formula.tools )


scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}




#' Adjusted design based estimator for ATE
#'
#'Given two covariates, X1 and X2, calculate adjusted ATE estimates for the
#'multisite trial.  This uses the formula from the Schochet RCT Yes technical
#'document.
#'
#'@param formula
#'@param control.formula
#'@param data
#'@param method
#'@param weight
#'
#'@return
#'@export
#'
estimate.ATE.design.based.adjusted = function( formula,
                                               control.formula,
                                               data,
                                               siteID = NULL,
                                               method = c( "finite", "superpop", "superpop.adj" ),
                                               weight = c( "individual", "site" ) ) {

    stopifnot( !is.null( control.formula ) )

    # Determine which of the 4 versions of estimator we are doing.
    method = match.arg( method )
    weight = match.arg( weight )

    # Figure out the covariates we are using
    if(length(formula.tools::lhs.vars(formula)) != 1 | length(formula.tools::rhs.vars(formula)) != 2){
        stop("The formula argument must be of the form outcome ~ treatment:block_id.")
    }
    if(length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1){
        stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
    }

    main.vars <- formula.tools::get.vars(formula, data=data)
    if(any(!(main.vars %in% colnames(data)))){
        browser()
        stop("Some variables in formula are not present in your data.")
    }
    control.vars <- formula.tools::get.vars(control.formula,data=data)
    if(any(!(control.vars %in% colnames(data)))){
        stop("Some variables in control.formula are not present in your data.")
    }

    out.name = formula.tools::lhs.vars(formula)[[1]]
    main.name = formula.tools::rhs.vars(formula)
    c.names = formula.tools::rhs.vars(control.formula)

    # Copy over the variables to our names
    data$Yobs = data[[out.name]]
    data$Z = data[[ main.name[[1]] ]]
    data$B = data[[ main.name[[2]] ]]

    if ( is.null( siteID ) ) {
        data$siteID = data$B
    } else {
        data$siteID = data[[siteID]]
    }

    # Count observations, etc.
    N = nrow( data )  # number of observations
    h = length( unique( data$B ) ) # number of sites/clusters
    v = length( c.names ) # number of covariates



    # Center the X variables around their overall grand means
    center = function(x) {
        x - mean(x)
    }
    data = data %>% group_by( B ) %>%
        mutate_at( c.names, center ) %>% ungroup()

    new.form = make.FE.int.formula( control.formula=control.formula, data=data )

    # Fit a linear model with indicators for each site, and each site by treatment
    # interaction.  Also have covariates entered in not interacted with treatment.
    M0 = lm( new.form, data = data )

    data$resid = resid( M0 )



    # Aggregate data by group, calculating various summary statistics.
    # Note: Not all of these stats are used by all the estimators.
    # Key: The MSE.T and MSE.C are from the Schochet formulas.
    sum_tab = data %>% group_by( B, siteID ) %>%
        dplyr::summarise( n = n(),
                          nT = sum( Z ),
                          nC = n - nT,
                          Ybar.C = mean( Yobs[Z==0] ),
                          Ybar.T = mean( Yobs[Z==1] ),
                          MSE.T = sum( resid[Z==1]^2 ) / ( (N-v)*(nT/N) - 1 ),
                          MSE.C = sum( resid[Z==0]^2 ) / ( (N-v)*(nC/N) - 1 ) )
    #                      X1.bar = mean( X ),
    #                      X2.bar = mean( X2 ) )

    # make sure we have one row of stats per randomization block.
    stopifnot( length( unique( data$B ) ) == nrow( sum_tab ) )

    # Copy over our block-level impact estimates
    sum_tab$tau.hat.b = coef( M0 )[(h+v+1):(2*h+v)]

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


    tau.hat = NA
    SE = NA

    # calculate overall ATE estimate by taking a weighted average of the
    # estimated block effects.
    tau.hat = with( sum_tab, sum( tau.hat.b * .weight ) / sum( .weight ) )
    tau.hat

    if ( method == "finite" ) {
        # Finite population SE calculation


        # finite pop (Neyman)
        w.tot = sum( sum_tab$.weight )

        # calculate SEs for each block by itself
        sum_tab = mutate( sum_tab,
                          block.var = MSE.T / nT + MSE.C / nC )

        # and then take a weighted sum of these
        var = with( sum_tab, sum ( .weight^2 * block.var ) / w.tot^2 )

        SE = sqrt( var )

    } else {
        # Superpopulation SE calculation

        if ( method =="superpop.adj" ) {
            # Center our averaged covariates around their site-weighted means.
            X.w = weighted.mean( sum_tab$X1.bar, w = sum_tab$.weight )
            X2.w = weighted.mean( sum_tab$X2.bar, w = sum_tab$.weight )

            sum_tab = mutate( sum_tab,
                              X1c = X1.bar - X.w,
                              X2c = X2.bar - X2.w )

            # Our second stage regression
            M.sp = lm( tau.hat.b ~ 1 + X1c + X2c, weights=.weight, data=sum_tab )

            # The intercept of this model should be our tau.hat
            stopifnot( abs( coef(M.sp)[[1]] - tau.hat ) < 0.00001 )

            # Now calculate the Superpopoulation SE
            sum_tab$tau.pred = predict( M.sp )

            # Equation 6.28, page 86 (using the residuals)
            w.bar = mean( sum_tab$.weight )
            SE.SP = with( sum_tab,
                          sum( .weight^2 * ( tau.hat.b - tau.pred )^2 ) / ( (h - v - 1)*h*w.bar^2 ) )
            SE = sqrt( SE.SP )
        } else {
            # First aggregate to get sites, if needed
            if ( !is.null( siteID ) ) {
                sum_tab = sum_tab %>% group_by_( siteID ) %>%
                    summarise( tau.hat.b = sum( tau.hat.b * .weight ) / sum( .weight ),
                               .weight = sum( .weight ) )
                h = nrow( sum_tab )
            }

            w.bar = mean( sum_tab$.weight )
            SE.SP = with( sum_tab,
                          sum( .weight^2 * ( tau.hat.b - tau.hat )^2 ) / ( (h - 1)*h*w.bar^2 ) )
            SE = sqrt( SE.SP )

        }
    }

    data.frame(  tau.hat = tau.hat,
                 SE = SE,
                 weight = weight,
                 method=method,
                 stringsAsFactors = FALSE)


} # end estimator function


all.adjusted.estimators = function( formula,
                                    control.formula,
                                    data ) {
    ests = expand.grid( method = c("finite", "superpop" ), #, "superpop.adj" ),
                        weight = c( "individual", "site" ),
                        stringsAsFactors = FALSE )
    ests = as_tibble(ests)

    ests$est = pmap( ests, estimate.ATE.design.based.adjusted,
                     formula = formula,
                     control.formula = control.formula,
                     data=data )

    ests$method = NULL
    ests$weight = NULL
    ests = unnest( ests, c(est) )
    ests$method = paste0( "DBadj-", ests$method, "-", ests$weight )
    ests = rename( ests, tau = tau.hat )
    ests
}

##### Test code on single dataset (adjusted versions) ######

if ( FALSE ) {


    library( blkvar )
    library( tidyverse )


    #set.seed( 1019)
    dat = gen.dat( n.bar = 20, J = 5, beta.X = 0.5 )
    nrow( dat )
    head( dat )

    # Add second X variable for kicks
    dat$X2 = dat$W + rnorm( nrow( dat ) )

    dat.bk = dat
    dat = rename( dat,
                  Tx = Z,
                  Y = Yobs,
                  X1 = X,
                  ID = sid )
    head( dat )

    #dat = sample_n( dat, size=nrow(dat) )
    #head( dat )

    estimate.ATE.design.based.adjusted( formula = Y ~ Tx*ID,
                                        control.formula = ~ X1 + X2,
                                        data=dat,
                                        method="superpop",
                                        weight="individual" )


    estimate.ATE.design.based( formula = Y ~ Tx*ID,
                               data=dat,
                               method="superpop",
                               weight="individual" )


    all.adjusted.estimators( formula = Y ~ Tx*ID,
                             control.formula = ~ X1 + X2,
                             data=dat )




    # With more vars
    dat$Xx = rnorm( nrow(dat) )
    dat$Xalso = rnorm( nrow(dat) )
    all.adjusted.estimators( formula = Y ~ Tx*ID,
                             control.formula = ~ X1 + X2 + Xx + Xalso,
                             data=dat )


    rs = compare_methods( Y, Tx, ID, data=dat,
                          include.MLM = FALSE )
    rs


    if ( FALSE ) {
        source( "covariate_adjusted_design_cleaned.R" )
        all.adjusted.estimators.old( data=dat.bk )

    }



}




