

#' Calculate ATE from block-level summary statistics.
#'
#' This method is an implementation of the formula found in RCT-YES
#' documentation.
#'
#' This can handle a superpopulation model, non-clustered but blocked.
#'
#' Formula used is taken from page 83 of Shochet RCT-YES paper (eq 6.25). The
#' `superpop` variant is a modification of the original 'superpop.original',
#' pulling the weights from inside the squared term to outside. This method was
#' suggested in personal correspondance with Schochet.  If the weights are not
#' all 1, this can make a difference.
#'
#' @param sum_tab Table of summary statistics by block, from, e.g.,
#'   `block.data()`
#' @inheritParams estimate_ATE_design_based
#' @param weight Individual weight (i.e., number of individuals in each block),
#'   site weight (average site estimates, which will be considered block
#'   estimates if siteID is null), tx weight (i.e., number of treated
#'   individuals in each block), or "passed" (with weight_col being a string
#'   name of what column has the weights to use).
#' @param weight_col Name of column that holds the weights to use for each block
#'   (NULL default).
#' @param method finite, superpop, or superpop2 to give SEs that either capture
#'   uncertainty due to targeting a superpopulation quantity or not.
#' @return dataframe with calculated impacts and standard errors.
#' @importFrom magrittr "%>%"
#' @export

estimate_ATE_design_based_from_stats <- function(sum_tab,
                                                 siteID = NULL,
                                                 method = c( "finite", "superpop", "superpop.original"),
                                                 weight = c("individual", "site", "tx", "passed"),
                                                 weight_col = NULL ) {
    ATE_hat <- NA
    stopifnot(is.data.frame(sum_tab))
    stopifnot(all(c("Ybar0","Ybar1","n", "var1","var0", "n1","n0" ) %in% names( sum_tab)))
    method <- match.arg(method)
    weight <- match.arg(weight)
    h <- nrow(sum_tab)
    if ( !is.null( weight_col ) && weight != "passed" ) {
        warning("Non-null weight_col typically means you have to specify 'weight=passed'" )
    }
    # Get our block weights depending on target estimand
    if (weight == "passed" ) {
        sum_tab$.weight = sum_tab[[weight_col]]
    } else if (weight == "individual") {
        sum_tab$.weight = sum_tab$n
    } else if (weight == "tx" ) {
        sum_tab$.weight = sum_tab$n1
    } else {
        if (!is.null(siteID)) {
            sum_tab <- sum_tab %>% dplyr::group_by(!!as.name(siteID)) %>%
                dplyr::mutate(.weight = n / sum(n)) %>% ungroup()
        } else {
            sum_tab$.weight <- rep(1, h)
        }
    }

    # calculate individual block treatment impact estimates
    sum_tab <- mutate( sum_tab, ATE_hat_b = Ybar1 - Ybar0)

    # calculate overall ATE estimate by taking a weighted average of the
    # individual
    ATE_hat <- with(sum_tab, sum(ATE_hat_b * .weight) / sum(.weight))
    # Now do the SEs.
    if (method == "finite") {
        # finite pop (Neyman)
        w.tot <- sum(sum_tab$.weight)
        # calculate SEs for each block by itself
        sum_tab <- mutate(sum_tab, block.vars = (var1 / n1 + var0 / n0 ))
        # and then take a weighted sum of these
        var <- sum (sum_tab$.weight ^ 2 * sum_tab$block.vars) / w.tot ^ 2
        SE <- sqrt(var)
    } else {  # superpopulation!
        # First aggregate to get sites, if needed
        if (!is.null(siteID)) {
            sum_tab <- sum_tab %>% group_by(!!as.name(siteID)) %>%
                dplyr::summarise(ATE_hat_b = sum(.data$ATE_hat_b * .data$.weight) / sum(.data$.weight),
                                 .weight = sum(.data$.weight))
            h <- nrow( sum_tab )
        }
        # Calculate average weight across sites
        wbar <- mean(sum_tab$.weight)
        if (method == "superpop.original") {
            # This is the formula 6.25
            asyVar <- sum((sum_tab$.weight * sum_tab$ATE_hat_b - wbar * ATE_hat) ^ 2) / ((h - 1) * h * wbar ^ 2)
        } else if (method == "superpop") {
            # This is based on the email chain with Weiss, Pashley, etc.
            asyVar <- sum(sum_tab$.weight ^ 2 * (sum_tab$ATE_hat_b - ATE_hat) ^ 2) / ((h -1 ) * h * wbar ^ 2)
        }
        SE <- sqrt(asyVar)
    }
    data.frame(ATE_hat = ATE_hat, SE = SE, weight = weight, method = method, stringsAsFactors = FALSE)
}




#' Implementation of the formula found in RCT-YES documentation.
#'
#' This can handle a superpopulation model, non-clustered but blocked.
#'
#' Taken from page 83 of Shochet RCT-YES paper (eq 6.25).
#'
#'The adjusted version, i.e., if control formula is passed, also uses
#'formula from the Schochet RCT Yes technical document.
#'
#'The `superpop` variant is a modification of the original 'superpop.original',
#'pulling the weights from inside the squared term to outside. This method was
#'suggested in personal correspondance with Schochet.  If the weights are not
#'all 1, this can make a difference.
#'
#'@param formula Input formula for analysis
#'@param control_formula What variables to control for, in the form of "~ X1 +
#'  X2".
#'@param data Dataframe with defined Yobs, Z, and B variables.
#'@param method finite, superpop, or superpop2 to give SEs that either capture
#'  uncertainty due to targeting a superpopulation quantity or not.
#'@param siteID Vector of site IDs if there are randomization blocks nested in
#'  site that should be aggregated (will change results for site weighting
#'  only).
#'@return dataframe with calculated impacts and standard errors.
#'@export

estimate_ATE_design_based <- function(formula,
                                      control_formula = NULL,
                                      data,
                                      siteID = NULL,
                                      method = c("finite", "superpop", "superpop.original"),
                                      weight = c("individual", "site", "tx")) {
    if (is.null(control_formula)) {
        data_table <- calc_summary_stats_formula(formula, data=data, siteID=siteID, add.neyman = TRUE)
        estimate_ATE_design_based_from_stats( data_table, siteID = siteID, method = method, weight = weight)
    } else {
        estimate_ATE_design_based_adjusted(formula = formula,
                                           control_formula = control_formula,
                                           siteID = siteID, data = data,
                                           method = method, weight = weight)
    }
}





#'@describeIn estimate_ATE_design_based  This directly implements the adjusted.  The main method will dispatch to this one if control_formula is not NULL.
#'@importFrom rlang .data
#'@export
estimate_ATE_design_based_adjusted <- function(formula,
                                               control_formula, data, siteID = NULL,
                                               method = c("finite", "superpop", "superpop.adj"),
                                               weight = c("individual", "site", "tx")) {
    stopifnot(!is.null(control_formula))

    # Determine which of the 4 versions of estimator we are doing.
    method <- match.arg(method)
    weight <- match.arg(weight)
    data <- make_canonical_data(formula, control_formula, siteID, data)
    # Get control variables
    c.names <- formula.tools::rhs.vars(control_formula)
    # Count observations, etc.
    N <- nrow(data)  # number of observations
    h <- length( unique( data$B ) ) # number of sites/clusters
    v <- length(c.names) # number of covariates

    # Center the X variables around their overall grand means
    center <- function(x) {
        x - mean(x)
    }
    data <- data %>% dplyr::group_by(B) %>% dplyr::mutate_at(c.names, center) %>% dplyr::ungroup()
    new.form <- make_FE_int_formula( control_formula = control_formula, data = data)
    # Fit a linear model with indicators for each site, and each site by treatment
    # interaction.  Also have covariates entered in not interacted with treatment.
    M0 <- lm( new.form, data = data )
    data$resid <- resid( M0 )

    # Aggregate data by group, calculating various summary statistics.
    # Note: Not all of these stats are used by all the estimators.
    # Key: The MSE.T and MSE.C are from the Schochet formulas.
    sum_tab <- data %>% group_by( B, siteID ) %>%
        dplyr::summarise(n = n(),
                         nT = sum(Z),
                         nC = n - nT,
                         Ybar.C = mean(Yobs[Z==0]),
                         Ybar.T = mean(Yobs[Z==1]),
                         MSE.T = sum(resid[Z==1] ^ 2) / ((N - v) * (nT / N) - 1),
                         MSE.C = sum(resid[Z == 0] ^ 2) / ((N - v) * (nC / N) - 1))
    #                      X1.bar = mean( X ),
    #                      X2.bar = mean( X2 ) )
    # make sure we have one row of stats per randomization block.
    stopifnot(length(unique(data$B)) == nrow(sum_tab))

    # Copy over our block-level impact estimates
    sum_tab$ATE_hat_b <- coef(M0)[(h + v + 1):(2 * h + v)]

    # Get our block weights depending on target estimand
    if (weight == "tx" ) {
        sum_tab$.weight <- sum_tab$nT
    } else if (weight == "individual") {
        sum_tab$.weight <- sum_tab$n
    } else {
        if (!is.null(siteID)) {
            sum_tab <- sum_tab %>%
                dplyr::group_by( !!as.name( siteID ) ) %>%
                dplyr::mutate( .weight = n / sum( n ) ) %>%
                ungroup()
        } else {
            sum_tab$.weight <- rep(1, h)
        }
    }

    ATE_hat <- NA
    SE <- NA

    # calculate overall ATE estimate by taking a weighted average of the
    # estimated block effects.
    ATE_hat <- sum( sum_tab$ATE_hat_b * sum_tab$.weight ) / sum( sum_tab$.weight )

    if (method == "finite") {
        # Finite population SE calculation
        # finite pop (Neyman)
        w.tot <- sum(sum_tab$.weight)

        # calculate SEs for each block by itself
        sum_tab = dplyr::mutate(sum_tab, block.var = MSE.T / nT + MSE.C / nC )

        # and then take a weighted sum of these
        var <- with(sum_tab, sum(.weight ^ 2 * block.var) / w.tot ^ 2)
        SE <- sqrt(var)
    } else {
        # Superpopulation SE calculation
        if (method == "superpop.adj") {
            # Center our averaged covariates around their site-weighted means.
            X.w <- weighted.mean(sum_tab$X1.bar, w = sum_tab$.weight)
            X2.w <- weighted.mean(sum_tab$X2.bar, w = sum_tab$.weight)
            sum_tab = mutate(sum_tab, X1c = X1.bar - X.w, X2c = X2.bar - X2.w)

            # Our second stage regression
            M.sp <- lm( ATE_hat_b ~ 1 + X1c + X2c, weights = .weight, data = sum_tab)
            # The intercept of this model should be our ATE_hat
            stopifnot(abs(coef(M.sp)[[1]] - ATE_hat) < 0.00001)

            # Now calculate the Superpopoulation SE
            sum_tab$ATE.pred <- predict(M.sp)

            # Equation 6.28, page 86 (using the residuals)
            w.bar <- mean(sum_tab$.weight)
            SE.SP <- sum(sum_tab$.weight ^ 2 * (sum_tab$ATE_hat_b - sum_tab$ATE.pred) ^ 2) / ((h - v - 1) * h * w.bar ^ 2)
            SE <- sqrt(SE.SP)
        } else {
            # First aggregate to get sites, if needed
            if (!is.null(siteID)) {
                sum_tab <- sum_tab %>% group_by(!!as.name(siteID)) %>%
                    dplyr::summarise( ATE_hat_b = sum(.data$ATE_hat_b * .data$.weight) / sum(.data$.weight),
                                      .weight = sum(.data$.weight))
                h <- nrow(sum_tab)
            }
            w.bar <- mean(sum_tab$.weight)
            SE.SP <- sum(sum_tab$.weight ^ 2 * (sum_tab$ATE_hat_b - ATE_hat) ^ 2) / ((h - 1) * h * w.bar ^ 2)
            SE <- sqrt(SE.SP)
        }
    }
    data.frame(ATE_hat = ATE_hat, SE = SE, weight = weight, method = method, stringsAsFactors = FALSE)
} # end estimator function






# # testing
# if ( FALSE ) {

    # dat = generate_blocked_data_obs_linear( method="big")
    # #dat = generate_blocked_data_obs( method="small")
    # head( dat )
    # #write_csv( dat, path="some_fake_data.csv" )

    # sdat = calc_summary_stats( dat )
    # sdat
    # #write_csv( sdat, path="summary_fake_data.csv" )

    # estimate.design.based( sdat, weight="individual", method="finite" )
    # estimate.design.based( sdat, weight="site", method="finite" )
    # estimate.design.based( sdat, weight="individual", method="superpop" )
    # estimate.design.based( sdat, weight="site", method="superpop" )

    # compare_methods( dat$Yobs, dat$Z, dat$B )
# }


# R implementation of methods described in Schochet paper
# Miratrix, 2019



all_adjusted_estimators <- function(formula, control_formula, data) {

  est <- ATE_hat <- NA
  ests <- expand.grid(method = c("finite", "superpop" ), #, "superpop.adj" ),
    weight = c( "individual", "site" ), stringsAsFactors = FALSE )
  ests <- tibble::as_tibble(ests)
  ests$est <- purrr::pmap( ests, estimate_ATE_design_based_adjusted,
                           formula = formula,
                           control_formula = control_formula,
                           data=data)
  ests$method <- NULL
  ests$weight <- NULL
  ests <- tidyr::unnest( ests, cols = est)
  ests$method <- paste0( "DBadj-", ests$method, "-", ests$weight )
  ests <- dplyr::rename( ests, ATE_hat = ATE_hat )
  ests
}

##### Test code on single dataset (adjusted versions) ######
# if ( FALSE ) {

    # #set.seed( 1019)
    # dat = generate_multilevel_data( n.bar = 20, J = 5, beta.X = 0.5 )
    # nrow( dat )
    # head( dat )

    # # Add second X variable for kicks
    # dat$X2 = dat$W + rnorm( nrow( dat ) )

    # dat.bk = dat
    # dat = rename( dat,
                  # Tx = Z,
                  # Y = Yobs,
                  # X1 = X,
                  # ID = sid )
    # head( dat )

    # #dat = sample_n( dat, size=nrow(dat) )
    # #head( dat )

    # estimate_ATE_design_based_adjusted( formula = Y ~ Tx*ID,
                                        # control_formula = ~ X1 + X2,
                                        # data=dat,
                                        # method="superpop",
                                        # weight="individual" )


    # estimate_ATE_design_based( formula = Y ~ Tx*ID,
                               # data=dat,
                               # method="superpop",
                               # weight="individual" )


    # all_adjusted_estimators( formula = Y ~ Tx*ID,
                             # control_formula = ~ X1 + X2,
                             # data=dat )




    # # With more vars
    # dat$Xx = rnorm( nrow(dat) )
    # dat$Xalso = rnorm( nrow(dat) )
    # all_adjusted_estimators( formula = Y ~ Tx*ID,
                             # control_formula = ~ X1 + X2 + Xx + Xalso,
                             # data=dat )


    # rs = compare_methods( Y, Tx, ID, data=dat,
                          # include_MLM = FALSE )
    # rs


    # if ( FALSE ) {
        # source( "covariate_adjusted_design_cleaned.R" )
        # all_adjusted_estimators.old( data=dat.bk )

    # }



# }
