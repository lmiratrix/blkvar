# # testing
# if ( FALSE ) {

    # dat = make_obs_data_linear( method="big")
    # #dat = make_obs_data( method="small")
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



all_adjusted_estimators <- function(formula, control.formula, data) {

  est <- tau_hat <- NA
  ests <- expand.grid(method = c("finite", "superpop" ), #, "superpop.adj" ),
    weight = c( "individual", "site" ), stringsAsFactors = FALSE )
  ests <- tibble::as_tibble(ests)
  ests$est <- purrr::pmap( ests, estimate_ATE_design_based_adjusted, formula = formula, control.formula = control.formula, data=data)
  ests$method <- NULL
  ests$weight <- NULL
  ests <- tidyr::unnest( ests, cols = est)
  ests$method <- paste0( "DBadj-", ests$method, "-", ests$weight )
  ests <- dplyr::rename( ests, tau = tau_hat )
  ests
}

##### Test code on single dataset (adjusted versions) ######
# if ( FALSE ) {

    # #set.seed( 1019)
    # dat = gen_dat( n.bar = 20, J = 5, beta.X = 0.5 )
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
                                        # control.formula = ~ X1 + X2,
                                        # data=dat,
                                        # method="superpop",
                                        # weight="individual" )


    # estimate_ATE_design_based( formula = Y ~ Tx*ID,
                               # data=dat,
                               # method="superpop",
                               # weight="individual" )


    # all_adjusted_estimators( formula = Y ~ Tx*ID,
                             # control.formula = ~ X1 + X2,
                             # data=dat )




    # # With more vars
    # dat$Xx = rnorm( nrow(dat) )
    # dat$Xalso = rnorm( nrow(dat) )
    # all_adjusted_estimators( formula = Y ~ Tx*ID,
                             # control.formula = ~ X1 + X2 + Xx + Xalso,
                             # data=dat )


    # rs = compare_methods( Y, Tx, ID, data=dat,
                          # include_MLM = FALSE )
    # rs


    # if ( FALSE ) {
        # source( "covariate_adjusted_design_cleaned.R" )
        # all_adjusted_estimators.old( data=dat.bk )

    # }



# }
