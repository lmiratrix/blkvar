##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking RCT YES estimators")


test_that("Design Based (RCT Yes) functions work", {
    set.seed( 1019 )
    dat = make_obs_data_linear( method="big")
    #dat = make_obs_data( method="small")
    head( dat )

    sdat = calc_summary_stats( dat )
    sdat

    a = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="finite" )
    b = estimate_ATE_design_based_from_stats( sdat, weight="site", method="finite" )
    c = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="superpop" )
    d = estimate_ATE_design_based_from_stats( sdat, weight="site", method="superpop" )

    expect_equal( a$tau_hat, c$tau_hat )
    expect_equal( b$tau_hat, d$tau_hat )



})





