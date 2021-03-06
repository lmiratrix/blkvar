##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking RCT YES estimators")


test_that("Design Based (RCT Yes) functions work", {
    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( method="big")
    #dat = generate_blocked_data_obs( method="small")
    head( dat )

    sdat = calc_summary_stats( dat )
    sdat

    a = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="finite" )
    a2 = estimate_ATE_design_based_from_stats( sdat, weight="tx", method="finite" )
    expect_equal( a$ATE_hat, a2$ATE_hat )
    b = estimate_ATE_design_based_from_stats( sdat, weight="site", method="finite" )
    c = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="superpop" )
    d = estimate_ATE_design_based_from_stats( sdat, weight="site", method="superpop" )

    expect_equal( a$ATE_hat, c$ATE_hat )
    expect_equal( b$ATE_hat, d$ATE_hat )

    dat = generate_blocked_data_obs( n_k = c(6, 9, 12), p = c( 1/2, 1/3, 1/4 ) )
    sdat = calc_summary_stats( dat )
    a = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="finite" )
    a2 = estimate_ATE_design_based_from_stats( sdat, weight="tx", method="finite" )
    expect_true( a$ATE_hat != a2$ATE_hat )

})





