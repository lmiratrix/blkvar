##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking regression estimators")


test_that("interacted regression runs on balanced data", {
    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( method="big")
    head( dat )

    rs = blkvar:::interacted_linear_estimators( Yobs, Z, B, data=dat )

    rs
    expect_equal( rs$ATE_hat[[1]], rs$ATE_hat[[2]] )
    expect_equal( rs$SE[[1]], rs$SE[[2]] )

})







test_that("Weighted regression matches db", {
    set.seed( 10191010 )
    dat = generate_blocked_data_obs_linear( method="small")
    #dat = generate_blocked_data_obs( method="small")
    head( dat )
    table( dat$B )
    sdat = calc_summary_stats( dat )
    sdat

    a = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="finite" )
    b = estimate_ATE_design_based_from_stats( sdat, weight="site", method="finite" )
    a
    b
    expect_true( is.na( b$SE ) )

    head( dat )
    r = blkvar:::weighted_linear_estimators( Yobs ~ Z*B, data=dat )
    r
    expect_equal( a$ATE_hat, r$ATE_hat[[1]] )
    expect_equal( b$ATE_hat, r$ATE_hat[[2]] )


})






test_that("Different Weighted regression flags work", {
    set.seed( 1019101010 )
    dat = generate_blocked_data_obs_linear( method="small")
    head( dat )

    A = blkvar:::weighted_linear_estimators( Yobs ~ Z*B, data=dat )
    B = blkvar:::weighted_linear_estimators( Yobs ~ Z*B, data=dat,
                                             scaled.weights = FALSE )

    C = blkvar:::weighted_linear_estimators( Yobs ~ Z*B, data=dat,
                                             weight.method = "precision" )
    D = blkvar:::weighted_linear_estimators( Yobs ~ Z*B, data=dat,
                                             scaled.weights = FALSE,
                                             weight.method = "precision" )

    A
    B
    C
    D

    expect_equal( A$ATE_hat, B$ATE_hat )
    expect_equal( C$ATE_hat, D$ATE_hat )

    expect_equal( A$ATE_hat, C$ATE_hat )



})


