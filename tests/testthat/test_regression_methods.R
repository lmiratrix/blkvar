##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking regression estimators")


test_that("interacted regression runs on balanced data", {
    set.seed( 1019 )
    dat = make.obs.data.linear( method="big")
    head( dat )

    rs = blkvar:::interacted.linear.estimators( Yobs, Z, B, data=dat )

    rs
    expect_equal( rs$tau[[1]], rs$tau[[2]] )
    expect_equal( rs$SE[[1]], rs$SE[[2]] )

})







test_that("Weighted regression matches db", {
    set.seed( 10191010 )
    dat = make.obs.data.linear( method="small")
    #dat = make.obs.data( method="small")
    head( dat )

    sdat = calc.summary.stats( dat )
    sdat

    a = estimate.ATE.design.based.from.stats( sdat, weight="individual", method="finite" )
    b = estimate.ATE.design.based.from.stats( sdat, weight="site", method="finite" )
    a
    b

    head( dat )
    r = blkvar:::weighted.linear.estimators( Yobs, Z, B, data=dat )
    r
    expect_equal( a$tau.hat, r$tau[[1]] )
    expect_equal( b$tau.hat, r$tau[[2]] )


})


