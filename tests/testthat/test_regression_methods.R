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





