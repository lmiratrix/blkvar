library( testthat )
library( blkvar )

context("Checking blocking functions (Nicole package parts)")

test_that("Dropping levels works", {

    set.seed( 1019 )
    dat = make.obs.data.linear( method="big")
    dat$Z[1:4] = 0

    expect_false( all( table( dat$Z, dat$B ) == 0 ) )

    expect_warning( rs <- calc.summary.stats( dat$Yobs, dat$Z, dat$B ) )

    expect_equal( nrow(rs), 3 )

})



