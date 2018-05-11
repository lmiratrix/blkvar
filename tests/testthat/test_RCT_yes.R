##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking RCT YES estimators")


test_that("RCT Yes functions work", {
    set.seed( 1019 )
    dat = make.obs.data.linear( method="big")
    #dat = make.obs.data( method="small")
    head( dat )

    sdat = calc.summary.stats( dat )
    sdat

    a = calc.RCT.Yes.SE( sdat, weight="individual", method="finite" )
    b = calc.RCT.Yes.SE( sdat, weight="site", method="finite" )
    c = calc.RCT.Yes.SE( sdat, weight="individual", method="superpop" )
    d = calc.RCT.Yes.SE( sdat, weight="site", method="superpop" )

    expect_equal( a$tau.hat, c$tau.hat )
    expect_equal( b$tau.hat, d$tau.hat )

    scat = function( str, ... ) {
        cat( sprintf( str, ... ) )
    }
   # scat( "A se = %.2f\nC se = %.2f\n", a$SE, c$SE )
   # expect_true( a$SE <= c$SE )
   # expect_true( b$SE <= d$SE )

    sdat = block.data( dat$Yobs, B = dat$blk, Z=dat$Z )
    sdat
    dt = blkvar:::convert.table.names( sdat )
    dt
    d2 = calc.RCT.Yes.SE( dt, weight="site", method="superpop" )
    d2
    expect_equal( d, d2 )

    # checking auto-conversion of table names
    d3 = calc.RCT.Yes.SE( sdat, weight="site", method="superpop" )
    expect_equal( d, d3 )
})




