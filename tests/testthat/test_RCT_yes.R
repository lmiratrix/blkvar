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

    a = estimate.ATE.design.based.from.stats( sdat, weight="individual", method="finite" )
    b = estimate.ATE.design.based.from.stats( sdat, weight="site", method="finite" )
    c = estimate.ATE.design.based.from.stats( sdat, weight="individual", method="superpop" )
    d = estimate.ATE.design.based.from.stats( sdat, weight="site", method="superpop" )

    expect_equal( a$tau.hat, c$tau.hat )
    expect_equal( b$tau.hat, d$tau.hat )



})



test_that("calc summary stats with site works", {
    set.seed( 1019 )
    dat = make.obs.data.linear( X=1:50, method="big" )
    #dat = make.obs.data( method="small")
    nrow( dat )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )
    table( dat$siteNo )
    table( dat$B )

    sdat = calc.summary.stats( dat, siteID="siteNo" )
    sdat
    expect_true( "siteID" %in% names(sdat) )
    expect_true( all( sdat$siteID %in% unique( dat$siteNo ) ) )
} )





