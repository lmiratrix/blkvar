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

    a = estimate.ATE.design.based( sdat, weight="individual", method="finite" )
    b = estimate.ATE.design.based( sdat, weight="site", method="finite" )
    c = estimate.ATE.design.based( sdat, weight="individual", method="superpop" )
    d = estimate.ATE.design.based( sdat, weight="site", method="superpop" )

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


test_that("RCT Yes functions work with nested randomization blocks", {
    set.seed( 1019 )
    dat = make.obs.data.linear( X=1:50, method="big" )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )

    sdat = calc.summary.stats( dat, siteID="siteNo" )
    sdat

    a = estimate.ATE.design.based( sdat, siteID="siteID", weight="individual", method="finite" )
    b = estimate.ATE.design.based( sdat, siteID="siteID", weight="site", method="finite" )

    a2 = estimate.ATE.design.based( sdat, weight="individual", method="finite" )
    expect_equal( a, a2 )

    b2 = estimate.ATE.design.based( sdat, weight="site", method="finite" )
    b2
    b
    expect_true( b$tau.hat != b2$tau.hat )


    a = make.data( c( 4, 4, 4 ), tau = c( 10, 20, 30 ), exact=TRUE )
    a$sssite = c( 1, 1, 1, 1, 1, 1, 1, 1, 2, 2 , 2, 2)
    a$Z = 1
    a$Yobs = a$Y1
    b = a
    b$Z = 0
    b$Yobs = b$Y0
    a = bind_rows( a, b )
    a
    a %>% group_by( B ) %>% summarise( tau = mean( Y1 ) - mean( Y0 ),
                                       ybar1 = mean( Y1 ),
                                       ybar0 = mean( Y0 ) )
    aa = calc.summary.stats( a, siteID="sssite" )
    aa


    expect_equal( 20, estimate.ATE.design.based( aa, weight="site", method="finite"  )$tau.hat )
    expect_equal( (10+20)/4 + 15, estimate.ATE.design.based( aa, siteID="siteID", weight="site", method="finite"  )$tau.hat )


})



