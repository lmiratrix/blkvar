##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking FIRC estimators")


test_that("FIRC functions work", {
    set.seed( 1019 )

    sdf = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv = analysis.FIRC( sdf )
    pv2 = estimate.ATE.FIRC(Yobs, Z, sid, data=sdf, include.testing = TRUE )

    expect_equal( pv, pv2$p.variation )

    pv.aov = analysis.FIRC( sdf, anova=TRUE )
    expect_equal( pv, pv.aov )

} )




test_that("RIRC functions work", {
    set.seed( 1019 )

    sdf = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv = estimate.ATE.RIRC( Yobs, Z, sid, data=sdf )

    pv = estimate.ATE.RIRC( Yobs, Z, sid, data=sdf, REML = TRUE, include.testing = FALSE )

    expect_error( pv = estimate.ATE.RIRC( Yobs, Z, sid, data=sdf, REML = TRUE, include.testing = TRUE ) )


} )




test_that("RICC functions work", {
    set.seed( 1019 )

    sdf = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv = estimate.ATE.RICC( Yobs, Z, sid, data=sdf )
    pv
    expect_true( is.na( pv$deviance ) )
    expect_true( pv$ATE > 0 )
} )
