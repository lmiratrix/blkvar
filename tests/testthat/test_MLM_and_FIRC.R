##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking FIRC estimators")


test_that("FIRC functions work", {
    set.seed( 1019 )

    sdf = generate_multilevel_data( n.bar=10, J=10,
                   rho2.0W = 0.3, rho2.1W = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv2 = estimate_ATE_FIRC(Yobs, Z, sid, data=sdf, include_testing = TRUE )


    pv.aov = estimate_ATE_FIRC( Yobs, Z, sid, data=sdf, include_testing = TRUE, anova=TRUE )
    expect_equal( pv2$p_variation, pv.aov$p_variation )


    # Check arb variable names works.
    sdf2 = rename( sdf, YY = Yobs, ZZ = Z, S.I.D = sid )

    pv2.B = estimate_ATE_FIRC(YY, ZZ, S.I.D, data=sdf2, include_testing = TRUE )

    expect_equal( pv2.B, pv2 )

} )





test_that("RIRC functions work", {
    set.seed( 1019 )

    sdf = generate_multilevel_data( n.bar=10, J=10,
                   rho2.0W = 0.3, rho2.1W = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv = estimate_ATE_RIRC( Yobs, Z, sid, data=sdf )

    pv = estimate_ATE_RIRC( Yobs, Z, sid, data=sdf, REML = TRUE, include_testing = FALSE )

    expect_error( pv = estimate_ATE_RIRC( Yobs, Z, sid, data=sdf, REML = TRUE, include_testing = TRUE ) )


} )




test_that("RICC functions work", {
    set.seed( 1019 )

    sdf = generate_multilevel_data( n.bar=10, J=10,
                   rho2.0W = 0.3, rho2.1W = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv = estimate_ATE_RICC( Yobs, Z, sid, data=sdf )
    pv
    expect_true( is.na( pv$deviance ) )


    # Check arb variable names works.
    sdf2 = rename( sdf, YY = Yobs, ZZ = Z, S.I.D = sid )
    pv2 = estimate_ATE_RICC( YY, ZZ, S.I.D, data=sdf2 )
    pv2

    expect_equal( pv, pv2 )

} )

