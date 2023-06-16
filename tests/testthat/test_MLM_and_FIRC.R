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


    pv.aov = estimate_ATE_FIRC( Yobs, Z, sid, data=sdf, include_testing = TRUE, anova=TRUE,
                                keep_EB_estimates = FALSE )
    expect_equal( pv2$p_variation, pv.aov$p_variation )


    # Check arb variable names works.
    sdf2 = rename( sdf, YY = Yobs, ZZ = Z, S.I.D = sid )

    pv2.B = estimate_ATE_FIRC(YY, ZZ, S.I.D, data=sdf2, include_testing = TRUE )

    expect_equal( pv2.B, pv2 )


    eb = get_EB_estimates( pv2.B )
    expect_true( nrow(eb) == 10 )

    expect_true( is.null( get_EB_estimates(pv.aov) ) )

    REML = estimate_ATE_FIRC(YY, ZZ, S.I.D, data=sdf2, REML = TRUE )
    expect_true( REML$tau_hat > pv2.B$tau_hat )
    expect_true( sd( get_EB_estimates(REML)$beta_hat ) > sd( eb$beta_hat ) )
} )





test_that("RIRC functions work", {
    set.seed( 1019 )

    sdf = generate_multilevel_data( n.bar=10, J=11,
                   rho2.0W = 0.3, rho2.1W = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )
    head( sdf )
    pv = estimate_ATE_RIRC( Yobs, Z, sid, data=sdf, keep_EB_estimates = FALSE )
    pv
    pv = estimate_ATE_RIRC( Yobs, Z, sid, data=sdf, REML = TRUE, include_testing = FALSE )

    expect_error( pv = estimate_ATE_RIRC( Yobs, Z, sid, data=sdf, REML = TRUE, include_testing = TRUE ) )

    eb = get_EB_estimates( pv )
    expect_true( nrow(eb) == 11 )

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

