library( testthat )
library( blkvar )

context("Checking blocking functions (Nicole package parts)")


test_that("Nicole block estimators provide reasonable answers", {

    set.seed( 1019 )
    dat = make.obs.data.linear( method="big")
    head( dat )
    rs <- calc.summary.stats( dat$Yobs, dat$Z, dat$B, add.neyman = TRUE)
    rs

    f1 = fitdata( Yobs, Z, B, data=dat )

    f2 = fitdata.sumtable( rs )

    expect_equal( f1, f2 )

    expect_true( f1$var_est > 0 )

    rctyes = estimate.ATE.design.based.from.stats( sum_tab = rs, method = "finite", weight = "individual" )
    rctyes
    expect_equal( f1$tau_est, rctyes$tau.hat )
    expect_equal( f1$var_est, rctyes$SE^2 )

})



test_that("With small block estimators Nicole estimators provide answers", {

    set.seed( 1019 )
    dat = make.obs.data.linear( method="small")
    head( dat )
    rs <- calc.summary.stats( dat$Yobs, dat$Z, dat$B, add.neyman = TRUE)
    rs

    expect_warning( f1 <- fitdata( Yobs, Z, B, data=dat ) )

    f1 <- fitdata( Yobs, Z, B, data=dat, method="hybrid_p" )
    f1
    expect_equal( f1$percent_small_blocks, 50 )
    expect_true( f1$var_est > 0 )

    f2 = fitdata.sumtable( rs, method="hybrid_p" )

    expect_equal( f1, f2 )


    rctyes = estimate.ATE.design.based.from.stats( sum_tab = rs, method = "finite", weight = "individual" )
    rctyes
    expect_true( is.na( rctyes$SE ) )

    expect_equal( f1$tau_est, rctyes$tau.hat )

})




test_that("Dropping levels works", {

    set.seed( 1019 )
    dat = make.obs.data.linear( method="big")
    dat$Z[1:4] = 0

    expect_false( all( table( dat$Z, dat$B ) == 0 ) )

    expect_warning( rs <- calc.summary.stats( dat$Yobs, dat$Z, dat$B ) )

    expect_equal( nrow(rs), 3 )

})



