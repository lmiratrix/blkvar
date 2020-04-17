library( testthat )
library( blkvar )

context("Checking blocking functions (Nicole package parts)")


test_that("Nicole block estimators provide reasonable answers", {

    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( method="big")
    head( dat )
    rs <- calc_summary_stats( dat$Yobs, dat$Z, dat$B, add.neyman = TRUE)
    rs

    dat2 = dat
    head( dat )
    dat2 = dplyr::rename( dat2, BBB=B )
    rs2 <- calc_summary_stats( Yobs, Z, BBB, data=dat2, add.neyman = TRUE)
    rs2
    expect_equal( rs, rs2 )

    f1 = block_estimator( Yobs, Z, B, data=dat )

    f2 = block_estimator_tabulated( rs )

    expect_equal( f1, f2 )

    expect_true( f1$var_est > 0 )

    rctyes = estimate_ATE_design_based_from_stats( sum_tab = rs, method = "finite", weight = "individual" )
    rctyes
    expect_equal( f1$tau_est, rctyes$tau_hat )
    expect_equal( f1$var_est, rctyes$SE^2 )

})


test_that("All options of method for block_estimator_tabulated run without error", {
    set.seed( 101910 )
    dat = generate_blocked_data_obs_linear( method="big")
    head( dat )
    rs <- calc_summary_stats( dat$Yobs, dat$Z, dat$B, add.neyman = TRUE)
    rs
    f2 = block_estimator_tabulated( rs )
    f2

    methods = list("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small")
    grabs = lapply( methods, function( x ) {
        block_estimator_tabulated( rs, method=x ) } )

    expect_true( length( grabs ) == length( methods ) )
} )





test_that("calc_summary_stats functions work", {
    datr = data.frame( TX = c(0,0,1,1,0,0,0,1,1,1,1,1),
                      BK = c(1,1,1,1,2,2,2,2,2,2,2,2),
                      Y = c(1,3,10,13,1,2,3,1,2,3,4,5) )
    datr

    #dat = generate_blocked_data_obs( method="small")
    nrow( datr )

    sdat = calc_summary_stats( Y, TX, BK, datr )
    sdat
    expect_equal( sdat$Ybar0, c(2,2) )
    expect_equal( sdat$Ybar1, c(11.5,3) )
    expect_equal( sdat$n0, c(2,3) )
    expect_equal( sdat$n1, c(2,5) )
    expect_equal( sdat$n, c(4,8) )
    expect_equal( sdat$var0,
                  c( var( c(1,3) ), var( c(1,2,3) ) ) )


    sdat2 = calc_summary_stats_formula( Y~TX*BK, data=datr )
    expect_equal( sdat, sdat2 )

} )





test_that("calc_summary_stats with site works", {
    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( X=1:50, method="big" )
    #dat = generate_blocked_data_obs( method="small")
    nrow( dat )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )
    table( dat$siteNo )
    table( dat$B )

    sdat = calc_summary_stats( dat, siteID="siteNo" )
    sdat
    expect_true( "siteID" %in% names(sdat) )
    expect_true( all( sdat$siteID %in% unique( dat$siteNo ) ) )
} )







test_that("With small block estimators Nicole estimators provide answers", {

    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( method="small")
    head( dat )
    rs <- calc_summary_stats( dat$Yobs, dat$Z, dat$B, add.neyman = TRUE)
    rs

    expect_warning( f1 <- block_estimator( Yobs, Z, B, data=dat ) )

    f1 <- block_estimator( Yobs, Z, B, data=dat, method="hybrid_p" )
    f1
    expect_equal( f1$percent_small_blocks, 50 )
    expect_true( f1$var_est > 0 )

    f2 = block_estimator_tabulated( rs, method="hybrid_p" )

    expect_equal( f1, f2 )


    rctyes = estimate_ATE_design_based_from_stats( sum_tab = rs, method = "finite", weight = "individual" )
    rctyes
    expect_true( is.na( rctyes$SE ) )

    expect_equal( f1$tau_est, rctyes$tau_hat )

})




test_that("Dropping levels works", {

    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( method="big")
    dat$Z[1:4] = 0

    expect_false( all( table( dat$Z, dat$B ) == 0 ) )

    expect_warning( rs <- calc_summary_stats( dat$Yobs, dat$Z, dat$B ) )

    expect_equal( nrow(rs), 3 )

})
