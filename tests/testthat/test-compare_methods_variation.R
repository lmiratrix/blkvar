


test_that( "component methods work", {
    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )


    # FIRC model (separate variances)
    FIRC <- estimate_ATE_FIRC(Yobs, Z, sid, data=dat, include_testing = TRUE)
    FIRC

    # FIRC model (with pooled residual variances)
    FIRC_pool <- estimate_ATE_FIRC(Yobs, Z, sid,  data=dat, include_testing = TRUE, pool = TRUE)

    # the random-intercept, random-coefficient (RIRC) model
    RIRC <- estimate_ATE_RIRC(Yobs, Z, sid,  data=dat, include_testing = TRUE)

    # the random-intercept, random-coefficient (RIRC) model
    RIRC_pool <- estimate_ATE_RIRC(Yobs, Z, sid,  data=dat, include_testing = TRUE,
                                   pool = TRUE)



    rsA =  analysis_Qstatistic( Yobs, Z, sid, data=dat, calc_CI = TRUE )
    rsA

    expect_true( all( !is.na( rsA ) ) )

})





test_that( "Comparing variation methods works", {
    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rsA =  compare_methods_variation( Yobs, Z, sid, data=dat )
    rsA
    expect_true( is.data.frame(rsA) )

    rs =  compare_methods_variation( Yobs, Z, sid, data=dat, include_testing = FALSE )
    rs
    expect_true( all( !is.na( rs) ) )

    expect_equal( ncol( rs ), 5 )

    rs =  compare_methods_variation( Yobs, Z, sid, data=dat, long_results = TRUE )
    rs
    expect_equal( nrow( rs ), 5 )
    expect_true( all( !is.na( rs$tau_hat ) ) )
    expect_true( all( !is.na( rs$pv ) ) )

    rs2 = compare_methods_variation( dat$Yobs, dat$Z, dat$sid, long_results = TRUE )
    expect_equal( rs, rs2 )

    head( dat )
    names( dat ) = c( "ID", "y.0", "y.1", "Tx", "outcome",  "W" )
    d2 = dat
    d2$y.0 = d2$y.1 = d2$W = NULL
    head( d2 )
    rs2 = compare_methods_variation( outcome, Tx, B=ID, data=dat, long_results = TRUE )
    expect_equal( rs, rs2 )

})







test_that( "Comparing variation with site and block works", {
    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )
    dat$siteNo = as.factor( round( as.numeric( dat$sid ) / 4 ) )

    # With site blocking
    rsA =  compare_methods_variation( Yobs, Z, sid, siteID = "siteNo", data=dat )
    rsA
    expect_true( is.data.frame(rsA) )

    # without site blocking
    rsA2 =  compare_methods_variation( Yobs, Z, sid, data=dat )
    rsA2

    expect_equal( rsA$tau_hat.RIRC, rsA2$tau_hat.RIRC )
    expect_equal( ncol( rsA ), ncol( rsA2 ) )
    expect_true( rsA$tau_hat_FIRC != rsA2$tau_hat_FIRC )
    expect_true( rsA$tau_hat_FIRC_pool != rsA2$tau_hat_FIRC_pool )

})

