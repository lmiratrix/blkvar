


test_that("Nice formatting of result object", {
    set.seed( 1010101 )

    dt = generate_multilevel_data( n.bar = 100, J = 40,
                                   ICC = 0.5, rho2.1W = 0,
                                   tau.11.star = 0.45, zero.corr = FALSE )

    head( dt )
    rs = analysis_Qstatistic( Yobs, Z, sid, data=dt )
    rs

    expect_true( class(rs) == "multisiteresult" )

    rs_FIRC = estimate_ATE_FIRC( Yobs, Z, sid, data=dt )
    rs_FIRC
    expect_true( class(rs_FIRC) == "multisiteresult" )

    rs_RIRC = estimate_ATE_RIRC( Yobs, Z, sid, data=dt,
                                 pool = TRUE,
                                 control_formula = ~ W )
    rs_RIRC
    expect_true( class(rs_RIRC) == "multisiteresult" )

    rs_RICC = estimate_ATE_RICC( Yobs, Z, sid, data=dt )
    rs_RICC

    expect_true( class(rs_RICC) == "multisiteresult" )

})



