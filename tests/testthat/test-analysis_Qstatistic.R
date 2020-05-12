test_that("Q stat test inversion works", {


    # More Demo of code working

    # Small amounts of variation, but very precise
    SE_hat = rep(0.005,30)
    ATE_hat = rnorm( 30, sd = 0.2 ) + rnorm( 30, sd=SE_hat )
    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat > 0 )
    expect_true( rs$p.value < 0.01 )

    # Even smaller amounts of variation, but very precise
    SE_hat = rep(0.01,300)
    ATE_hat = rnorm( 300, sd = 0.01 ) + rnorm( 300, sd=SE_hat )
    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat > 0 )
    expect_true( rs$p.value < 0.01 )


    # Not much variation but lots of thought to exist variation
    A = rnorm( 40 )
    rs = analysis_Qstatistic_stat( A, rep( 1.1, 40 ), calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat == 0 )
    expect_true( rs$CI_low == 0 )
    expect_true( rs$p.value > 0.05 )

    # Not much variation but lots and lots of thought to exist variation
    A = rnorm( 100 )
    rs = analysis_Qstatistic_stat( A, rep( 10, 100 ), calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat == 0 )
    expect_true( rs$CI_high == 0 )

})
