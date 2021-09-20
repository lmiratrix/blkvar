


test_that("Q stat test inversion works", {

    # Small amounts of variation, but very precise
    set.seed( 1010101 )
    SE_hat = rep(0.005,30)
    ATE_hat = rnorm( 30, sd = 0.2 ) + rnorm( 30, sd = SE_hat )
    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat > 0 )
    expect_true( rs$p.value < 0.01 )
    expect_true( rs$CI_low <= 0.2 )
    expect_true( rs$CI_high >= 0.2 )

    # Even smaller amounts of variation, but very precise
    SE_hat = rep(0.01, 300)
    ATE_hat = rnorm( 300, sd = 0.01 ) + rnorm( 300, sd=SE_hat )
    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat > 0 )
    expect_true( abs( rs$tau_hat - 0.01 ) < 0.01 )
    expect_true( rs$p.value < 0.01 )
    expect_true( rs$CI_high >= 0.01 )
    expect_true( rs$CI_low <= 0.01 )
    expect_true( rs$CI_low > 0.01^2 )

    # Not much variation but lots of thought to exist variation (so we would be
    # up against a boundary of no variation)
    A = rnorm( 40 )
    rs = analysis_Qstatistic_stat( A, rep( 1.1, 40 ), calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat == 0 )
    expect_true( rs$CI_low == 0 )
    expect_true( rs$p.value > 0.05 )

    # Not much variation but lots and lots of thought to exist variation (so we
    # are sure there is no possible cross site variation).
    A = rnorm( 100 )
    rs = analysis_Qstatistic_stat( A, rep( 10, 100 ), calc_CI = TRUE )
    rs
    expect_true( rs$tau_hat == 0 )
    expect_true( rs$CI_high == 0 )

})





test_that("Q stat with massive cross site variation works", {

    set.seed( 1010101 )
    SE_hat = rnorm( 100, 10, 2 ) / 5
    ATE_hat = rnorm( 100, sd = 10 ) + rnorm( 100, sd = SE_hat )

    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE ) #, verbose = TRUE )
    rs

    expect_true( abs(rs$tau_hat - 10) <= 0.5 )
    expect_true( rs$CI_high >= 10 )
    expect_true( rs$CI_low <= 10 )


    rs2 = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE,
                                   calc_optim_pvalue = TRUE ) #, verbose = TRUE )
    rs2

    expect_true( abs(rs2$tau_hat - 10) <= 0.5 )
    expect_true( rs2$CI_high >= 10 )
    expect_true( rs2$CI_low <= 10 )

    expect_true( abs( rs$tau_hat - rs2$tau_hat ) <= 0.25 )

})




test_that("Q stat weighting methods work", {


    # Variation correlated with precision
    N = 500
    SE_hat = seq( 0.01, 0.4, length.out = N )
    ATE_hat = sort( rnorm( N, mean = 1, sd = 0.2 ) ) + rnorm( N, sd=SE_hat )
    #ATE_hat = rnorm( N, mean = 1, sd = 0.2 ) + rnorm( N, sd=SE_hat )

    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE ) #, verbose = TRUE )
    rs
    expect_true( rs$ATE_hat == weighted.mean( ATE_hat, w=1/SE_hat^2 ) )

    rs2 = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE, mean_method = "raw" ) #, verbose = TRUE )
    rs2
    expect_true( rs2$ATE_hat == mean(ATE_hat) )


    rs3 = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE, mean_method = 1 ) #, verbose = TRUE )
    rs3
    expect_true( rs3$ATE_hat == 1 )

    alrs = bind_rows( as.data.frame( rs ), as.data.frame( rs2 ), as.data.frame( rs3 ) )
    alrs

    expect_true( nrow( alrs ) == 3 )


})
