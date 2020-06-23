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




test_that("Q stat weighting methods work", {



    # Variation correlated with precision
    N = 500
    SE_hat = seq( 0.01, 0.4, length.out = N )
    ATE_hat = sort( rnorm( N, mean = 1, sd = 0.2 ) ) + rnorm( N, sd=SE_hat )
    #ATE_hat = rnorm( N, mean = 1, sd = 0.2 ) + rnorm( N, sd=SE_hat )

    rs = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE ) #, verbose = TRUE )
    rs
    rs$ATE_hat = weighted.mean( ATE_hat, w=1/SE_hat^2 )

    rs2 = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE, mean_method = "raw" ) #, verbose = TRUE )
    rs2
    rs2$ATE_hat = mean(ATE_hat)


    rs3 = analysis_Qstatistic_stat( ATE_hat, SE_hat, calc_CI = TRUE, mean_method = 1 ) #, verbose = TRUE )
    rs3
    rs3$ATE_hat = 1

    alrs = bind_rows( as.data.frame( rs ), as.data.frame( rs2 ), as.data.frame( rs3 ) )
    alrs

    expect_true( nrow( alrs ) == 3 )



})
