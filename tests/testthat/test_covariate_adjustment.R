library( testthat )
library( blkvar )
context("Checking covariate adjustment works")


test_that( "make_canonical_data works", {
    set.seed( 221019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )
    dat$X1 = dat$W + rnorm( nrow(dat) )
    dat$X2 = dat$Y0 + rnorm( nrow( dat ) )

    rs = compare_methods( Yobs, Z, sid, data=dat, control_formula = ~ X1 + X2 )

    dd = blkvar:::make_canonical_data( formula = Yobs ~ Z * sid,
                                       control_formula = ~ X1 + X2,
                                       data = dat )
    head( dd )
    expect_true( ncol(dd) == 6 )
    expect_true( all( dd$X2 == dat$X2 ) )
    expect_true( all( dd$siteID == dd$B ) )
    dd = blkvar:::make_canonical_data( formula = Yobs ~ Z * sid,
                                       data = dat )
    head( dd )
    expect_true( ncol(dd) == 4 )
})



test_that("Check compare_methods with covariate adjustment", {

    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )
    dat$X1 = dat$W + rnorm( nrow(dat) )
    dat$X2 = dat$Y0 + rnorm( nrow( dat ) )

    rs = compare_methods( Yobs, Z, sid, data=dat, control_formula = ~ X1 + X2 )

    rs2 = compare_methods( Yobs, Z, sid, data=dat )

    expect_equal( ncol( rs ), ncol( rs2 ) )
    expect_equal( nrow( rs ), nrow( rs2 ) )

    if ( FALSE ) {
        head( rs )
        qplot( rs$ATE_hat, rs2$ATE_hat )
        rs2$ATE_hat.adj = rs$ATE_hat
        rs2 = mutate( rs2, delta = ATE_hat - ATE_hat.adj )
        rs2

    }
})






test_that( "Covariate adjusted Design based works through compare_methods", {
    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 4 )
    nrow( dat )
    head( dat )
    dat$X1 = dat$W + rnorm( nrow(dat) )
    dat$X2 = dat$Y0 + rnorm( nrow( dat ) )

    blocks = dat %>% group_by( sid ) %>%
        summarise( ATE.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   n=n() )
    blocks

    ATE = weighted.mean( blocks$ATE.hat, blocks$n )
    ATE

    rsA =  compare_methods( Yobs, Z, sid, data=dat, include_MLM = FALSE, include_block = FALSE,
                            control_formula = ~ X1 + X2)
    rsA

    rsA$ATE_hat[ rsA$method == "DB-FP-Persons-adj" ]
    ATE
    expect_false( rsA$ATE_hat[ rsA$method == "DB-FP-Persons-adj" ] == ATE )

    expect_true( is.data.frame(rsA) )



})






test_that( "Covariate adjustment with categorical covariates works", {
    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 4 )
    nrow( dat )
    head( dat )
    dat$X1 = dat$W + rnorm( nrow(dat) )
    dat$X2 = cut( dat$Y0 + rnorm( nrow( dat ) ), 3 )
    table( dat$X2 )
    class(dat$X2 )


    rsAdb =  estimate_ATE_design_based_adjusted( Yobs ~ Z*sid, data=dat,
                            control_formula = ~ X1 + X2)

    rsA =  compare_methods( Yobs, Z, sid, data=dat, include_MLM = FALSE, include_block = FALSE,
                            control_formula = ~ X1 + X2)
    rsA
    expect_true( rsAdb$ATE_hat == rsA$ATE_hat[ rsA$method == "DB-FP-Persons-adj" ] )


    expect_true( is.data.frame(rsA) )



})


