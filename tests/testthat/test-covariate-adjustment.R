library( testthat )
library( blkvar )
context("Checking covariate adjustment works")

test_that("Check compare_methods with covariate adjustment", {

    set.seed( 1019 )
    dat = gen.dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )
    dat$X1 = dat$W + rnorm( nrow(dat) )
    dat$X2 = dat$Y0 + rnorm( nrow( dat ) )

    rs = compare_methods( Yobs, Z, sid, data=dat, control.formula = ~X1 + X2 )

    rs2 = compare_methods( Yobs, Z, sid, data=dat )

    expect_equal( ncol( rs ), ncol( rs2 ) )
    expect_equal( nrow( rs ), nrow( rs2 ) )

    if ( FALSE ) {
        head( rs )
        library( tidyverse )
        qplot( rs$tau, rs2$tau )
        rs2$tau.adj = rs$tau
        rs2 = mutate( rs2, delta = tau - tau.adj )
        rs2

    }
})






test_that( "Covariate adjusted Design based works through compare_methods", {
    set.seed( 1019 )
    dat = gen.dat( n.bar = 30, J = 4 )
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

    rsA =  compare_methods( Yobs, Z, sid, data=dat, include.MLM = FALSE, include.block = FALSE,
                            control.formula = ~ X1 + X2)
    rsA

    rsA$tau[ rsA$method == "DB-FP-Persons-adj" ]
    ATE
    expect_false( rsA$tau[ rsA$method == "DB-FP-Persons-adj" ] == ATE )

    expect_true( is.data.frame(rsA) )



})



