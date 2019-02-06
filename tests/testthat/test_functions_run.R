library( testthat )
library( blkvar )
context("Checking analysis functions run")

test_that("Summary core function works", {

    dat = make.obs.data.linear( method="big")

    rs = calc.summary.stats( dat$Yobs, dat$Z, dat$B )
    rs
    expect_equal( nrow(rs), 4 )
    expect_equal( rs$n0, c(2,2,2,2) )

    rs2 = calc.summary.stats( data =  dat )
    rs2
    expect_equal( rs, rs2 )
})



test_that("compare functions doesn't crash when called on big blocks", {
    dat = make.obs.data.linear( method="big")
    res <- calc.summary.stats(dat$Yobs, dat$Z, dat$B)

    res <- fitdata(dat$Yobs, dat$Z, dat$B, method=c("hybrid_m") )
    res
    expect_equal( class( res ), "var_dat" )

    res <- compare_methods( dat$Yobs, dat$Z, dat$B, include.MLM = FALSE )
    res
    expect_true( sd( res$tau ) <= 0.0001 )
})



test_that("compare functions doesn't crash when called on small blocks", {
    dat = make.obs.data.linear( method="small")
    dat

    res <- calc.summary.stats(dat$Yobs, dat$Z, dat$B)

    expect_warning( res <- fitdata(dat$Yobs, dat$Z, dat$B, method=c("hybrid_m") ) )

    res <- compare_methods( dat$Yobs, dat$Z, dat$B, include.MLM = FALSE )
    res
})



test_that( "further tests of output, etc", {
    B.a<-c(rep(1,3), rep(2,6), rep(4,7),rep(5,2), rep(6,2))
    Z.a<-rep(c(1,0),10)
    Y.a<-Z.a*10+B.a+rnorm(20,0,1)

    A  = calc.summary.stats(Y=Y.a, Z=Z.a, B=B.a)
    A

    data.a<-cbind(Y.a, Z.a, B.a)
    data.a
    head( data.a )
    class( data.a )
    B = calc.summary.stats( Y.a, Z.a, B.a, data=as.data.frame( data.a ) )

    expect_equal( A, B )

    nsmall = sum( A$n[ A$n0 == 1 | A$n1 == 1 ] )
    nsmall / sum( A$n )

    expect_warning( method.hybrid.m<-fitdata(Y.a, Z.a, B.a, method="hybrid_m") )

    method.hybrid.p<-fitdata(Y.a, Z.a, B.a, method="hybrid_p")
    method.hybrid.p
    expect_equal( class( method.hybrid.p ), "var_dat" )

    # We can pull out whatever we want
    expect_true( !is.na( method.hybrid.p$tau_est ) )

    expect_true( !is.na( method.hybrid.p$se_est ) )

    expect_equal( method.hybrid.p$percent_small_blocks, 35 )
    expect_equal( dim( method.hybrid.p$block_sizes ), c( 5,3) )

    comp = compare_methods(Y.a, Z.a, B.a, include.MLM = FALSE)
    comp
    expect_equal( ncol( comp ), 3 )
#    expect_equal( nrow( comp ), 5 )
} )


