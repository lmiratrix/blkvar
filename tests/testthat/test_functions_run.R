library( testthat )
library( blkvar )
context("Checking analysis functions run")

test_that("Summary core function works", {

    dat = make.obs.data.linear( method="big")

    rs = block.data( dat$Yobs, dat$Z, dat$blk )
    rs
    expect_equal( nrow(rs), 4 )
    expect_equal( rs$Num_Ctrl, c(2,2,2,2) )


})



test_that("compare functions doesn't crash when called on big blocks", {
    dat = make.obs.data.linear( method="big")
    res <- block.data(dat$Yobs, dat$Z, dat$blk)

    res <- fitdata(dat$Yobs, dat$Z, dat$blk, method=c("hybrid_m") )
    res
    expect_equal( class( res ), "var_dat" )

    res <- compare_methods( dat$Yobs, dat$Z, dat$blk, include.MLM = FALSE )
    res
    expect_true( sd( res$tau ) <= 0.0001 )
})



test_that("compare functions doesn't crash when called on small blocks", {
    dat = make.obs.data.linear( method="small")
    dat

    res <- block.data(dat$Yobs, dat$Z, dat$blk)

    expect_warning( res <- fitdata(dat$Yobs, dat$Z, dat$blk, method=c("hybrid_m") ) )

    res <- compare_methods( dat$Yobs, dat$Z, dat$blk, include.MLM = FALSE )
    res
})



test_that( "further tests of output, etc", {
    B.a<-c(rep(1,3), rep(2,6), rep(4,7),rep(5,2), rep(6,2))
    Z.a<-rep(c(1,0),10)
    Y.a<-Z.a*10+B.a+rnorm(20,0,1)

    A  = block.data(Y=Y.a, Z=Z.a, B=B.a)

    data.a<-cbind(Y.a, Z.a, B.a)
    B = block.data(data=data.a)

    expect_equal(A, B )

    expect_warning( method.hybrid.m<-fitdata(Y.a, Z.a, B.a, method="hybrid_m") )

    method.hybrid.p<-fitdata(Y.a, Z.a, B.a, method="hybrid_p")
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


