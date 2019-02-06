##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking DGP functions")

test_that("data generation functions work", {
    dat = make.data.linear()
    expect_equal( nrow(dat), 16 )

    dat = make.data.linear()
    dat$blk = make.blocks( dat$X, method="pair" )
    expect_true( all( table( dat$blk ) == 2 ) )

    dat$blk = make.blocks( dat$X, method="none" )
    expect_true( length( unique( dat$blk ) ) == 1 )

    dat$blk = make.blocks( dat$X, method="small" )
    expect_true( all( table( dat$blk ) > 1 ) )


    dat = add.obs.data( dat, blockvar="blk" )
    expect_true( all( table( dat$blk, dat$Z ) > 0 ) )

    dat$blk = make.blocks( dat$X, method="big" )
    expect_true( all( table( dat$blk ) >= 4 ) )

    dat = add.obs.data( dat, blockvar="blk" )
    expect_true( all( table( dat$blk, dat$Z ) >= 2 ) )

})




test_that( "testing equispaced X", {

    dat = make.obs.data.linear( method="small", X=1:10 )
    dat
    table( dat$B )
    expect_equal( nrow(dat), 10 )
} )



test_that( "out of order X", {

    dat = make.blocks( method="small", X=c(1.1,1.2,1.3,4.1,4.2,4.3,1.4,1.5), num.blocks = 2 )
    expect_equal( dat, c( "B1", "B1", "B1", "B2", "B2", "B2", "B1", "B1" ) )

    dat = make.blocks( method="small", X=c(1.1,1.7,1.3,4.1,4.2,4.3,1.2,1.8) )
    dat
    expect_equal( dat, c( "B1", "B2", "B1", "B3","B3","B3", "B1", "B2" ) )

} )


test_that( "Pashley paper data generators work", {
    dat = make.data( c( 2,3,4) )
    dat
    expect_equal( as.numeric( table( dat$B ) ), c( 2, 3, 4 ) )

    dat = make.obs.data( c( 2,3,4, 20), p=0.001 )
    dat
    expect_equal( sum(dat$Z), 4 )

})


test_that( "Block specific generator works", {
    dt = generate.individuals.from.blocks( c( 4, 3 ), exact=TRUE )
    dt
    expect_equal( nrow( dt ), 7 )
    expect_equal( length( unique( dt$B ) ), 2 )
    expect_true( is.factor( dt$B ) )
    head( dt )
    dt = add.obs.data( dt )

    ss = calc.summary.stats.oracle( dt )
    ss
    expect_equal( ss$sd0, c( 1, 1 ) )
    expect_equal( ss$sd1, c( 1, 1 ) )
    expect_equal( ss$corr, c( 1, 1 ) )
    expect_equal( ss$tau, c( 0, 0 ) )
    expect_equal( ss$mu0, c( 0, 0 ) )

    dt = generate.individuals.from.blocks( c( 4, 8 ), c( 0, 10 ),  c( 10, 1 ), c(1, 3), c( 3, 1 ), c( 0, 1 ), TRUE )
    dt
    expect_equal( nrow( dt ), 12 )
    expect_equal( length( unique( dt$B ) ), 2 )
    expect_true( is.factor( dt$B ) )

    dt = add.obs.data( dt )
    ss = calc.summary.stats.oracle( dt )
    ss
    expect_equal( ss$tau, c( 10, 1 ) )
    expect_equal( ss$mu0, c( 0, 10 ) )

    dt = generate.individuals.from.blocks( c( 4, 8 ), c( 0, 10 ),c( 10, 1 ),  c(1, 3), c( 3, 1 ), c( 0, 1 ), FALSE )
    dt = add.obs.data( dt )
    ss = calc.summary.stats.oracle( dt )
    ss
    expect_true( ss$tau[[1]] > ss$tau[[2]] )

})


test_that( "Block factors are factors and in increasing order even if we hit 10+ blocks", {
    dat = make.data( rep( c( 5, 10, 20, 40, 60 ), each=2 ), sigma_alpha = 2, sigma_tau=1 )
    head( dat )
    expect_true( is.factor( dat$B ) )
    expect_equal( nlevels( dat$B ), 10 )
    levels( dat$B )
    expect_equal( levels( dat$B ), paste( "B", 1:10, sep="" ) )
} )



