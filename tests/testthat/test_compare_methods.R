library( testthat )
library( blkvar )
context("Checking compare methods all work")

test_that("Check call options of the compare_methods", {

    set.seed( 1019 )
    dat = gen.dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rs = compare_methods( Yobs, Z, sid, data=dat, include.MLM = FALSE )

    rs2 = compare_methods( Y=dat$Yobs, Z=dat$Z, B=dat$sid, include.MLM = FALSE )

    rs3 = compare_methods( data=dat[c("Yobs","Z","sid")], include.MLM = FALSE )

    expect_equal( rs, rs2 )
    expect_equal( rs, rs3 )


})



test_that("Check for lack warnings, etc., from compare_methods", {

    set.seed( 1019 )
    dat = gen.dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )


    expect_warning( compare_methods( Yobs, Z, sid, data=dat ), regexp = NA)

})




test_that("Check for asking for different parts", {

    set.seed( 1019 )
    dat = gen.dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include.block = FALSE, include.LM = FALSE, include.RCTYes = FALSE )
    rs
    expect_equal( nrow( rs ), 2 )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include.block = FALSE, include.LM = FALSE, include.RCTYes = FALSE, include.MLM = FALSE )
    expect_equal( nrow( rs ), 0 )

})


test_that( "Comparing variation methods works", {
    set.seed( 1019 )
    dat = gen.dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rsA =  compare_methods_variation( Yobs, Z, sid, data=dat )
    rsA
    expect_true( is.data.frame(rsA) )

    rs =  compare_methods_variation( Yobs, Z, sid, data=dat, include.testing = FALSE )
    rs

    expect_equal( ncol( rs ), 4 )

    rs =  compare_methods_variation( Yobs, Z, sid, data=dat, long.results = TRUE )
    rs
    expect_equal( nrow( rs ), 5 )

    rs2 = compare_methods_variation( dat$Yobs, dat$Z, dat$sid, long.results = TRUE )
    expect_equal( rs, rs2 )

    names( dat ) = c( "ID", "W", "y.0", "y.1", "Tx", "outcome" )
    d2 = dat
    d2$y.0 = d2$y.1 = d2$W = NULL
    head( d2 )
    rs2 = compare_methods_variation( outcome, Tx, B=ID, data=dat, long.results = TRUE )
    expect_equal( rs, rs2 )

})
