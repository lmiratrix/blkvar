library( testthat )
library( blkvar )
context("Checking compare methods all work")

test_that("Check call options of the compare_methods", {

    set.seed( 1019 )
    dat = gen_dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rs = compare_methods( Yobs, Z, sid, data=dat, include.MLM = FALSE )
    rs
    expect_true( all( rs$SE > 0 ) )

    rs2 = compare_methods( Y=dat$Yobs, Z=dat$Z, B=dat$sid, include.MLM = FALSE )

    rs3 = compare_methods( data=dat[c("Yobs","Z","sid")], include.MLM = FALSE )

    expect_equal( rs, rs2 )
    expect_equal( rs, rs3 )


})



test_that("Check for lack warnings, etc., from compare_methods", {

    set.seed( 1019 )
    dat = gen_dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )


    expect_warning( compare_methods( Yobs, Z, sid, data=dat ), regexp = NA)

})



test_that("Check method_characteristics works", {

    set.seed( 101974 )
    dat = gen_dat( n.bar = 30, J = 30 )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, include.method_characteristics = TRUE )

    expect_true( all( c( "weight","population","biased") %in% names(fulltab) ) )
    expect_true( all( !is.na( fulltab$weight ) ) )

    ptab = compare_methods( Yobs, Z, sid, data=dat, include.method_characteristics = FALSE )

    mc = method_characteristics()
    expect_true( all( fulltab$method %in% mc$method ) )

    expect_equal( nrow( fulltab ), nrow( ptab ) )
})




test_that("weighting lm flags gets passed", {

    set.seed( 102030 )
    dat = gen_dat( n.bar = 20, J = 20 )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, weight.LM.scale.weights = FALSE,
                               include.MLM = FALSE, include.block = FALSE, include.DBBlended = FALSE, include.DB = FALSE )
    expect_true( "FE-IPTW(n)" %in% fulltab$method )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, weight.LM.scale.weights = TRUE,
                               include.MLM = FALSE, include.block = FALSE, include.DBBlended = FALSE, include.DB = FALSE )
    fulltab
    expect_false( "FE-IPTW(n)" %in% fulltab$method )

})




test_that("Check method works on tibbles", {

    set.seed( 101974 )
    dat = tibble::as_tibble( gen_dat( n.bar = 30, J = 30 ) )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, include.method_characteristics = TRUE )
    fulltab
    expect_true( all( c( "weight","population","biased") %in% names(fulltab) ) )
    expect_true( all( !is.na( fulltab$weight ) ) )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, include.method_characteristics = FALSE )
    mc = method_characteristics()
    expect_true( all( fulltab$method %in% mc$method ) )


})


test_that("Check for asking for different parts", {

    set.seed( 1019 )
    dat = gen_dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include.block = FALSE, include.LM = FALSE, include.DB = FALSE )
    rs
    # Three MLM methods
    expect_equal( nrow( rs ), 3 )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include.block = FALSE, include.LM = FALSE, include.DB = FALSE, include.MLM = FALSE )
    expect_equal( nrow( rs ), 0 )

})


test_that( "Comparing variation methods works", {
    set.seed( 1019 )
    dat = gen_dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rsA =  compare_methods_variation( Yobs, Z, sid, data=dat )
    rsA
    expect_true( is.data.frame(rsA) )

    rs =  compare_methods_variation( Yobs, Z, sid, data=dat, include.testing = FALSE )
    rs

    expect_equal( ncol( rs ), 5 )

    rs =  compare_methods_variation( Yobs, Z, sid, data=dat, long.results = TRUE )
    rs
    expect_equal( nrow( rs ), 5 )

    rs2 = compare_methods_variation( dat$Yobs, dat$Z, dat$sid, long.results = TRUE )
    expect_equal( rs, rs2 )

    head( dat )
    names( dat ) = c( "ID", "y.0", "y.1", "Tx", "outcome",  "W" )
    d2 = dat
    d2$y.0 = d2$y.1 = d2$W = NULL
    head( d2 )
    rs2 = compare_methods_variation( outcome, Tx, B=ID, data=dat, long.results = TRUE )
    expect_equal( rs, rs2 )

})







test_that( "Comparing variation with site and block works", {
    set.seed( 1019 )
    dat = gen_dat( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )
    dat$siteNo = as.factor( round( as.numeric( dat$sid ) / 4 ) )

    rsA =  compare_methods_variation( Yobs, Z, sid, siteID = "siteNo", data=dat )
    rsA
    expect_true( is.data.frame(rsA) )

    rsA2 =  compare_methods_variation( Yobs, Z, sid, data=dat )
    rsA2

    expect_equal( rsA$tau.hat.RIRC, rsA2$tau.hat.RIRC )
    expect_equal( ncol( rsA ), ncol( rsA2 ) )
    expect_true( rsA$tau.hat.FIRC != rsA2$tau.hat.FIRC )
    expect_true( rsA$tau.hat.FIRC.pool != rsA2$tau.hat.FIRC.pool )

})






test_that( "Design based works through compare_methods", {
    set.seed( 1019 )
    dat = gen_dat( n.bar = 30, J = 4 )
    nrow( dat )
    head( dat )

    blocks = dat %>% group_by( sid ) %>%
        summarise( ATE.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   n=n() )
    blocks

    ATE = weighted.mean( blocks$ATE.hat, blocks$n )
    ATE

    rsA =  compare_methods( Yobs, Z, sid, data=dat, include.MLM = FALSE,
                            include.block = FALSE )
    rsA

    rsA$tau[ rsA$method == "DB-FP-Persons" ]
    ATE
    expect_equal( rsA$tau[ rsA$method == "DB-FP-Persons" ], ATE )

    sdat = calc_summary_stats( Yobs, Z, sid, data=dat )
    sdat
    sdat = mutate( sdat, ATE.hat = Ybar1-Ybar0 )
    sdat
    a = estimate_ATE_design_based_from_stats( sdat, weight="individual",
                                              method="finite" )
    a
    expect_equal( a$tau.hat, ATE )

    expect_true( is.data.frame(rsA) )



})



