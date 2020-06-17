library( testthat )
library( blkvar )
context("Checking compare methods all work")

test_that("Check different call options of the compare_methods", {

    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rs = compare_methods( Yobs, Z, sid, data=dat, include_MLM = FALSE )
    rs

    expect_true( all( rs$SE > 0 ) )

    rs2 = compare_methods( Y=dat$Yobs, Z=dat$Z, B=dat$sid, include_MLM = FALSE )

    rs3 = compare_methods( data=dat[c("Yobs","Z","sid")], include_MLM = FALSE )

    expect_equal( rs, rs2 )
    expect_equal( rs, rs3 )


})



test_that("Check for lack warnings, etc., from compare_methods", {

    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )


    expect_warning( compare_methods( Yobs, Z, sid, data=dat ), regexp = NA)

})



test_that("Check method_characteristics works", {

    set.seed( 101974 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, include_method_characteristics = TRUE )

    expect_true( all( c( "weight","population","biased") %in% names(fulltab) ) )
    expect_true( all( !is.na( fulltab$weight ) ) )

    ptab = compare_methods( Yobs, Z, sid, data=dat, include_method_characteristics = FALSE )

    mc = method_characteristics()
    expect_true( all( fulltab$method %in% mc$method ) )

    expect_equal( nrow( fulltab ), nrow( ptab ) )
})




test_that("weighting lm flags gets passed", {

    set.seed( 102030 )
    dat = generate_multilevel_data( n.bar = 20, J = 20 )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, weight_LM_scale_weights = FALSE,
                               include_MLM = FALSE, include_block = FALSE, include_DBBlended = FALSE, include_DB = FALSE )
    expect_true( "FE-IPTW(n)" %in% fulltab$method )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, weight_LM_scale_weights = TRUE,
                               include_MLM = FALSE, include_block = FALSE, include_DBBlended = FALSE, include_DB = FALSE )
    fulltab
    expect_false( "FE-IPTW(n)" %in% fulltab$method )

})




test_that("Check method works on tibbles", {

    set.seed( 101974 )
    dat = tibble::as_tibble( generate_multilevel_data( n.bar = 30, J = 30 ) )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, include_method_characteristics = TRUE )
    fulltab
    expect_true( all( c( "weight","population","biased") %in% names(fulltab) ) )
    expect_true( all( !is.na( fulltab$weight ) ) )

    fulltab = compare_methods( Yobs, Z, sid, data=dat, include_method_characteristics = FALSE )
    mc = method_characteristics()
    expect_true( all( fulltab$method %in% mc$method ) )


})


test_that("Check for asking for different parts", {

    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 30 )
    nrow( dat )
    head( dat )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include_block = FALSE, include_LM = FALSE, include_DB = FALSE )
    rs
    # Three MLM methods
    expect_equal( nrow( rs ), 3 )
    expect_true( names(rs)[[1]] == "method" )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include_block = FALSE, include_LM = TRUE, include_DB = FALSE,
                           include_MLM = FALSE)
    rs
    expect_true( names(rs)[[1]] == "method" )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include_block = FALSE,
                           include_MLM = FALSE)
    rs
    expect_true( names(rs)[[1]] == "method" )

    rs =  compare_methods( Yobs, Z, sid, data=dat, include_block = FALSE,
                           include_LM = FALSE, include_DB = FALSE, include_MLM = FALSE )
    expect_equal( nrow( rs ), 0 )

})






test_that( "Design based works through compare_methods", {
    set.seed( 1019 )
    dat = generate_multilevel_data( n.bar = 30, J = 4 )
    nrow( dat )
    head( dat )

    blocks = dat %>% group_by( sid ) %>%
        summarise( ATE.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   n=n() )
    blocks

    ATE = weighted.mean( blocks$ATE.hat, blocks$n )
    ATE

    rsA =  compare_methods( Yobs, Z, sid, data=dat, include_MLM = FALSE,
                            include_block = FALSE )
    rsA

    rsA$ATE_hat[ rsA$method == "DB-FP-Persons" ]
    ATE
    expect_equal( rsA$ATE_hat[ rsA$method == "DB-FP-Persons" ], ATE )

    sdat = calc_summary_stats( Yobs, Z, sid, data=dat )
    sdat
    sdat = mutate( sdat, ATE.hat = Ybar1-Ybar0 )
    sdat
    a = estimate_ATE_design_based_from_stats( sdat, weight="individual",
                                              method="finite" )
    a
    expect_equal( a$ATE_hat, ATE )

    expect_true( is.data.frame(rsA) )



})



