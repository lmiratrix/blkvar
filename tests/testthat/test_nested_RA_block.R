

library( testthat )
library( blkvar )
context("Checking nested randomization blocks")


# Make a dataset with nested blocks and full balance of tx and co assignment so
# we can "estimate" the true parameters.
make.balanced.dataset = function(  ) {
    # Make sure site weighting is correct
    # (Generate dataset with known tau to verify.)
    a = make.data( c( 4, 4, 4, 3, 3, 6 ), tau = c( 10, 20, 30, 40, 50, 60 ), exact=TRUE )
    a$sssite = c( 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 )
    a$Z = 1
    a$Yobs = a$Y1
    b = a
    b$Z = 0
    b$Yobs = b$Y0
    a = bind_rows( a, b )
    a = bind_rows( a, filter( b, B == "B6" ) )

    a  # perfectly balanced dataset.  Three sites, one with two blocks, one with three blocks.
    # one of those blocks has p!=0.5 treated
    table( a$sssite, a$B )
    table( a$B, a$Z )

    a

}

make.big.balanced.dataset = function(rps=5) {
    rps = plyr::ldply( 1:rps, function(l) {
        df = make.balanced.dataset()
        df$sssite = paste( l, df$sssite, sep="-" )
        df$B = paste( l, df$B, sep="-" )
        df
        })
    rps
}

get.params = function( a ) {
    ss = a %>% group_by( B ) %>% summarise( tau = mean( Y1 ) - mean( Y0 ),
                                                  ybar1 = mean( Y1 ),
                                                  ybar0 = mean( Y0 ) )
    ss

    ss2 = a %>% group_by( sssite ) %>% summarise( tau = mean( Y1 ) - mean( Y0 ),
                                                  ybar1 = mean( Y1 ),
                                                  ybar0 = mean( Y0 ) )
    ss2
    true.site.tau = mean( ss2$tau )
    true.site.tau

    true.indiv.tau = mean( a$Y1 - a$Y0 )

    true.block.tau = mean( ss$tau )

    list( true.site.tau = true.site.tau,
          true.block.tau = true.block.tau,
          true.indiv.tau = true.indiv.tau,
          n.site = nrow(ss2) )
}

test_that("DB estimators work with nested randomization blocks", {

    set.seed( 1019 )
    dat = make.obs.data.linear( X=1:50, method="big" )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )

    # Our dataset with blocks in sites
    sdat = calc.summary.stats( dat, siteID="siteNo" )
    sdat

    # Finite population, person weighted
    a = estimate.ATE.design.based.from.stats( sdat, siteID="siteID", weight="individual", method="finite" )

    a2 = estimate.ATE.design.based.from.stats( sdat, weight="individual", method="finite" )
    a2
    expect_equal( a, a2 )

    # Finite population, site weighted
    # Should be different due to weighting the RA blocks differently.
    b = estimate.ATE.design.based.from.stats( sdat, siteID="siteID", weight="site", method="finite" )
    b2 = estimate.ATE.design.based.from.stats( sdat, weight="site", method="finite" )
    b2
    b
    expect_true( b$tau.hat != b2$tau.hat )
    expect_true( b$SE != b2$SE )


    # Make sure site weighting is correct
    # (Generate dataset with known tau to verify.)
    a = make.balanced.dataset()

    ss = a %>% group_by( B ) %>% summarise( tau = mean( Y1 ) - mean( Y0 ),
                                       ybar1 = mean( Y1 ),
                                       ybar0 = mean( Y0 ) )
    ss
    aa = calc.summary.stats( a, siteID="sssite" )
    aa
    expect_equal( aa$Ybar1, ss$ybar1 )

    site.tau.hat = estimate.ATE.design.based.from.stats( aa, weight="site", method="finite"  )$tau.hat
    site.tau.hat
    expect_equal( 35, site.tau.hat )

    ss2 = a %>% group_by( sssite ) %>% summarise( tau = mean( Y1 ) - mean( Y0 ),
                                            ybar1 = mean( Y1 ),
                                            ybar0 = mean( Y0 ) )
    ss2
    true.site.tau = mean( ss2$tau )

    corrtauhat = estimate.ATE.design.based.from.stats( aa, siteID="siteID", weight="site", method="finite"  )$tau.hat
    corrtauhat
    expect_equal( true.site.tau,
                  corrtauhat )



    # Superpopulation, individual
    est1 = estimate.ATE.design.based.from.stats( aa, siteID="siteID", weight="individual", method="superpop"  )
    est2 = estimate.ATE.design.based.from.stats( aa, weight="individual", method="superpop"  )
    est1
    est2
    expect_true( est1$tau.hat == est2$tau.hat )
    expect_true( est1$SE != est2$SE )

    # Superpopulation, site
    est1 = estimate.ATE.design.based.from.stats( aa, siteID="siteID", weight="site", method="superpop"  )
    est2 = estimate.ATE.design.based.from.stats( aa, weight="site", method="superpop"  )
    est1
    est2
    expect_true( est1$tau.hat == true.site.tau )
    expect_true( est2$tau.hat == site.tau.hat )
    expect_true( est1$SE != est2$SE )

    # If we pass siteID but don't actually have nesting
    est3 = estimate.ATE.design.based.from.stats( aa, siteID="B", weight="site", method="superpop"  )
    expect_true( est3$tau.hat == site.tau.hat )
    expect_true( est3$SE == est2$SE )

})





test_that("linear regression works with nested randomization blocks", {

    a = make.balanced.dataset()
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = blkvar:::interacted.linear.estimators( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    expect_equal( est1$tau[[1]], params$true.site.tau )
    expect_equal( est1$tau[[2]], params$true.indiv.tau )

    est1b = blkvar:::interacted.linear.estimators( outcome, Tx, BB, data=a )
    est1b
    expect_equal( est1b$tau[[1]], params$true.block.tau )
    expect_equal( est1b$tau[[2]], params$true.indiv.tau )

    expect_true( est1b$SE[[1]] != est1$SE[[1]] )
    expect_true( est1b$SE[[2]] == est1$SE[[2]] )


    est2 = blkvar:::fixed.effect.estimators( outcome, Tx, BB, siteID = "sssite", data=a )
    est2
    est2b = blkvar:::fixed.effect.estimators( outcome, Tx, BB, data=a )
    est2b
    expect_true( est2$SE[[3]] > est2b$SE[[3]])
    expect_equal( est2$SE[1:2], est2b$SE[1:2])


})






test_that("weighted linear regression works with nested randomization blocks", {

    a = make.balanced.dataset()
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = blkvar:::weighted.linear.estimators( outcome ~ Tx * BB, siteID = "sssite", data=a, scaled.weights = TRUE )
    est1

    est1b = blkvar:::weighted.linear.estimators( outcome ~ Tx * BB, data=a, scaled.weights = FALSE )
    est1b

    estDB = compare_methods( outcome, Tx, BB, data=a, siteID = "sssite", include.MLM = FALSE, include.LM = TRUE, include.block = FALSE )
    estDB

    expect_equal( est1$tau[[2]], estDB$tau[[2]] )
    expect_equal( est1$tau[[2]], estDB$tau[estDB$method=="DB-FP-Sites"] )
    expect_equal( est1$tau[[2]], estDB$tau[estDB$method=="FE-IPTW-Sites"] )

    expect_equal( est1$tau[[1]], params$true.indiv.tau )
    expect_equal( est1$tau[[2]], params$true.site.tau )

    expect_equal( est1b$tau[[1]], params$true.indiv.tau )
    expect_equal( est1b$tau[[2]], params$true.block.tau )


    # Another test--match design based?
    set.seed( 1019 )
    dat = make.obs.data.linear( X=1:50, method="big" )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )
    head( dat )

    # Our dataset with blocks in sites
    sdat = calc.summary.stats( dat, siteID="siteNo" )
    sdat

    # Finite population, person weighted
    a = estimate.ATE.design.based.from.stats( sdat, siteID="siteID", weight="site", method="finite" )
    a

    a2 = blkvar:::weighted.linear.estimators( Yobs ~ Z*B, siteID = "siteNo", data=dat )
    a2

    expect_equal( a$tau.hat, a2$tau[[2]] )

})


test_that("multilevel regression works with nested randomization blocks", {

    a = make.big.balanced.dataset( rps = 5 )
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = blkvar:::compare.MLM.methods( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    est1b = blkvar:::compare.MLM.methods( outcome, Tx, BB, data=a )
    est1b

    expect_true( all( est1[c(1,3),] == est1b[c(1,3),] ) )

    expect_true( est1$tau[[2]] != est1b$tau[[2]] )
    expect_true( est1$SE[[2]] > est1b$SE[[2]] )

})


test_that("all linear regression works with nested randomization blocks", {

    a = make.big.balanced.dataset( 5 )
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = linear.model.estimators( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    est1b = linear.model.estimators( outcome, Tx, BB, data=a )
    est1b

    res = merge( est1, est1b, by="method", suffixes=c(".site",".block" ) )

    res$change = res$tau.site - res$tau.block
    res$changeSE = res$SE.site - res$SE.block

    res

    ff = filter( res, change != 0 | changeSE != 0 )
    ff

    expect_true( nrow( ff ) == 4 )

})




test_that("compare_methods works with nested randomization blocks", {

    a = make.big.balanced.dataset( 5 )
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )
    head( a )

    est1 = compare_methods( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    est1b = compare_methods( outcome, Tx, BB, data=a )
    est1b

    res = merge( est1, est1b, by="method", suffixes=c(".site",".block" ) )

    res$change = res$tau.site - res$tau.block
    res$changeSE = res$SE.site - res$SE.block

    res

    ff = filter( res, change != 0 | changeSE != 0 )
    ff

    expect_true( nrow( ff ) == 8 )

})


test_that("compare_methods with adjustment works with nested randomization blocks", {
    a = gen.dat( n.bar=10, J=10 )
    head( a )
    a$X1 = a$Y0 + rnorm( nrow(a), sd=1 )
    a$X2 = a$Y1 + rnorm( nrow(a), sd=1 )

    # cluster to make sites
    a$sssite = 10 + round( as.numeric( a$sid )  / 3 )
    a$Yobs = a$Yobs + a$sssite * 0.3
    table( a$sssite )

    a = rename( a, outcome = Yobs, Tx = Z, BB = sid )
    head( a )

    est0b = compare_methods( outcome, Tx, BB, data=a )
    est0b

    est0 = compare_methods( outcome, Tx, BB, siteID = "sssite", data=a )
    est0


    est1b = compare_methods( outcome, Tx, BB, data=a,
                             control.formula = ~ X1 + X2 )
    est1b


    est1 = compare_methods( outcome, Tx, BB, data=a,
                            siteID = "sssite",
                             control.formula = ~ X1 + X2 )
    est1


    res = merge( est1, est1b, by="method", suffixes=c(".site",".block" ) )

    res$change = res$tau.site - res$tau.block
    res$changeSE = res$SE.site - res$SE.block

    options( digits= 3 )
    res

    ff = filter( res, change == 0 & changeSE == 0 )
    ff

    ff = filter( res, change != 0 | changeSE != 0 )
    ff

    expect_true( nrow( ff ) == 8 )

})





