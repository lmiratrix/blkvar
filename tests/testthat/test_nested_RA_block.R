library( testthat )
library( blkvar )
context("Checking nested randomization blocks")




# Make a dataset with nested blocks and full balance of tx and co assignment so
# we can "estimate" the true parameters.
make.balanced.dataset = function(  ) {
    # Make sure site weighting is correct
    # (Generate dataset with known ATE to verify.)
    a = generate_blocked_data( c( 4, 4, 4, 3, 3, 6 ),
                               beta = c( 10, 20, 30, 40, 50, 60 ), exact=TRUE )
    a$sssite = c( 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 )
    a$Z = 1
    a$Yobs = a$Y1
    b = a
    b$Z = 0
    b$Yobs = b$Y0
    a = rbind( a, b )


    a = rbind( a,  b[which(b[, 1] == "B6"), ])

    # a  # perfectly balanced dataset.  Three sites, one with two blocks, one with three blocks.
    # one of those blocks has p!=0.5 treated
    table( a$sssite, a$B )
    table( a$B, a$Z )

    a

}

make.big.balanced.dataset = function(rps=5) {
    # rps = plyr::ldply( 1:rps, function(l) {
        # df = make.balanced.dataset()
        # df$sssite = paste( l, df$sssite, sep="-" )
        # df$B = paste( l, df$B, sep="-" )
        # df
        # })
    # rps
     full <- c()
   for (l in 1:rps) {
     df = make.balanced.dataset()
     df$sssite = paste( l, df$sssite, sep="-" )
     df$B = paste( l, df$B, sep="-" )
     full  <- rbind(full, df)
  }
  return(full)

}


get.params = function( a ) {
    ss = a %>% group_by( B ) %>% summarise( ATE = mean( Y1 ) - mean( Y0 ),
                                                  ybar1 = mean( Y1 ),
                                                  ybar0 = mean( Y0 ) )
    ss

    ss2 = a %>% group_by( sssite ) %>% summarise( ATE = mean( Y1 ) - mean( Y0 ),
                                                  ybar1 = mean( Y1 ),
                                                  ybar0 = mean( Y0 ) )
    ss2
    true.site.ATE = mean( ss2$ATE )
    true.site.ATE

    true.indiv.ATE = mean( a$Y1 - a$Y0 )

    true.block.ATE = mean( ss$ATE )

    list( true.site.ATE = true.site.ATE,
          true.block.ATE = true.block.ATE,
          true.indiv.ATE = true.indiv.ATE,
          n.site = nrow(ss2) )
}

test_that("DB estimators work with nested randomization blocks", {

    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( X=1:50, method="big" )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )

    # Our dataset with blocks in sites
    sdat = calc_summary_stats( dat, siteID="siteNo" )
    sdat

    # Finite population, person weighted
    a = estimate_ATE_design_based_from_stats( sdat, siteID="siteID", weight="individual", method="finite" )

    a2 = estimate_ATE_design_based_from_stats( sdat, weight="individual", method="finite" )
    expect_equal( a, a2 )

    # Finite population, site weighted
    # Should be different due to weighting the RA blocks differently.
    b = estimate_ATE_design_based_from_stats( sdat, siteID="siteID", weight="site", method="finite" )
    b2 = estimate_ATE_design_based_from_stats( sdat, weight="site", method="finite" )

    expect_true( b$ATE_hat != b2$ATE_hat )
    expect_true( b$SE != b2$SE )


    # Make sure site weighting is correct
    # (Generate dataset with known ATE to verify.)
    a = make.balanced.dataset()

    ss = a %>% group_by( B ) %>% summarise( ATE = mean( Y1 ) - mean( Y0 ),
                                       ybar1 = mean( Y1 ),
                                       ybar0 = mean( Y0 ) )
    ss
    aa = calc_summary_stats( a, siteID="sssite" )
    aa
    expect_equal( aa$Ybar1, ss$ybar1 )

    site.ATE_hat = estimate_ATE_design_based_from_stats( aa, weight="site", method="finite"  )$ATE_hat
    site.ATE_hat
    expect_equal( 35, site.ATE_hat )

    ss2 = a %>% group_by( sssite ) %>% summarise( ATE = mean( Y1 ) - mean( Y0 ),
                                            ybar1 = mean( Y1 ),
                                            ybar0 = mean( Y0 ) )
    ss2
    true.site.ATE = mean( ss2$ATE )

    corrATEhat = estimate_ATE_design_based_from_stats( aa, siteID="siteID", weight="site", method="finite"  )$ATE_hat
    corrATEhat
    expect_equal( true.site.ATE,
                  corrATEhat )



    # Superpopulation, individual
    est1 = estimate_ATE_design_based_from_stats( aa, siteID="siteID", weight="individual", method="superpop"  )
    est2 = estimate_ATE_design_based_from_stats( aa, weight="individual", method="superpop"  )
    est1
    est2
    expect_true( est1$ATE_hat == est2$ATE_hat )
    expect_true( est1$SE != est2$SE )

    # Superpopulation, site
    est1 = estimate_ATE_design_based_from_stats( aa, siteID="siteID", weight="site", method="superpop"  )
    est2 = estimate_ATE_design_based_from_stats( aa, weight="site", method="superpop"  )
    est1
    est2
    expect_equal( est1$ATE_hat, true.site.ATE )
    expect_equal( est2$ATE_hat, site.ATE_hat )
    expect_true( est1$SE != est2$SE )

    # If we pass siteID but don't actually have nesting
    est3 = estimate_ATE_design_based_from_stats( aa, siteID="B", weight="site", method="superpop"  )
    expect_true( est3$ATE_hat == site.ATE_hat )
    expect_true( est3$SE == est2$SE )

})





test_that("linear regression works with nested randomization blocks", {

    a = make.balanced.dataset()
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = blkvar:::interacted_linear_estimators( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    expect_equal( est1$ATE[[1]], params$true.site.ATE )
    expect_equal( est1$ATE[[2]], params$true.indiv.ATE )

    est1b = blkvar:::interacted_linear_estimators( outcome, Tx, BB, data=a )
    est1b
    expect_equal( est1b$ATE[[1]], params$true.block.ATE )
    expect_equal( est1b$ATE[[2]], params$true.indiv.ATE )

    expect_true( est1b$SE[[1]] != est1$SE[[1]] )
    expect_true( est1b$SE[[2]] == est1$SE[[2]] )


    est2 = blkvar:::fixed_effect_estimators( outcome, Tx, BB, siteID = "sssite", data=a )
    est2
    est2b = blkvar:::fixed_effect_estimators( outcome, Tx, BB, data=a )
    est2b
    expect_true( est2$SE[[3]] > est2b$SE[[3]])
    expect_equal( est2$SE[1:2], est2b$SE[1:2])


})






test_that("weighted linear regression works with nested randomization blocks", {

    a = make.balanced.dataset()
    head( a )
    params = get.params(a)
    params
    table( a$sssite )
    table( a$B, a$sssite )
    mean( a$Y1 - a$Y0 )

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = weighted_linear_estimators( outcome ~ Tx * BB, siteID = "sssite", data=a, scaled.weights = TRUE )

    est1b = weighted_linear_estimators( outcome ~ Tx * BB, data=a, scaled.weights = FALSE )

    estDB = compare_methods( outcome, Tx, BB, data=a, siteID = "sssite", include_MLM = FALSE, include_LM = TRUE,
                             include_block = FALSE )
    est1
    est1b
    estDB

    expect_equal( est1$ATE[[2]], estDB$ATE[[2]] )
    expect_equal( est1$ATE[[2]], estDB$ATE[estDB$method=="DB-FP-Sites"] )
    expect_equal( est1$ATE[[2]], estDB$ATE[estDB$method=="FE-IPTW-Sites"] )

    expect_equal( est1$ATE[[1]], params$true.indiv.ATE )
    expect_equal( est1$ATE[[2]], params$true.site.ATE )

    expect_equal( est1b$ATE[[1]], params$true.indiv.ATE )
    expect_equal( est1b$ATE[[2]], params$true.block.ATE )


    # Another test--match design based?
    set.seed( 1019 )
    dat = generate_blocked_data_obs_linear( X=1:50, method="big" )
    dat$siteNo = round( 1 + as.numeric( dat$B ) / 3 )

    # Our dataset with blocks in sites
    sdat = calc_summary_stats( dat, siteID="siteNo" )

    # Finite population, person weighted
    a = estimate_ATE_design_based_from_stats( sdat, siteID="siteID", weight="site", method="finite" )
    a2 = blkvar:::weighted_linear_estimators( Yobs ~ Z*B, siteID = "siteNo", data=dat )
    expect_equal( a$ATE_hat, a2$ATE[[2]] )

})


test_that("multilevel regression works with nested randomization blocks", {

    a = make.big.balanced.dataset( rps = 5 )
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = blkvar:::compare_MLM_methods( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    est1b = blkvar:::compare_MLM_methods( outcome, Tx, BB, data=a )
    est1b

    expect_true( all( est1[c(1,3),] == est1b[c(1,3),] ) )

    expect_true( est1$ATE[[2]] != est1b$ATE[[2]] )
    expect_true( est1$SE[[2]] > est1b$SE[[2]] )

})


test_that("all linear regression works with nested randomization blocks", {

    a = make.big.balanced.dataset( 5 )
    head( a )
    params = get.params(a)
    params

    a = rename( a, outcome = Yobs, Tx = Z, BB = B )

    est1 = linear_model_estimators( outcome, Tx, BB, siteID = "sssite", data=a )
    est1

    est1b = linear_model_estimators( outcome, Tx, BB, data=a )
    est1b

    res = merge( est1, est1b, by="method", suffixes=c(".site",".block" ) )
    res

    res$change = res$ATE_hat.site - res$ATE_hat.block
    res$changeSE = res$SE.site - res$SE.block

    res

    ff = res [sort(unique(c(which(res$change != 0), which(res$changeSE != 0)))), ]
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

    res$change = res$ATE_hat.site - res$ATE_hat.block
    res$changeSE = res$SE.site - res$SE.block

    res


    ff = res [sort(unique(c(which(res$change != 0), which(res$changeSE != 0)))), ]
    ff

    expect_true( nrow( ff ) == 8 )

})


test_that("compare_methods with adjustment works with nested randomization blocks", {
    a = generate_multilevel_data( n.bar=10, J=10 )
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
                             control_formula = ~ X1 + X2 )
    est1b


    est1 = compare_methods( outcome, Tx, BB, data=a,
                            siteID = "sssite",
                             control_formula = ~ X1 + X2 )
    est1


    res = merge( est1, est1b, by="method", suffixes=c(".site",".block" ) )

    res$change = res$ATE_hat.site - res$ATE_hat.block
    res$changeSE = res$SE.site - res$SE.block

    options( digits = 3 )
    res

    # ff = filter( res, change == 0 & changeSE == 0 )
    # ff

    ff = res [sort(unique(c(which(res$change != 0), which(res$changeSE != 0)))), ]
    ff

    expect_true( nrow( ff ) == 8 )

})





