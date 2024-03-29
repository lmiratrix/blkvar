##
## Testing the data generation code
##

library( testthat )
library( blkvar )

context("Checking multisite DGP functions")

test_that( "Multisite DGP works", {

    # exploring sites
    sdf = generate_multilevel_data( n.bar=10, J=10,
                   rho2.0W = 0.3, rho2.1W = 0.1,
                   tau.11.star = 0.3, return.sites=TRUE )

    expect_equal( nrow(sdf), 10 )

    head(sdf)
    nrow( sdf )
    cov( sdf$beta.0, sdf$beta.1 )
    cov( sdf$u0, sdf$u1 )


    dat = generate_multilevel_data( n.bar=10, J=10,
                   rho2.0W = 0.3, rho2.1W = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )

    head( dat )
    expect_equal( length( unique( dat$sid ) ), 10 )

} )

test_that( "Variable p and n works", {
    set.seed( 1019 )
    df = generate_multilevel_data_no_cov( n.bar=200, J=20,
                         tau.11.star = 0.1^2,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         finite.model = FALSE )

    sites = df %>% group_by( sid ) %>%
        summarise( n = n(),
                   p.Z = mean( Z ),
                   Y.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ) )

    head( sites )
    # ggplot2::qplot( sites$n )
    range( sites$n )
    rat = max( sites$n ) / min( sites$n )
    rat
    # ggplot2::qplot( sites$p.Z )
    mean( sites$p.Z )
    rat.Z = max( sites$p.Z ) / min( sites$p.Z )
    rat.Z
    expect_true( rat > 2 )
    expect_true( rat.Z > 2 )

} )

test_that( "Site impact correlation works", {
    set.seed( 1019 )
    df = generate_multilevel_data_no_cov( n.bar=200, J=200,
                         tau.11.star = 0.3^2,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         size.impact.correlate = TRUE,
                         finite.model = FALSE )
    head( df )
    sites = df %>% group_by( sid ) %>%
        summarise( n = n(),
                   p.Z = mean( Z ),
                   Y.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   Y.true = mean( Y1 ) - mean( Y0 ) )

    head( sites )
    # ggplot2::qplot( sites$n, sites$Y.true )
    cor(  sites$n, sites$Y.true )

    ATE.site = mean( sites$Y.true )
    ATE.site
    ATE.indiv = mean( df$Y1 - df$Y0 )
    ATE.indiv
    expect_true( cor(  sites$n, sites$Y.true ) > 0.50 )

} )



test_that( "prop treated by impact correlation works", {
    set.seed( 1019 )
    df = generate_multilevel_data_no_cov( n.bar=20, J=400,
                         tau.11.star = 0.3^2,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         proptx.impact.correlate = 1,
                         finite.model = FALSE )
    head( df )
    sites = df %>% group_by( sid ) %>%
        summarise( n = n(),
                   p.Z = mean( Z ),
                   Y.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   Y.true = mean( Y1 ) - mean( Y0 ) )

    mean( sites$n == 4 )
    mean( sites$n )
    summary( sites$n )
    # ggplot2::qplot( sites$n )

    head( sites )
    # ggplot2::qplot( sites$n, sites$Y.true )
    cor(  sites$n, sites$Y.true )

    # ggplot2::qplot( sites$p.Z, sites$Y.true )
    cor( sites$p.Z, sites$Y.true )

    expect_true( cor(  sites$p.Z, sites$Y.true ) > 0.50 )
    expect_true( abs( cor(  sites$n, sites$Y.true ) ) < 0.10 )




    df = generate_multilevel_data_no_cov( n.bar=20, J=200,
                         tau.11.star = 0.3^2,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         proptx.impact.correlate = -1,
                         finite.model = FALSE )
    head( df )
    sites = df %>% group_by( sid ) %>%
        summarise( n = n(),
                   p.Z = mean( Z ),
                   Y.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   Y.true = mean( Y1 ) - mean( Y0 ) )

    head( sites )
    # ggplot2::qplot( sites$n, sites$Y.true )
    cor(  sites$n, sites$Y.true )

    # ggplot2::qplot( sites$p.Z, sites$Y.true )
    cor( sites$p.Z, sites$Y.true )

    expect_true( cor(  sites$p.Z, sites$Y.true ) < -0.50 )
    expect_true( abs( cor(  sites$n, sites$Y.true ) ) < 0.15 )
    expect_true( abs( cor(  sites$n, sites$Y.true ) ) > -0.15 )





    df = generate_multilevel_data_no_cov( n.bar=200, J=200,
                         tau.11.star = 0.3^2,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         size.impact.correlate = TRUE,
                         proptx.impact.correlate = TRUE,
                         finite.model = FALSE )
    sites = df %>% group_by( sid ) %>%
        summarise( n = n(),
                   p.Z = mean( Z ),
                   Y.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   Y.true = mean( Y1 ) - mean( Y0 ) )

    head( sites )
    # ggplot2::qplot( sites$n, sites$Y.true )
    cor(  sites$n, sites$Y.true )

    # ggplot2::qplot( sites$p.Z, sites$Y.true )
    cor( sites$p.Z, sites$Y.true )

    expect_true( cor(  sites$p.Z, sites$Y.true ) > 0.50 )
    expect_true( cor(  sites$n, sites$Y.true ) > 0.50 )



    df = generate_multilevel_data_no_cov( n.bar=25, J=80,
                         tau.11.star = 0,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         size.impact.correlate = TRUE,
                         proptx.impact.correlate = TRUE,
                         finite.model = FALSE )
    compare_methods( Yobs, Z, sid, data= df )
    tt = table( df$sid, df$Z )
    head( tt )
    tail( tt )
    sizes = unlist( table( df$sid, df$Z ))
    expect_true( all( sizes >= 2 ) )

} )



test_that( "Bounding of 2 tx and 2 co units works", {
    set.seed( 1019 )
    df = generate_multilevel_data_no_cov( n.bar=16, J=200,
                         tau.11.star = 0.3^2,
                         ICC = 0.20,
                         variable.n = TRUE,
                         variable.p = TRUE,
                         size.impact.correlate = TRUE,
                         finite.model = FALSE )
    head( df )
    sites = df %>% group_by( sid ) %>%
        summarise( n = n(),
                   nT = sum( Z ),
                   nC = sum( 1-Z ),
                   p.Z = mean( Z ),
                   Y.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
                   Y.true = mean( Y1 ) - mean( Y0 ) )

    table( sites$n )
    table( sites$nT )
    table( sites$nC )
    expect_true( all( sites$nT >= 2 ) )
    expect_true( all( sites$nC >= 2 ) )

} )


test_that( "Other DGP calls work", {
    set.seed( 101010 )

    df = generate_multilevel_data_model( 10, J=300, 0.5, 0, 0, 0, 0, 0.3, 0, 0.3, 1, variable.n=FALSE)
    expect_equal( nrow(df), 3000 )
} )

test_that( "Other DGP calls work (#2)", {
    set.seed( 101010 )

    df = generate_multilevel_data( n.bar=10, J=300,
                  tau.11.star = 0.3,
                  verbose=FALSE)

    var( df$Y0 )
    var( df$Y1 )

    M0 = lmer( Yobs ~ 1 + Z + (Z|sid), data=df )
    #display( M0 )

    M1 = lmer( Yobs ~ 1 + W*Z + (Z|sid), data=df )
    #display( M1 )


    sites = df %>% group_by( sid, W ) %>% summarise( Y0.bar = mean( Y0 ),
                                                     Y1.bar = mean( Y1 ),
                                                     beta = mean( Y1 - Y0 ),
                                                     n = n(),
                                                     p.Tx = mean( Z ) )
    nrow( sites )
    rst <- sites %>% ungroup() %>% summarise( mean.Y0 = mean( Y0.bar ),
                                       mean.Y1 = mean( Y1.bar ),
                                       cor.Ys = cor( Y0.bar, Y1.bar ),
                                       cov.Ys = cov( Y0.bar, Y1.bar ),
                                       mean.beta = mean( beta ),
                                       n.bar = mean( n ),
                                       p = mean( p.Tx ),
                                       W.bar = mean( W ) )
    tt = t.test( sites$W )
    tt
    expect_true( tt$p.value > 0.05 )

} )


test_that( "Cluster randomization options work", {
    df = generate_multilevel_data( n.bar=10, J=20,
                  tau.11.star = 0.3,
                  verbose=FALSE,
                  cluster.rand=TRUE)
    tb = table( df$Z, df$sid)
    tb
    expect_true( all( tb[1,] * tb[2,] == 0 ) )
    df = generate_multilevel_data( n.bar=10, J=300,
                  tau.11.star = 0.3,
                  verbose=FALSE,
                  cluster.rand=FALSE)
    tb = table( df$Z, df$sid)
    expect_true( all( tb[1,] * tb[2,] != 0 ) )


    dat = generate_multilevel_data_model( n.bar = 20, J = 4,
                                          p = 0.5,
                                          beta.X = 0.5,
                                          tau.00 = 0.5,
                                          tau.01 = 0,
                                          tau.11 = 0.2,
                                          sigma2.e = 1,
                                          cluster.rand= TRUE,
                                          gamma.00 = 0,
                                          gamma.10 = 0.5,
                                          gamma.01 = NULL,
                                          gamma.11 = NULL )
    expect_true( is.null( dat$W ) )

} )





test_that( "Individual covariate options work", {
    set.seed( 1020 )
    df = generate_multilevel_data_model( n.bar = 20, J = 300,
                        gamma.00 = 1, gamma.01 = 1, gamma.10 = 0.3, gamma.11 = 0,
                        tau.00 = 1, tau.01 = 0, tau.11 = 0,
                        sigma2.e = 1, sigma2.W = 1,
                        beta.X = 0.8, sigma2.mean.X = 0.5, variable.n = FALSE, variable.p = FALSE )
    #head( df )

    sd( df$Y0 )
    sd( df$Y1 )

    d1 = d2 = df
    d1$Yobs = d1$Y0
    d1$Z = 0
    d2$Yobs = d2$Y1
    d2$Z = 1
    dd = bind_rows( d1, d2 )
    M0 = lmer( Yobs ~ 1 + Z + W + X + (1|sid), data=dd )
    #summary( M0 )

    params = c( 1, 0.3, 1, 0.8 )
    CI = confint(M0, method="Wald")[3:6,]
    CI
    CI = cbind( CI, params )
    CI

    expect_true( all( CI[,1] <= params ) )
    expect_true( all( CI[,2] >= params ) )


    M0 = lm( Yobs ~ 1 + Z + W + X, data=df )
    #summary( M0 )

    gp = df %>% group_by( sid ) %>%
        summarise( mean.X = mean( X ) )
    M1 = lmer( X ~ 1 + (1|sid), data=df )
    M1
    #arm::display( M1 )

} )


