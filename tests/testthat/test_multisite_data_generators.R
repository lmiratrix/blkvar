##
## Testing the data generation code
##

library( testthat )
library( blkvar )
context("Checking multisite DGP functions")

test_that( "Multisite DGP works", {

    # exploring sites
    sdf = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=TRUE )

    expect_equal( nrow(sdf), 10 )

    head(sdf)
    nrow( sdf )
    cov( sdf$beta.0, sdf$beta.1 )
    cov( sdf$u0, sdf$u1 )


    dat = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )

    head( dat )
    expect_equal( length( unique( dat$sid ) ), 10 )

} )

test_that( "Other DGP calls work", {
    set.seed( 101010 )

    df = gen.dat.model( 10, J=300, 0.5, 0, 0, 0, 0, 0.3, 0, 0.3, 1, variable.n=FALSE)
    expect_equal( nrow(df), 3000 )

    df = gen.dat( n.bar=10, J=300,
                  tau.11.star = 0.3,
                  verbose=FALSE)

    var( df$Y0 )
    var( df$Y1 )

    M0 = lmer( Yobs ~ 1 + Z + (Z|sid), data=df )
    #display( M0 )

    M1 = lmer( Yobs ~ 1 + X*Z + (Z|sid), data=df )
    #display( M1 )


    sites = df %>% group_by( sid, X ) %>% summarise( Y0.bar = mean( Y0 ),
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
                                       X.bar = mean( X ) )
    expect_true( abs( rst$X.bar ) < 0.05 )

} )

