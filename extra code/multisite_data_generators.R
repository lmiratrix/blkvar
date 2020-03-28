

## Code for generating multi-site randomized trial data with cross-site
## impact variation
##
## Core functions:
##     gen_dat_model - generate according to a model with the coefficients specified
##     gen_dat - generate data according to some higher level parameters such as ICC, which will
##        automatically calculate model coefficients and then call gen_dat_model





block_distn <- function(J, n.bar, size.ratio) {
    N <- 1 + 3 * size.ratio
    p <- (N - 1) / N
    small <- rbinom(J, 1, p)
    Y <- runif(J)
    Y <- n.bar * ifelse(small, Y, Y * (N - 1) + 1)
    Y
}



# For pretty-printing output.  Can do things like, e.g.,
# scat( "This is %d and %f\n", 10, 3.333 )
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


# # # # Some initial variables that one might consider
# if ( FALSE ) {

    # n.bar = 10
    # J = 30
    # p = 0.5
    # tau.11.star = 0.3
    # rho2.0W = 0.8
    # rho2.1W = 0.5
    # ICC = 0
    # gamma.00 = 0
    # gamma.10 = 0.2  # ATE


    # # for debugging
    # gamma.00 = gamma.01 = gamma.10 = gamma.11 = 0
    # tau.00 = tau.11 = 0.3
    # tau.01 = 0
    # sigma2.e = 1
# }






##### Another DGP from Catherine ######

# Catherine's create data
#     a -- average treatment effect
#     t -- tau
#     s -- number of sites
#     n -- number of students per site
catherine_gen_dat <- function(a, t, s, n) {
 #create site ates
 raw_ate_vec <- rnorm(n=s, mean=a,sd=t)
 tau_fp <- sd(raw_ate_vec)
 ate_vec <- rep(raw_ate_vec,each=n)
 #create data frame with sids
 dat <- data.frame(sid=as.factor(rep(1:s,each=n)))
 # randomize within blocks
 dat <- dat %>% group_by(sid) %>% mutate( Z = as.numeric( sample( n() ) <= n()/2 ) ) #50% treatment
 dat$Y0 <- rnorm(n * s, mean = 0, sd = 1)
 dat$Y1 <- dat$Y0 + ate_vec
 dat$Yobs <- dat$Y0 * (1 - dat$Z) + dat$Y1 * (dat$Z)
 dat$tau_fp <- tau_fp
 return(dat)
}


##### Seeing how code works  #####

# if ( FALSE ) {
    # # exploring sites
    # sdf = gen_dat( n.bar=10, J=10,
                   # rho2.0W = 0.3, rho2.1W = 0.1,
                   # tau.11.star = 0.3, return.sites=TRUE )

    # head(sdf)
    # nrow( sdf )
    # cov( sdf$beta.0, sdf$beta.1 )
    # cov( sdf$u0, sdf$u1 )


    # dat = gen_dat( n.bar=10, J=10,
                   # rho2.0W = 0.3, rho2.1W = 0.1,
                   # tau.11.star = 0.3, return.sites=FALSE )

    # head( dat )


    # sdf = gen_dat_no_cov( n.bar=10, J=10,
                          # tau.11.star = 0.3, return.sites=TRUE )
    # head(sdf)
    # nrow( sdf )
    # cov( sdf$beta.0, sdf$beta.1 )
    # cov( sdf$u0, sdf$u1 )
# }

# if ( FALSE ) {
    # df = gen_dat_model( 10, J=300, 0.5, 0, 0, 0, 0, 0.3, 0, 0.3, 1 )

    # df = gen_dat( n.bar=10, J=300,
                  # tau.11.star = 0.3,
                  # verbose=TRUE)

    # var( df$Y0 )
    # var( df$Y1 )

    # M0 = lmer( Yobs ~ 1 + Z + (Z|sid), data=df )
    # display( M0 )

    # M1 = lmer( Yobs ~ 1 + W*Z + (Z|sid), data=df )
    # display( M1 )


    # sites = df %>% group_by( sid, W ) %>% dplyr::( Y0.bar = mean( Y0 ),
                                                     # Y1.bar = mean( Y1 ),
                                                     # beta = mean( Y1 - Y0 ),
                                                     # n = n(),
                                                     # p.Tx = mean( Z ) )
    # nrow( sites )
    # sites %>% ungroup() %>% dplur::summarise( mean.Y0 = mean( Y0.bar ),
                                       # mean.Y1 = mean( Y1.bar ),
                                       # cor.Ys = cor( Y0.bar, Y1.bar ),
                                       # cov.Ys = cov( Y0.bar, Y1.bar ),
                                       # mean.beta = mean( beta ),
                                       # n.bar = mean( n ),
                                       # p = mean( p.Tx ),
                                       # W.bar = mean( W ) )

# }

# ######## A small simulation study  ########

# if ( FALSE ) {

    # single.rho = function( rho2.1W, R = 5 ) {
        # cat( "Doing rho = ", rho2.1W, "\n" )
        # res = plyr::rdply( R, {
            # df = gen_dat( n.bar=10, J=20,
                          # rho2.1W = rho2.1W,
                          # tau.11.star = 0.1 )

            # c( ideosyncratic = analysis.1( df ),
               # systematic = analysis.2( df ),
               # combination = analysis.3( df  ) )
        # })
        # res = reshape2::melt( res, id=".n" )
        # res$rho2.1W = rho2.1W
        # res
    # }

    # rhos = seq( 0, .4, by = 0.1/4 )

    # res = map_df( rhos, single.rho, R=100 )
    # head( res )

    # powers = res %>% group_by( rho2.1W, variable ) %>% dplyr::summarise( power = mean( value <= 0.05 ) )
    # head( powers )

    # ggplot( powers, aes(x=rho2.1W, y=power, col=variable ) ) +
        # geom_smooth() +
        # geom_hline( yintercept = 0.05 )



    # ggplot( subset( res, rho2.1W == 0.2 ), aes( x=value ) ) +
        # facet_wrap( ~ variable, ncol= 1 ) +
        # geom_histogram()

    # print( powers )

# }


# ###### Testing the finite.sample code ##########

# # Testing the finite model code.
# if ( FALSE ) {
    # CANONICAL <- NULL

    # df = gen_dat_no_cov( n.bar=10, J=J,
                         # gamma.10 = 0,
                         # tau.11.star = 0.2^2,
                         # ICC = 0.85,
                         # variable.n = FALSE,
                         # return.sites=TRUE,
                         # finite.model=TRUE
    # )
    # head( df )
    # sd( df$u0 )
    # sd( df$u1 )

    # df2 = gen_dat_no_cov( n.bar=10, J=J,
                          # gamma.10 = 0,
                          # tau.11.star = 0.2^2,
                          # ICC = 0.85,
                          # variable.n = FALSE,
                          # return.sites=TRUE,
                          # finite.model=TRUE

    # )
    # head( df2 )
    # df$u1 - df2$u1

# }

