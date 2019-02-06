

## Code for generating multi-site randomized trial data with cross-site
## impact variation
##
## Core functions:
##     gen.dat.model - generate according to a model with the coefficients specified
##     gen.dat - generate data according to some higher level parameters such as ICC, which will
##        automatically calculate model coefficients and then call gen.dat.model


library( tidyverse )
library( MASS )
library( lme4 )
library( arm )
library( lmtest )

# This is a hack so if the 'finite.model' flag is set to true, it will save the
# randomly generated units in this variable and simply rerandomize them with
# subsequent calls to gen.dat.
CANONICAL = NULL


# For pretty-printing output.  Can do things like, e.g.,
# scat( "This is %d and %f\n", 10, 3.333 )
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


# Some initial variables that one might consider
if ( FALSE ) {

    n.bar = 10
    J = 30
    p = 0.5
    tau.11.star = 0.3
    rho2.0X = 0.8
    rho2.1X = 0.5
    ICC = 0
    gamma.00 = 0
    gamma.10 = 0.2  # ATE


    # for debugging
    gamma.00 = gamma.01 = gamma.10 = gamma.11 = 0
    tau.00 = tau.11 = 0.3
    tau.01 = 0
    sigma2.e = 1
}


# Utility function to describe different characteristics of the distribution of
# sites (given all the potential outcomes).
describe.data = function( data, Y0="Y0", Y1="Y1", Z="Z", sid="sid" ) {
    #old param list: Y0, Y1, Z, sid, data = NULL
    data = rename( data, Y0 = !!rlang::sym(Y0),
                   Y1 = !!rlang::sym(Y1),
                   Z = !!rlang::sym(Z),
                   sid = !!rlang::sym(sid) )

    sites = data %>% group_by( !!rlang::sym(sid) ) %>% summarise( Y0.bar = mean( Y0 ),
                                                  Y1.bar = mean( Y1 ),
                                                  beta = Y1.bar - Y0.bar,
                                                  n = n(),
                                                  p.Tx = mean( Z ) )

    site.var = sd( sites$beta )
    scat( "\n\nDescription of site distribution:" )
    scat( "\n\tsite average variation:\tsd(Y0.bar) = %.3f\tsd(Y1.bar) = %.3f", sd( sites$Y0.bar ), sd( sites$Y1.bar ) )
    scat( "\n\tmarginal variation:\tsd(Y0) = %.3f\t\tsd(Y1) = %.3f\t\tcor(Y0,Y1) = %.3f", sd( data$Y0 ), sd( data$Y1 ), cor( data$Y0, data$Y1 ) )
    scat( "\n\tsd( beta ) = %.3f\n\tcorr( Y0.bar, Y1.bar ) = %.3f\tcov = %.3f", site.var, cor( sites$Y0.bar, sites$Y1.bar ),
          cov( sites$Y0.bar, sites$Y1.bar ))
    scat( "\n\tcor( Y0.bar, beta ) = %.3f", cor( sites$Y0.bar, sites$beta ) )
    scat( "\n\tmean( n ) = %.3f\tsd( n ) = %.3f", mean( sites$n ), sd( sites$n ) )
    scat( "\n\tmean( Z.bar ) = %.3f\t(%.3f)\n", mean( sites$p.Tx ), sd( sites$p.Tx) )


    invisible( sites )
}


#' @title Generate multilevel data from a model
#'
#' @description Given a 2-level model, generate data to specifications
#'
#' @param n.bar average site size, Default: 10
#' @param J number sites, Default: 30
#' @param p prop treated, Default: 0.5
#' @param tau.11.star Total amount of cross site treatment variation, both
#'   explained by covariate and not, Default: 0.3
#' @param rho2.0X Explanatory power of X for control outcomes, Default: 0.1
#' @param rho2.1X Explanatory power of X for average treatment impact, Default:
#'   0.5
#' @param ICC The ICC, Default: 0.7
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param zero.corr TRUE means treatment impact and mean site outcome are not
#'   correlated.  TRUE means they are negatively correlated to make the variance
#'   of the treatment group 1, Default: FALSE
#' @param gamma.01 Coefficient for X to site control mean
#' @param gamma.11 Coefficient for X to treatment impact
#' @param tau.00 Variance of site conrol means
#' @param tau.01 Correlation of treatment impact and mean site outcome under control
#' @param tau.11 Treatment impact variance
#' @param sigma2.e Residual standard error
#' @param variable.n Allow n to vary around n.bar, Default: TRUE
#' @param return.sites Return sites, not individual students, Default: FALSE
#' @param cluster.rand TRUE means cluster-randomized.  FALSE means randomized within site.
#'
#' @return Dataframe of data!
#' @rdname gen.dat.model
#' @export
gen.dat.model  = function( n.bar = 10,
                             J = 30,
                             p = 0.5,
                             gamma.00, gamma.01, gamma.10, gamma.11,
                             tau.00, tau.01, tau.11,
                             sigma2.e,
                             variable.n = TRUE,
                             cluster.rand = FALSE,
                             return.sites=FALSE,
                             verbose = FALSE ) {
    require( tidyverse )

    if ( verbose ) {
        scat( "gammas:\t%.2f\t%.2f (%.2f)\n\t%.2f\t%.2f (%.2f)\n",
              gamma.00, gamma.01, gamma.01^2, gamma.10, gamma.11, gamma.11^2 )
        scat( "taus:\t%.2f\t%.2f\tsd=%.2f\n\t%.2f\t%.2f\tsd=%.2f\tcor=%.2f\n",
              tau.00, tau.01, sqrt(tau.00), tau.01, tau.11, sqrt(tau.11), tau.01 / sqrt( tau.00 * tau.11) )
    }

    # site level
    if ( variable.n ) {
        nj = rpois( J, n.bar)
    } else {
        nj = rep( n.bar, J )
    }
    Xj = rnorm( J )
    Sigma = matrix( c( tau.00, tau.01, tau.01, tau.11 ), nrow=2 )
    mv = MASS::mvrnorm( J, c( 0, 0 ), Sigma )
    beta.0j = gamma.00 + gamma.01 * Xj + mv[,1]
    beta.1j = gamma.10 + gamma.11 * Xj + mv[,2]

    if ( cluster.rand ) {
        Zj = 0 + (sample( J ) <= J * p)
    }

    if ( return.sites ) {
        if ( cluster.rand ) {
            data.frame( n = nj, X = Xj, beta.0 = beta.0j, beta.1 = beta.1j, u0 = mv[,1], u1 = mv[,2], Zj = Zj )
        } else {
            data.frame( n = nj, X = Xj, beta.0 = beta.0j, beta.1 = beta.1j, u0 = mv[,1], u1 = mv[,2] )
        }
    } else {    # generate individual level outcomes
        sid = rep( 1:J, nj )

        N = sum( nj )

        if ( cluster.rand ) {
            Zij = rep( Zj, nj )
        } else {
            #Zij = as.numeric( sample( N ) <= N*p )
            dd = data.frame( sid = as.factor(sid) )
            dd = dd %>% group_by( sid ) %>%
                mutate( Z = as.numeric( sample( n() ) <= p*n() ) )
            Zij = dd$Z
        }

        e = rnorm( N, mean=0, sd=sqrt( sigma2.e ) )
        Y0 = beta.0j[sid] + e
        Y1 = beta.0j[sid] + beta.1j[sid] + e

        data.frame( sid=as.factor(sid), X = Xj[sid], Y0 = Y0, Y1 = Y1, Z = Zij, Yobs = ifelse( Zij, Y1, Y0 ) )
    }
}




#' Generate multilevel data with no covariates.  A simplification of function
#' above.
#'
#' @inheritParams gen.dat
#' @param size.impact.correlate    TRUE/FALSE: Are site impacts correlated with site size?
#' @param proptx.impact.correlate  TRUE/FALSE: Are proportion of units treated correlated with site size?
#'
#' @return multisite data with at least 2 tx and 2 co units in each block.
#' @export
gen.dat.model.no.cov = function( n.bar = 16,
                                 J = 30,
                                 p = 0.5,
                                 gamma.00, gamma.10,
                                 tau.00, tau.01, tau.11,
                                 sigma2.e,
                                 variable.n = TRUE,
                                 variable.p = FALSE,
                                 return.sites=FALSE,
                                 finite.model=FALSE,
                                 size.impact.correlate = FALSE,
                                 proptx.impact.correlate = FALSE,
                                 verbose = FALSE ) {
    require( tidyverse )

    # generate site sizes (all the same or different sizes)
    if ( variable.n ) {
        #nj = rpois( J, n.bar)
        nj = round( n.bar * runif( J, 0.25, 1.75 ) )
        if ( any( nj < 4 ) ) {
            warning( "Some sites have fewer than 4 units, disallowing 2 tx and 2 co units" )
        }
    } else {
        nj = rep( n.bar, J )
    }

    # Generate average control outcome and average ATE for all sites
    Sigma = matrix( c( tau.00, tau.01, tau.01, tau.11 ), nrow=2 )

    if ( finite.model ) {
        # Make a canonical set of site charactaristics and then never change
        # them until J changes.
        if ( is.null( CANONICAL ) || nrow( CANONICAL ) != J ) {
            CC <- MASS::mvrnorm( J, c( 0, 0 ), Sigma )
            if ( var( CC[,2] ) > 0 ) {
                CC[,2] = sqrt( tau.11 ) * CC[,2] / sd( CC[,2] )
            }
            if ( var( CC[,1]) > 0 ) {
                CC[,1] = sqrt( tau.00 ) * CC[,1] / sd( CC[,1] )
            }
            CANONICAL <<- CC
        }
        mv = CANONICAL
    } else {
        mv <- MASS::mvrnorm( J, c( 0, 0 ), Sigma )
    }

    if ( size.impact.correlate || proptx.impact.correlate ) {
        mv = mv[ order( mv[,2] + rnorm( J, sd = 0.75 * sqrt( tau.11 ) )), ]
    }

    if ( size.impact.correlate ) {
        if ( !variable.n ) {
            warning( "Can't have correlated site size with site impact without variable sized sites" )
        } else {
            nj = sort( nj )
        }
    }

    # Calculate site intercept and average impacts
    beta.0j = gamma.00 + mv[,1]
    beta.1j = gamma.10 + mv[,2]

    if ( return.sites ) {
        data.frame( n = nj, beta.0 = beta.0j, beta.1 = beta.1j, u0 = mv[,1], u1 = mv[,2] )
    } else {
        # generate individual level outcomes

        sid = as.factor( rep( 1:J, nj ) )

        N = sum( nj )

        # randomize units within each site (proportion p to treatment)
        dd = data.frame( sid = sid  )

        if ( variable.p ) {

            ths = min( p, 1-p ) * 0.75
            ps =  runif( J, p - ths, p + ths )

            if ( proptx.impact.correlate ) {
                ps = sort( ps )
            }

            #threshold to ensure we always have 2 units in tx and co for all blocks
            # regardless of assigned p
            ps = pmin( (nj-2)/nj, pmax( 2/nj, ps ) )

            dd$p = ps[ as.numeric( dd$sid ) ]
        } else if ( proptx.impact.correlate ) {
            warning( "Can't have proportion treated correlated with impact if the proportion treated does not vary." )
        }
        dd = dd %>% group_by( sid ) %>%
            mutate( Z = as.numeric( sample( n() ) <= n()*p ) )
        Zij = dd$Z
        dd$p = NULL

        # individual residuals
        e = rnorm( N, mean=0, sd=sqrt( sigma2.e ) )
        Y0 = beta.0j[sid] + e
        Y1 = beta.0j[sid] + beta.1j[sid] + e

        df <- data.frame( sid=sid,
                          Y0 = Y0, Y1 = Y1,
                          Z = Zij,
                          Yobs = ifelse( Zij, Y1, Y0 ) )
        attr( df, "tau.S" ) <- sd( beta.1j )
        df
    }
}

#' Rerandomize a given multisite simulated dataset and recalulate observed
#' outcomes
#'
#' Useful for simulations of finite sample inference where the dataset should be
#' held static and only randomization is necessary.
#'
#' @param dat Dataframe from, e.g., gen.dat()
#'
#' @return Same dataframe with treatment shuffled and Yobs recalculated.
#' @export
rerandomize.data = function( dat ) {
    dat = dat %>% group_by( sid ) %>%
        mutate( Z = sample( Z ) ) %>% ungroup()
    dat = mutate( dat, Yobs = ifelse( Z, Y1, Y0 ) )
    dat
}


#' @title gen.dat
#' @description Generate data for a multisite trial with a given collection of
#'   features.
#' @param n.bar average site size, Default: 10
#' @param J number sites, Default: 30
#' @param p prop treated, Default: 0.5
#' @param tau.11.star Total amount of cross site treatment variation, both
#'   explained by covariate and not, Default: 0.3
#' @param rho2.0X Explanatory power of X for control outcomes, Default: 0.1
#' @param rho2.1X Explanatory power of X for average treatment impact, Default:
#'   0.5
#' @param ICC The ICC, Default: 0.7
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param zero.corr TRUE means treatment impact and mean site outcome are not
#'   correlated.  TRUE means they are negatively correlated to make the variance
#'   of the treatment group 1, Default: FALSE
#' @param ... Further parameters passed to gen.dat.model()
#' @return Dataframe of data!
#' @rdname gen.dat
#' @export
gen.dat = function( n.bar = 10,
                    J = 30,
                    p = 0.5,
                    tau.11.star = 0.3,
                    rho2.0X = 0.1,
                    rho2.1X = 0.5,
                    ICC = 0.7,
                    gamma.00 = 0,
                    gamma.10 = 0.2,
                    verbose = FALSE,
                    zero.corr = FALSE,
                    ... ) {
    sigma2.X = 1

    gamma.01 = sqrt( rho2.0X * ICC / sigma2.X )
    gamma.11 = sqrt( rho2.1X * tau.11.star )
    tau.00 = (1 - rho2.0X) * ICC
    if ( zero.corr ) {
        tau.01 = 0
    } else {
        tau.01 = - tau.11.star / 2 - gamma.01 * gamma.11 * sigma2.X
    }
    tau.11 = (1 - rho2.1X) * tau.11.star
    sigma2.e = 1 - ICC

    if ( verbose ) {
        scat( "tau.11* = %.2f\tICC = %.2f\trho2.Xs = %.2f, %.2f\n", tau.11.star, ICC, rho2.0X, rho2.1X )
        scat( "tau.00* = %.2f\n", gamma.01^2 * sigma2.X + tau.00 )
        scat( "tau.11* = %.2f\n", gamma.11^2 * sigma2.X + tau.11 )
        scat( "sigma2.e* = %.2f\n", sigma2.e )
    }
    gen.dat.model( n.bar=n.bar, J=J, p=p,
                   gamma.00, gamma.01, gamma.10, gamma.11,
                   tau.00, tau.01, tau.11,
                   sigma2.e,
                   verbose = verbose,
                   ... )
}


#'
#' @title  Simplified version of gen.dat() with no X covariate.
#' @description Generate fake data for simulation studies
#' @inheritParams gen.dat
#' @param control.sd.Y1 Make correlation of random intercept and random slope
#'   such that the variance of the Y1s is 1.0, Default: TRUE
#' @return Dataframe of individual data from a MLM DGP.
#' @rdname gen.dat.no.cov
#' @export
gen.dat.no.cov = function( n.bar = 10,
                           J = 30,
                           p = 0.5,
                           tau.11.star = 0.3,
                           ICC = 0.7,
                           gamma.00 = 0,
                           gamma.10 = 0.2,
                           verbose = FALSE,
                           variable.n = TRUE,
                           control.sd.Y1 = TRUE,
                           ... ) {
    sigma2.X = 1

    tau.00 =  ICC
    tau.11 = tau.11.star
    if ( control.sd.Y1 ) {
        tau.01 = -tau.11 / 2
    } else {
        tau.01 = 0
    }
    sigma2.e = 1 - ICC

    if ( verbose ) {
        scat( "tau.11* = %.2f\tICC = %.2f\n",
              tau.11.star, ICC )
        scat( "tau.00* = %.2f\n",  tau.00 )
        scat( "tau.11* = %.2f\n",  tau.11 )
        scat( "sigma2.e* = %.2f\n", sigma2.e )
    }

    gen.dat.model.no.cov( n.bar=n.bar, J=J, p=p,
                          gamma.00, gamma.10,
                          tau.00, tau.01, tau.11,
                          sigma2.e,
                          verbose = verbose,
                          variable.n = variable.n,
                          ... )
}

##### Another DGP from Catherine ######

# Catherine's create data
#     a -- average treatment effect
#     t -- tau
#     s -- number of sites
#     n -- number of students per site
catherine.gen.dat <- function(a,t,s,n){
    #create site ates
    raw_ate_vec <- rnorm(n=s, mean=a,sd=t)
    tau_fp <- sd(raw_ate_vec)
    ate_vec <- rep(raw_ate_vec,each=n)

    #create data frame with sids
    dat <- data.frame(sid=as.factor(rep(1:s,each=n)))
    #randomize within blocks
    dat = dat %>% group_by(sid) %>%
        mutate( Z = as.numeric( sample( n() ) <= n()/2 ) ) #50% treatment
    dat$Y0 <- rnorm(n*s,mean=0,sd=1)
    dat$Y1 <- dat$Y0+ate_vec
    dat$Yobs <- dat$Y0*(1-dat$Z) + dat$Y1*(dat$Z)

    dat$tau_fp <- tau_fp

    return(dat)
}



##### Seeing how code works  #####

if ( FALSE ) {
    # exploring sites
    sdf = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=TRUE )

    head(sdf)
    nrow( sdf )
    cov( sdf$beta.0, sdf$beta.1 )
    cov( sdf$u0, sdf$u1 )


    dat = gen.dat( n.bar=10, J=10,
                   rho2.0X = 0.3, rho2.1X = 0.1,
                   tau.11.star = 0.3, return.sites=FALSE )

    head( dat )


    sdf = gen.dat.no.cov( n.bar=10, J=10,
                          tau.11.star = 0.3, return.sites=TRUE )
    head(sdf)
    nrow( sdf )
    cov( sdf$beta.0, sdf$beta.1 )
    cov( sdf$u0, sdf$u1 )
}

if ( FALSE ) {
    df = gen.dat.model( 10, J=300, 0.5, 0, 0, 0, 0, 0.3, 0, 0.3, 1 )

    df = gen.dat( n.bar=10, J=300,
                  tau.11.star = 0.3,
                  verbose=TRUE)

    var( df$Y0 )
    var( df$Y1 )

    M0 = lmer( Yobs ~ 1 + Z + (Z|sid), data=df )
    display( M0 )

    M1 = lmer( Yobs ~ 1 + X*Z + (Z|sid), data=df )
    display( M1 )


    sites = df %>% group_by( sid, X ) %>% summarise( Y0.bar = mean( Y0 ),
                                                     Y1.bar = mean( Y1 ),
                                                     beta = mean( Y1 - Y0 ),
                                                     n = n(),
                                                     p.Tx = mean( Z ) )
    nrow( sites )
    sites %>% ungroup() %>% summarise( mean.Y0 = mean( Y0.bar ),
                                       mean.Y1 = mean( Y1.bar ),
                                       cor.Ys = cor( Y0.bar, Y1.bar ),
                                       cov.Ys = cov( Y0.bar, Y1.bar ),
                                       mean.beta = mean( beta ),
                                       n.bar = mean( n ),
                                       p = mean( p.Tx ),
                                       X.bar = mean( X ) )

}

######## A small simulation study  ########

if ( FALSE ) {

    single.rho = function( rho2.1X, R = 5 ) {
        cat( "Doing rho = ", rho2.1X, "\n" )
        res = plyr::rdply( R, {
            df = gen.dat( n.bar=10, J=20,
                          rho2.1X = rho2.1X,
                          tau.11.star = 0.1 )

            c( ideosyncratic = analysis.1( df ),
               systematic = analysis.2( df ),
               combination = analysis.3( df  ) )
        })
        res = reshape2::melt( res, id=".n" )
        res$rho2.1X = rho2.1X
        res
    }

    rhos = seq( 0, .4, by = 0.1/4 )

    res = map_df( rhos, single.rho, R=100 )
    head( res )

    powers = res %>% group_by( rho2.1X, variable ) %>% summarise( power = mean( value <= 0.05 ) )
    head( powers )

    ggplot( powers, aes(x=rho2.1X, y=power, col=variable ) ) +
        geom_smooth() +
        geom_hline( yintercept = 0.05 )



    ggplot( subset( res, rho2.1X == 0.2 ), aes( x=value ) ) +
        facet_wrap( ~ variable, ncol= 1 ) +
        geom_histogram()

    print( powers )

}


###### Testing the finite.sample code ##########


# Testing the finite model code.
if ( FALSE ) {
    CANONICAL <- NULL

    df = gen.dat.no.cov( n.bar=10, J=J,
                         gamma.10 = 0,
                         tau.11.star = 0.2^2,
                         ICC = 0.85,
                         variable.n = FALSE,
                         return.sites=TRUE,
                         finite.model=TRUE
    )
    head( df )
    sd( df$u0 )
    sd( df$u1 )

    df2 = gen.dat.no.cov( n.bar=10, J=J,
                          gamma.10 = 0,
                          tau.11.star = 0.2^2,
                          ICC = 0.85,
                          variable.n = FALSE,
                          return.sites=TRUE,
                          finite.model=TRUE

    )
    head( df2 )
    df$u1 - df2$u1

}

