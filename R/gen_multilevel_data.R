


# This global variable is a hack so if the 'finite.model' flag in
# generate_multilevel_data_model() is set to TRUE, that method will save the
# first randomly generated units in this variable and then, in subsequent calls
# to generate_multilevel_data_model() simply rerandomize these units.
# (generate_multilevel_data_model will check if J has changed and reset if the
# parameter J is different).
.MULTISITE_CANONICAL = NULL

# TODO: BETTER VERSION - Make a "site seed" that gets passed instead.
# Then get a new seed from the current setting, set the site seed and
# generate the finite model, and then set and use the new seed for the
# randomization.



#' Generate site sizes for simulation of multisite data
#'
#' Generate piecewise uniform distribution with a mean of n.bar and a
#' 'size.ratio' that controls the variance of site sizes.
#'
#' @param J number of sites to generate.
#' @param n.bar Average site size.
#' @param size.ratio Index of
#' @return Numeric vector of site sizes.
#' @examples
#' block_distn( 4, 10, 1 )
#' @export
block_distn <- function(J, n.bar, size.ratio) {
    N <- 1 + 3 * size.ratio
    p <- (N - 1) / N
    small <- rbinom(J, 1, p)
    Y <- runif(J)
    Y <- n.bar * ifelse(small, Y, Y * (N - 1) + 1)
    round( Y )
}



#' Rerandomize a given multisite simulated dataset and recalulate observed
#' outcomes
#'
#' Permute the treatment vector Z within site ID.  Useful for simulations of
#' finite sample inference where the dataset should be held static and only
#' randomization is necessary.
#'
#' @param dat Dataframe from, e.g., generate_multilevel_data() that has two columns, 'sid', and 'Z'.
#'
#' @return Same dataframe with treatment shuffled and Yobs recalculated.
#' @export
rerandomize_data <- function(dat) {
    dat <- dat %>% group_by(sid) %>% mutate(Z = sample(Z)) %>% ungroup()
    dat <- mutate(dat, Yobs = ifelse(Z, Y1, Y0))
    dat
}




#' Generate individual data given a dataframe of site level
#' characteristics.
#'
#' @param sdat Dataframe of site level characteristics
#'
#' @describeIn generate_multilevel_data_model Part of data generation
#'   that generates individual level covariates.
#'
#' @return Dataframe of indiivdual level data
#' @export
#'
generate_individual_data = function( sdat, p = 0.5,
                                     sigma2.e = 1,
                                     sigma2.X = 1,
                                     beta.X = NULL,
                                     variable.p = FALSE,
                                     cluster.rand = FALSE,
                                     sigma2.mean.X = 0,
                                     proptx.impact.correlate = FALSE,
                                     verbose = FALSE) {

    include_X <- FALSE
    if (!is.null(beta.X)) {
        include_X <- TRUE
    }

    J = nrow( sdat )
    nj = sdat$n
    beta.0j = sdat$beta.0
    beta.1j = sdat$beta.1
    Wj = sdat$W

    Zj = NULL
    if ( cluster.rand ) {
        stopifnot( !is.null( sdat$Z ) )
        Zj = sdat$Z
    }

    # generate individual level outcomes
    sid <- as.factor(rep(1:J, nj))
    N <- sum(nj)
    # randomize units within each site (proportion p to treatment)
    dd <- data.frame(sid = sid)
    if (variable.p) {
        ths <- min(p, 1 - p) * 0.75
        ps <- runif(J, p - ths, p + ths)

        if (proptx.impact.correlate != 0) {
            ps <- sort( ps, decreasing = (proptx.impact.correlate < 0))
        }

        # threshold to ensure we always have 2 units in tx and co for all blocks
        # regardless of assigned p
        ps <- pmin((nj - 2) / nj, pmax(2 / nj, ps))
        dd$p <- ps[as.numeric(dd$sid)]
    } else if (proptx.impact.correlate != 0) {
        warning( "Can't have proportion treated correlated with impact if the proportion treated does not vary." )
    }

    if (cluster.rand) {
        dd$Z <- Zj[ dd$sid ]
    } else {
        dd <- dd %>% group_by( sid ) %>%
            mutate(Z = as.numeric(sample(n()) <= n() * p))
    }

    Zij <- dd$Z
    dd$p <- NULL

    # individual covariates and residuals
    if (include_X) {
        stopifnot(sigma2.mean.X < sigma2.X)
        stopifnot(sigma2.e - beta.X ^ 2 > 0)
        Xbar <- rnorm(J, mean = 0, sd = sqrt(sigma2.mean.X))
        X <- Xbar[sid] + rnorm(N, mean = 0, sd = sqrt(sigma2.X - sigma2.mean.X))
        e <- rnorm(N, mean = 0, sd = sqrt(sigma2.e - beta.X^2))
        Y0 <- beta.0j[sid] + beta.X * X + e
    } else {
        e <- rnorm(N, mean = 0, sd = sqrt(sigma2.e))
        Y0 <- beta.0j[sid] + e
    }
    Y1 <- Y0 + beta.1j[sid]
    df <- data.frame(sid = as.factor(sid),
                     Y0 = Y0, Y1 = Y1, Z = Zij,
                     Yobs = ifelse(Zij, Y1, Y0))
    if (include_X) {
        df$X <- X
    }
    if ( !is.null( Wj ) ) {
        df$W <- Wj[df$sid]
    }
    df
}



#' @title Generate multilevel data from a model
#'
#' @description Given a 2-level model, generate data to specifications
#'
#'   Model has site-level covariate W and individual-level covariate
#'   X.
#'
#' @param n.bar average site size, Default: 10
#' @param J number sites, Default: 30
#' @param p prop treated, Default: 0.5
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param gamma.01 Coefficient for W to site control mean
#' @param gamma.11 Coefficient for W to treatment impact
#' @param tau.00 Variance of site control means
#' @param tau.01 Covariance of treatment impact and mean site outcome
#'   under control
#' @param tau.11 Treatment impact variance
#' @param sigma2.e Residual standard error
#' @param sigma2.W The variation of the site-level covariate.
#' @param beta.X Coefficient for the individual-level X covariate.  NA
#'   means no covariate.
#' @param sigma2.mean.X How much the individual-level X covariate
#'   means vary across site.
#' @param variable.n Allow n to vary around n.bar, Default: TRUE
#' @param return.sites Return sites, not individual students, Default:
#'   FALSE
#' @param verbose Say stuff while making data?, Default: FALSE
#' @param cluster.rand TRUE means cluster-randomized.  FALSE means
#'   randomized within site.
#' @param size.impact.correlate    Takes values of -1, 0, or 1: Are
#'   site impacts negatively correlated, uncorrelated, or positively
#'   correlated with site size?
#' @param proptx.impact.correlate  Takes values of -1, 0, or 1: Are
#'   proportion of units treated negatively correlated, uncorrelated,
#'   or positively correlated with site size?
#' @param correlate.strength In [0,1], and describes how correlated
#'   the ranking of site impacts will be with proptx and site size, if
#'   they are set to be correlated.
#' @param variable.p Should the proportion of units treated in each
#'   site vary?  Yes/No.
#' @param finite.model If TRUE use a canonical set of random site
#'   effects.  When TRUE this method will save the multivariate normal
#'   draw and reuse it in subsequent calls to
#'   generate_multilevel_data_model until a call with a different J is
#'   made.  Recommended to use FALSE.
#' @param size.ratio The degree to which the site sizes should vary,
#'   if they should vary.
#' @param site.sizes (Optional) vector of manually specified site sizes.
#'   If not specified, use n.bar and variable.n to generate site sizes.
#'
#' @return Dataframe of individual level data (unless
#'   return.sites=TRUE)!  Dataframe has treatment column, outcome
#'   column, covariates, and block IDs.
#' @importFrom magrittr '%>%'
#' @export
generate_multilevel_data_model <- function(n.bar = 10, J = 30, p = 0.5,
                                           gamma.00, gamma.01, gamma.10, gamma.11,
                                           tau.00, tau.01, tau.11, sigma2.e, sigma2.W = 1,
                                           beta.X = NULL, sigma2.mean.X = 0,
                                           variable.n = TRUE, variable.p = FALSE,
                                           site.sizes = NULL,
                                           cluster.rand = FALSE, return.sites = FALSE,
                                           finite.model = FALSE,
                                           size.impact.correlate = 0, proptx.impact.correlate = 0,
                                           correlate.strength = 0.75, size.ratio = 1 / 3, verbose = FALSE) {

    stopifnot(size.impact.correlate %in% c(-1, 0, 1))
    stopifnot(proptx.impact.correlate %in% c(-1, 0, 1))
    if (verbose) {
        scat( "gammas:\t%.2f\t%.2f (%.2f)\n\t%.2f\t%.2f (%.2f)\n",
              gamma.00, gamma.01, gamma.01 ^ 2, gamma.10, gamma.11, gamma.11 ^ 2)
        scat( "taus:\t%.2f\t%.2f\tsd=%.2f\n\t%.2f\t%.2f\tsd=%.2f\tcor=%.2f\n",
              tau.00, tau.01, sqrt(tau.00), tau.01, tau.11, sqrt(tau.11), tau.01 / sqrt(tau.00 * tau.11))
        scat( "beta.X:\t%.2f\n", beta.X)
    }

    # generate site sizes (all the same or different sizes)
    if (is.null(site.sizes)) {
      if (variable.n) {
        stopifnot(n.bar > 4)
        # nj = rpois( J, n.bar)
        # nj = round( n.bar * runif( J, 0.25, 1.75 ) )
        nj <- 4 + block_distn(J, n.bar - 4, size.ratio)
        nj[ nj < 4 ] <- 4
        # if ( any( nj < 4 ) ) {
        #    warning( "Some sites have fewer than 4 units, disallowing 2 tx and 2 co units" )
        #}
      } else {
        nj <- rep(n.bar, J)
      }
    } else {
      stopifnot(length(site.sizes) == J)                # specify size for each site
      stopifnot(min(site.sizes) >= 4)                   # sites must all be large enough
      stopifnot(all(round(site.sizes) == site.sizes))   # integer-valued site sizes

      nj <- site.sizes
    }

    include_W = !is.null( gamma.10 ) && !is.null(gamma.11)
    if ( include_W ) {
        Wj <- rnorm(J, mean = 0, sd = sqrt(sigma2.W))
    } else {
        Wj <- rep( 0, J )
        gamma.01 = 0
        gamma.11 = 0
    }



    # Generate average control outcome and average ATE for all sites
    Sigma <- matrix(c(tau.00, tau.01, tau.01, tau.11), nrow = 2)
    if (finite.model) {
        # Make a canonical set of site charactaristics and then never change
        # them until J changes.
        if (is.null(.MULTISITE_CANONICAL) || nrow(.MULTISITE_CANONICAL) != J) {
            CC <- MASS::mvrnorm(J, c(0, 0), Sigma)
            if (var(CC[, 2]) > 0) {
                CC[, 2] <- sqrt(tau.11) * CC[, 2] / sd(CC[, 2])
            }
            if (var(CC[, 1]) > 0 ) {
                CC[, 1] <- sqrt(tau.00) * CC[, 1] / sd(CC[, 1])
            }
            .MULTISITE_CANONICAL <<- CC
        }
        mv <- .MULTISITE_CANONICAL
    } else {
        mv <- MASS::mvrnorm(J, c(0, 0), Sigma)
    }

    # If we want a relationship of impact to site size or proportion treated,
    # sort our random effects from small to large (with some noise so the ordering is not perfect).
    if ( (size.impact.correlate != 0) || (proptx.impact.correlate != 0) ) {
        mv <- mv[ order( correlate.strength ^ 2 * mv[, 2] + (1 - correlate.strength ^ 2) * rnorm(J, sd = sqrt(tau.11))), ]
    }
    if (size.impact.correlate != 0) {
        if (!variable.n) {
            warning( "Can't have correlated site size with site impact without variable sized sites" )
        } else {
            nj <- sort(nj, decreasing = (size.impact.correlate < 0))
        }
    }

    # Calculate site intercept and average impacts
    beta.0j <- gamma.00 + gamma.01 * Wj + mv[, 1]
    beta.1j <- gamma.10 + gamma.11 * Wj + mv[, 2]

    if (cluster.rand) {
        stopifnot( !variable.p )
        Zj <- 0 + (sample(J) <= J * p)
    }

    if (cluster.rand) {
        df <- data.frame(n = nj, W = Wj, beta.0 = beta.0j, beta.1 = beta.1j,
                         u0 = mv[,1], u1 = mv[,2], Zj = Zj)
    } else {
        df <- data.frame(n = nj, W = Wj, beta.0 = beta.0j, beta.1 = beta.1j,
                         u0 = mv[,1], u1 = mv[,2])
    }

    if (!include_W) {
        df$W <- NULL
    }

    if (return.sites) {
        return( df )
    } else {
        df <- generate_individual_data( df, p = p, sigma2.e = sigma2.e, beta.X = beta.X,
                                        proptx.impact.correlate = proptx.impact.correlate,
                                        sigma2.mean.X = sigma2.mean.X,
                                        variable.p = variable.p, cluster.rand = cluster.rand,
                                        verbose = verose )

        attr(df, "tau.S") <- sd( beta.1j)

        return( df)
    }
}






#' @title generate_multilevel_data
#'
#' @describeIn generate_multilevel_data_model Wrapper for
#'   generate_multilevel_data_model that rescales parameters to make
#'   standardization easier.
#'
#' @param tau.11.star Total amount of cross site treatment variation,
#'   both explained by covariate and not, Default: 0.3
#' @param rho2.0W Explanatory power of W for control outcomes,
#'   Default: 0.1
#' @param rho2.1W Explanatory power of W for average treatment impact,
#'   Default: 0.5
#' @param ICC The ICC, Default: 0.7
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param zero.corr TRUE means treatment impact and mean site outcome
#'   are not correlated.  TRUE means they are negatively correlated to
#'   make the variance of the treatment group 1, Default: FALSE
#' @param ... Further parameters passed to
#'   generate_multilevel_data_model()
#' @export
generate_multilevel_data <- function(n.bar = 10, J = 30, p = 0.5,
                                     tau.11.star = 0.3, rho2.0W = 0.1, rho2.1W = 0.5, ICC = 0.7,
                                     gamma.00 = 0, gamma.10 = 0.2,
                                     verbose = FALSE, zero.corr = FALSE, ... ) {
    sigma2.W <- 1
    gamma.01 <- sqrt(rho2.0W * ICC / sigma2.W)
    gamma.11 <- sqrt(rho2.1W * tau.11.star)
    tau.00 <- (1 - rho2.0W) * ICC
    if (zero.corr) {
        tau.01 <- 0
    } else {
        tau.01 <- -tau.11.star / 2 - gamma.01 * gamma.11 * sigma2.W
    }
    tau.11 <- (1 - rho2.1W) * tau.11.star
    sigma2.e <- 1 - ICC
    if (verbose) {
        scat( "tau.11* <- %.2f\tICC <- %.2f\trho2.Ws <- %.2f, %.2f\n", tau.11.star, ICC, rho2.0W, rho2.1W)
        scat( "tau.00* <- %.2f\n", gamma.01 ^ 2 * sigma2.W + tau.00)
        scat( "tau.11* <- %.2f\n", gamma.11 ^ 2 * sigma2.W + tau.11)
        scat( "sigma2.e* <- %.2f\n", sigma2.e)
    }
    generate_multilevel_data_model(n.bar = n.bar, J = J, p = p,
                                   gamma.00 = gamma.00, gamma.01 = gamma.01, gamma.10 = gamma.10, gamma.11 = gamma.11,
                                   tau.00 = tau.00, tau.01 = tau.01, tau.11 = tau.11, sigma2.e = sigma2.e,
                                   verbose = verbose, ...)
}



#' @describeIn generate_multilevel_data_model   Simplified version of generate_multilevel_data() with no W covariate.
#'
#' @param control.sd.Y1 Make correlation of random intercept and random slope
#'   such that the variance of the Y1s is 1.0, Default: TRUE
#' @param tau.11.star Total amount of cross site treatment variation
#' @param ICC The ICC, Default: 0.7
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param variable.n Allow n to vary around n.bar, Default: TRUE
#' @param control.sd.Y1 Make correlation of random intercept and random slope
#' @export
generate_multilevel_data_no_cov <- function(n.bar = 10, J = 30, p = 0.5,
                           tau.11.star = 0.3, ICC = 0.7, gamma.00 = 0, gamma.10 = 0.2, verbose = FALSE,
                           variable.n = TRUE, control.sd.Y1 = TRUE, ... ) {
    tau.00 <- ICC
    tau.11 <- tau.11.star
    if (control.sd.Y1) {
        tau.01 <- -tau.11 / 2
    } else {
        tau.01 <- 0
    }
    sigma2.e <- 1 - ICC

    if (verbose) {
        scat( "tau.11* = %.2f\tICC = %.2f\n", tau.11.star, ICC)
        scat( "tau.00* = %.2f\n",  tau.00)
        scat( "tau.11* = %.2f\n",  tau.11)
        scat( "sigma2.e* = %.2f\n", sigma2.e)
    }
    generate_multilevel_data_model(n.bar = n.bar, J = J, p = p, gamma.00 = gamma.00, gamma.10 = gamma.10, gamma.01 = 0, gamma.11 = 0, tau.00 = tau.00, tau.01 = tau.01, tau.11 = tau.11,
                  sigma2.e = sigma2.e, verbose = verbose, variable.n = variable.n, ...)
}


