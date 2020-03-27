#' @title Generate multilevel data from a model
#'
#' @description Given a 2-level model, generate data to specifications
#'
#' @param n.bar average site size, Default: 10
#' @param J number sites, Default: 30
#' @param p prop treated, Default: 0.5
#' @param gamma.00 The mean control outcome, Default: 0
#' @param gamma.10 The ATE, Default: 0.2
#' @param verbose Say stuff while maing data?, Default: FALSE
#' @param gamma.01 Coefficient for W to site control mean
#' @param gamma.11 Coefficient for W to treatment impact
#' @param tau.00 Variance of site conrol means
#' @param tau.01 Correlation of treatment impact and mean site outcome under control
#' @param tau.11 Treatment impact variance
#' @param sigma2.e Residual standard error
#' @param beta.X Coefficient for the individual-level X covariate.  NA means no covariate.
#' @param sigma2.mean.X How much the individual-level X covariate means vary across site.
#' @param variable.n Allow n to vary around n.bar, Default: TRUE
#' @param return.sites Return sites, not individual students, Default: FALSE
#' @param cluster.rand TRUE means cluster-randomized.  FALSE means randomized within site.
#' @param size.impact.correlate    Takes values of -1, 0, or 1: Are site impacts
#'   negatively correlated, uncorrelated, or positively correlated with site
#'   size?
#' @param proptx.impact.correlate  Takes values of -1, 0, or 1: Are proportion
#'   of units treated negatively correlated, uncorrelated, or positively
#'   correlated with site size?
#' @param correlate.strength In [0,1], and describes how correlated the ranking
#'   of site impacts will be with proptx and site size, if they are set to be
#'   correlated.

#' @param variable.p TOADD
#' @param sigma2.W TOADD
#' @param finite.model TOADD
#' @param size.ratio TOADD


#' @return Dataframe of data!
#' @rdname gen_dat_model
#' @export
gen_dat_model <- function(n.bar = 10, J = 30, p = 0.5, gamma.00, gamma.01, gamma.10, gamma.11, tau.00, tau.01, tau.11, sigma2.e, sigma2.W = 1,
  beta.X = NULL, sigma2.mean.X = 0, variable.n = TRUE, variable.p = FALSE, cluster.rand = FALSE, return.sites = FALSE, finite.model = FALSE,
  size.impact.correlate = 0, proptx.impact.correlate = 0, correlate.strength = 0.75, size.ratio = 1 / 3, verbose = FALSE) {
  
  stopifnot(size.impact.correlate %in% c(-1, 0, 1))
  stopifnot(proptx.impact.correlate %in% c(-1, 0, 1))
  if (verbose) {
    scat( "gammas:\t%.2f\t%.2f (%.2f)\n\t%.2f\t%.2f (%.2f)\n", gamma.00, gamma.01, gamma.01 ^ 2, gamma.10, gamma.11, gamma.11 ^ 2)
    scat( "taus:\t%.2f\t%.2f\tsd=%.2f\n\t%.2f\t%.2f\tsd=%.2f\tcor=%.2f\n", tau.00, tau.01, sqrt(tau.00), tau.01, tau.11, sqrt(tau.11), tau.01 / sqrt(tau.00 * tau.11))
    scat( "beta.X:\t%.2f\n", beta.X)
  }

  # generate site sizes (all the same or different sizes)
  if (variable.n) {
    stopifnot(n.bar > 4)
    # nj = rpois( J, n.bar)
    # nj = round( n.bar * runif( J, 0.25, 1.75 ) )
    nj <- 4 + round(block_distn(J, n.bar - 4, size.ratio))
    nj[ nj < 4 ] <- 4
    # if ( any( nj < 4 ) ) {
    #    warning( "Some sites have fewer than 4 units, disallowing 2 tx and 2 co units" )
    #}
  } else {
    nj <- rep(n.bar, J)
  }
  
  if (!is.null( gamma.10 ) || !is.null(gamma.11)) {
    Wj <- rnorm(J, mean = 0, sd = sqrt(sigma2.W))
    include.W <- TRUE
  } else {
    Wj <- 0
    include.W <- FALSE
  }

  if (!is.null(beta.X)) {
    include.X <- TRUE
  } else {
    include.X <- FALSE
  }

  # Generate average control outcome and average ATE for all sites
  Sigma <- matrix(c(tau.00, tau.01, tau.01, tau.11), nrow = 2)
  if (finite.model) {
    # Make a canonical set of site charactaristics and then never change
    # them until J changes.
    if (is.null(CANONICAL) || nrow(CANONICAL) != J) {
      CC <- MASS::mvrnorm(J, c(0, 0), Sigma)
      if (var(CC[, 2]) > 0) {
        CC[, 2] <- sqrt(tau.11) * CC[, 2] / sd(CC[, 2])
      }
      if (var(CC[, 1]) > 0 ) {
        CC[, 1] <- sqrt(tau.00) * CC[, 1] / sd(CC[, 1])
      }
      CANONICAL <<- CC
    }
    mv <- CANONICAL
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

  if (return.sites) {
    if (cluster.rand) {
      df <- data.frame(n = nj, W = Wj, beta.0 = beta.0j, beta.1 = beta.1j, u0 = mv[,1], u1 = mv[,2], Zj = Zj)
    } else {
      df <- data.frame(n = nj, W = Wj, beta.0 = beta.0j, beta.1 = beta.1j, u0 = mv[,1], u1 = mv[,2])
    }
    if (!include.W) {
      df$W <- NULL
    }
    return( df )
  } else {
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
      dd <- dd %>% group_by( sid ) %>% mutate(Z = as.numeric(sample(n()) <= n() * p))
    }
    Zij <- dd$Z
    dd$p <- NULL
    
    # individual covariates and residuals
    if (include.X) {
      stopifnot(sigma2.mean.X < 1)
      Xbar <- rnorm(J, mean = 0, sd = sqrt(sigma2.mean.X))
      X <- rnorm(N, mean = Xbar[sid], sd = sqrt(1 - sigma2.mean.X))
      stopifnot(sigma2.e - beta.X ^ 2 > 0)
      e <- rnorm(N, mean = 0, sd = sqrt(sigma2.e - beta.X ^ 2))
      Y0 <- beta.0j[sid] + beta.X * X + e
    } else {
      e <- rnorm(N, mean = 0, sd = sqrt(sigma2.e))
      Y0 <- beta.0j[sid] + e
    }
    Y1 <- Y0 + beta.1j[sid]
    df <- data.frame(sid = as.factor(sid), Y0 = Y0, Y1 = Y1, Z = Zij, Yobs = ifelse(Zij, Y1, Y0))
    if (include.X) {
      df$X <- X
    }
    if (include.W) {
      df$W <- Wj[df$sid]
    }
    attr(df, "tau.S") <- sd( beta.1j)
    df
  }
}
