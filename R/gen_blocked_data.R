##
## Functions to make fake data for testing and simulations
##

#' Make data from a linear model specification
#'
#' Function that, given a covariate vector X, returns a dataframe of potential
#' outcomes and blocks.
#'
#' This uses the model: Y = a + bX + ATE Z + d X ATE + epsilon
#' (It will standardize X for this model, and then standardize Y0, Y1 so sd(Y0) = 1.)
#'
#' @param X vector of indiviual level covariates
#' @param a Intercept of Y0
#' @param b Main effect of X
#' @param ATE Average Tx
#' @param d Interaction effect term.
#' @export
generate_blocked_data_linear = function(X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65),
                            a = 0, b = 0, ATE = 0.2, d = 0) {
  X.sd <- round((X - mean(X)) / sd(X), digits = 1)
  # Quadratic relationship, OLS no help
  # Y0 = a + b*X^2 + rnorm( length(X), 0, 1 )

  # Linear relationship, OLS help!
  Y0 <- a + b * X.sd + rnorm(length(X), 0, 1)
  Y1 <- Y0 + ATE + d * X.sd
  Y1 <- Y1 / sd(Y0)
  Y0 <- Y0 / sd(Y0)
  data.frame(Y0 = Y0, Y1 = Y1, X = X)
}


#' Randomize and calculate observed outcome
#'
#' Given a passed potential outcomes schedule, randomize within block and generate
#' observed potential outcomes and add them to the schedule.
#'
#' @param dat Dataframe with a named Y0, Y1, and block column
#' @param p Proportion of units treated.  Can be vector and will treat that
#'   proportion in each block, rounded to nearest and with at least 1 treatment
#'   and control unit.
#' @param Y0 name of Y0 column
#' @param Y1 name of Y1 column
#' @param blockvar name of blocking column.  This column will be converted to a
#'   factor, if it is not already, and the order of p corresponds to the levels
#'   of this factor.
#'
#' @return augmented `dat` with Z and Yobs columns
#' @export
add_obs_data <- function(dat, p = 0.5, Y0 = "Y0", Y1 = "Y1", blockvar = "B") {
  N <- nrow(dat)
  B <- as.factor(dat[[blockvar]])
  K <- nlevels(B)

  if (length(p) == 1) {
    p <- rep(p, K)
  }

  # Make initial treatment assignment vector to shuffle in subsequent step
  Z <- rep(NA, N)
  for (i in 1:nlevels(B)) {
    nk <- sum( B == levels(B)[[i]])
    stopifnot(nk > 1)
    ntx <- round(nk * p[[i]])
    if (ntx == 0) {
      ntx <- 1
    }
    if (ntx == nk) {
      ntx <- nk - 1
      }
    Z[ B == levels(B)[[i]] ] <- sample( nk ) <= ntx
  }
  #    dat$Z = randomizationInference::blockRand(Z, 1, dat[[blockvar]])[[1]]
  dat$Z <- as.numeric(Z)
  dat$Yobs <- ifelse(Z, dat[[Y1]], dat[[Y0]])
  dat
}

if ( FALSE ) {
    dat = generate_blocked_data( c( 2, 5, 10 ) )
    debug( add_obs_data )
    add_obs_data( dat, p=0.2 )
}

#' Generate individual-level data from a list of block sizes and
#' block characteristics.
#'
#' This generates potential outcomes by sampling from the specified bivariate
#' normal distributions within each block.
#'
#' @param n_k Vector of block sizes
#' @param alpha Vector of (expected) means of the control potential outcome
#' @param beta Vector of the block-level Average Treatment Effects
#' @param sigma_c Standard deviation of the control potential outcomes.
#' @param sigma_t Standard deviation of the treatment potential outcomes
#' @param corr Correlation of the potential outcomes.
#' @param exact If TRUE generate data so we match the desired means and moments
#'   exactly.  False means pull from bivariate normal.
#' @return Matrix of the potential outcomes and block ids.
#' @export
generate_individuals_from_blocks <- function(n_k, alpha = 0, beta = 0,
                                             sigma_c = 1, sigma_t = 1, corr = 1, exact = FALSE) {
  K <- length(n_k)
  if (length(alpha) == 1) {
    alpha <- rep(alpha, K)
  } else {
    stopifnot(length(alpha) == K)
  }
  if (length(beta) == 1) {
    beta <- rep(beta, K)
  } else {
    stopifnot(length(beta) == K)
  }
  if (length(sigma_c) == 1) {
    sigma_c <- rep(sigma_c, K)
  } else {
    stopifnot(length(sigma_c) == K)
  }
  if (length(sigma_t) == 1) {
    sigma_t <- rep(sigma_t, K)
  } else {
    stopifnot(length(sigma_t) == K)
  }
  if (length(corr) == 1) {
    corr <- rep(corr, K)
  } else {
    stopifnot(length(corr) == K)
  }
  Y <- matrix(nrow = sum(n_k), ncol = 2)
  j <- 1
  for (i in 1:K) {
    mu <- c(alpha[i], alpha[i] + beta[i])
    Sigma <- matrix(c(sigma_c[i], corr[i] * sqrt(sigma_c[i]) * sqrt(sigma_t[i]), corr[i] * sqrt(sigma_c[i]) * sqrt(sigma_t[i]), sigma_t[i]), ncol = 2)
    Y[j:sum(n_k[1:i]), ] <- MASS::mvrnorm(n_k[i], mu, Sigma, empirical = exact)
    j <- j + n_k[i]
  }
  B <- rep(1:K, n_k)
  B <- factor(B, levels = 1:K, labels = paste("B", 1:K, sep="" ))
  data.frame(B = B, Y0 = Y[, 1], Y1 = Y[, 2])
}

if ( FALSE ) {
    dt = generate_individuals_from_blocks( c( 4, 8 ), c( 0, 10 ), c(1, 3), c( 10, 1 ), c( 3, 1 ), c( 0, 1 ), TRUE )
    dt
    levels( dt$B )
    calc_summary_stats( data=dt )
}


#' Make simulated dataset from a list of block sizes
#'
#' This method is the one that generates the simulation data used in Pashley &
#' Miratrix.
#'
#' The block means are sampled from a multivariate normal distribution.  This
#' can be controlled so the variances are exact using the `exact` flag.
#'
#' @param n_k Vector of block sizes. Let K be length of this vector.
#' @param sigma_alpha Standard deviation of the block mean Y0s.
#' @param sigma_beta Standard deviation of the block mean treatment effects
#'   (Y1-Y0)s.
#' @param beta Block Average ATE.
#' @param sigma_0 Standard deviation of residual Y0 added to block means (can be
#'   vector for individual variances per block).
#' @param sigma_1 As `sigma_0` but for Y1s.
#' @param corr Correlation of Y0, Y1 within a block (can be vector of length K
#'   for different blocks).
#' @param exact Passed to mvrnorm to control how block means are generated.
#'
#' @return Dataframe with block indicators, Y0, and Y1.
#'
#' @export
generate_blocked_data <- function(n_k, sigma_alpha = 1, sigma_beta = 0,
                                  beta = 5, sigma_0 = 1,
                                  sigma_1 = 1, corr = 0.5, exact = FALSE) {
  K <- length(n_k)
  percents <- seq(from = (1 - 1 / (K + 1)), to = 1 / (K + 1), by = -1 /(K + 1))
  alpha <- qnorm(percents, 0, sigma_alpha)
  beta <- qnorm(percents, beta, sigma_beta)
  generate_individuals_from_blocks(n_k = n_k,
                                   alpha = alpha, beta = beta,
                                   sigma_c = sigma_0, sigma_t = sigma_1, corr = corr,
                                   exact = exact)
}



#' Make a random simulated dataset from a linear model.
#'
#' Generate data, form_blocks_from_continuous, and randomize within block and generate
#' observed potential outcomes
#'
#' @rdname generate_blocked_data_linear
#' @param method How to block
#' @param X vector of indiviual level covariates
#' @param a Intercept of Y0
#' @param b Main effect of X
#' @param ATE Average Tx
#' @param d Interaction effect term.
#' @param p Proportion of units treated (as close as possible given block sizes)
#' @return Dataframe with original potential outcomes and observed outcome based on random assigment.
#' @export
generate_blocked_data_obs_linear <- function(X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65), p = 0.5, a = 0, b = 0, ATE = 0.2, d = 0,
  method = c("small", "pair", "big", "none")) {
  dat <- generate_blocked_data_linear(X, a, b, ATE, d)
  dat$B <- form_blocks_from_continuous(dat$X, method = method)
  dat <- add_obs_data(dat, p = p)
  dat
}

#' Make a random simulated dataset from a list of block sizes
#'
#' Generate data, form_blocks_from_continuous, and randomize within block and generate observed
#' potential outcomes
#'
#' @rdname generate_blocked_data
#' @param n_k List of block sizes
#' @param p Proportion of units treated (as close as possible given block
#'   sizes).  This can be a vector with a probability for each block.
#' @param ... Parameters to be passed to generate_blocked_data()
#' @return Dataframe with original potential outcomes and observed outcome based
#'   on random assigment.
#' @export
generate_blocked_data_obs = function(n_k = c(2, 3, 4, 8), p = 0.5, ... ) {
  dat <- generate_blocked_data( n_k = n_k, ... )
  dat <- add_obs_data(dat, p = p)
  dat
}




#' Calculate residual sum of squares
#'
#' Function that calculates scaled variance to help with bias calculation.
#'
#' @param Y  vector of potential outcomes
#' @importFrom stats aggregate lm quantile rnorm sd var
#'
#' @noRd
within_blk_var <- function(Y) {
  sum((Y - mean(Y)) ^ 2)
}




#' Cut continuous covariate into similar-valued blocks
#'
#' Make blocks out of a continuous covariate by cutting it into pieces that are
#' ideally relatively homogenous. This will create a bunch of blocks based on
#' passed covariate and return a factor vector of those blocks.
#'
#' @param X covariate vector to block on
#' @param method How to block.  "small"  "pair" makes matched pairs. "big" makes as many
#'   blocks of at least 4 as possible, trying to keep blocks small. "none" makes
#'   a single block.
#' @param num.blocks If method is small, how many blocks to attempt to make
#' @return Vector with one element per element of `X`
#' @export
#' @examples
#' table( form_blocks_from_continuous( 1:22, method="small" ) )
#' table( form_blocks_from_continuous( 1:22, method="big", num.blocks=3 ) )
#' table( form_blocks_from_continuous( 1:22, method="pair" ) )
#' table( form_blocks_from_continuous( 1:22, method="none" ) )
#' table( form_blocks_from_continuous( 1:22, method="big" ) )
form_blocks_from_continuous <- function(X, method = c("small", "pair", "big", "none"), num.blocks) {
    method <- match.arg(method)
    X.orig = X
    X.order <- rank(X, ties.method = "first")
    X <- sort(X)
    if (method == "small") {
        dels <- diff(X)
        dels <- jitter(dels)
        N <- length(X)
        if (!missing(num.blocks)) {
            ct <- sort(dels, decreasing = TRUE)[num.blocks - 1]
        } else {
            ct <- quantile(dels, 2 / 3 + 1 / N)
        }
        cuts <- (dels >= ct) & (c(dels[-1], Inf) < ct) & (c(Inf, dels[ - (N - 1)]) < ct)
        cuts <- 1 + cumsum(cuts)
        B <- paste("B", c(1, cuts), sep = "")
    } else if (method == "pair") {
        B <- paste("B", rep(1:(length(X) / 2), each = 2), sep = "")
    } else if (method == "big") {
        minB <- trunc(length(X) / 4)
        B <- cut(X, quantile(X, seq(0, 1, length.out = minB + 1)), include.lowest = TRUE)
    } else {  # "none"
        B <- paste("B", rep(1, length(X)), sep = "")
    }
    B[X.order]
}

