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
make_data_linear = function(X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65), a = 0, b = 0, ATE = 0.2, d = 0) {
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


#' Given potential outcomes schedule, randomize within block and generate
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
    dat = make_data( c( 2, 5, 10 ) )
    debug( add_obs_data )
    add_obs_data( dat, p=0.2 )
}

#' Function to generate individual-level data from a list of block sizes and
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
generate_individuals_from_blocks <- function(n_k, alpha = 0, beta = 0, sigma_c = 1, sigma_t = 1, corr = 1, exact = FALSE) {
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
#' @param sigma_alpha Standard deviation of the separation of the block mean Y0
#' @param sigma_tau Standard deviation of the separation of the block mean
#'   treatment effects (Y1-Y0)
#' @param tau Cross site VARIANCE of site-level ATEs.
#' @param sigma_0 Standard deviation of residual Y0 added to block means (can be
#'   vector for individual variances per block)
#' @param sigma_1 As `sigma_0` but for Y1s.
#' @param corr Correlation of Y0, Y1 within a block (can be vector of length K
#'   for different blocks)
#' @param exact Passed to mvrnorm to control how block means are generated.
#'
#' @return Dataframe with block indicators, Y0, and Y1.
#'
#' @export
make_data <- function(n_k, sigma_alpha = 1, sigma_tau = 0, tau = 5, sigma_0 = 1, sigma_1 = 1, corr = 0.5, exact = FALSE) {
  K <- length(n_k)
  percents <- seq(from = (1 - 1 / (K + 1)), to = 1 / (K + 1), by = -1 /(K + 1))
  alpha <- qnorm(percents, 0, sigma_alpha)
  beta <- qnorm(percents, tau, sigma_tau)
  generate_individuals_from_blocks(n_k, alpha, beta = beta, sigma_c = sigma_0, sigma_t = sigma_1, corr = corr, exact = exact)
}

#' Make a random simulated dataset from a linear model.
#'
#' Generate data, form_blocks_from_continuous, and randomize within block and generate
#' observed potential outcomes
#'
#' @rdname make_data_linear
#' @param method How to block
#' @param X vector of indiviual level covariates
#' @param a Intercept of Y0
#' @param b Main effect of X
#' @param ATE Average Tx
#' @param d Interaction effect term.
#' @param p Proportion of units treated (as close as possible given block sizes)
#' @return Dataframe with original potential outcomes and observed outcome based on random assigment.
#' @export
make_obs_data_linear <- function(X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65), p = 0.5, a = 0, b = 0, ATE = 0.2, d = 0,
  method = c("small", "pair", "big", "none")) {
  dat <- make_data_linear(X, a, b, ATE, d)
  dat$B <- form_blocks_from_continuous(dat$X, method = method)
  dat <- add_obs_data(dat, p = p)
  dat
}

#' Make a random simulated dataset from a list of block sizes
#'
#' Generate data, form_blocks_from_continuous, and randomize within block and generate observed
#' potential outcomes
#'
#' @rdname make_data
#' @param n_k List of block sizes
#' @param p Proportion of units treated (as close as possible given block
#'   sizes).  This can be a vector with a probability for each block.
#' @param ... Parameters to be passed to make_data()
#' @return Dataframe with original potential outcomes and observed outcome based
#'   on random assigment.
#' @export
make_obs_data = function(n_k = c(2, 3, 4, 8), p = 0.5, ... ) {
  dat <- make_data( n_k = n_k, ... )
  dat <- add_obs_data(dat, p = p)
  dat
}

#' Table of data for simulations
#'
#' Function that returns a summary of block level true values for sims.
#'
#' @param Y vector of all outcomes
#' @param Z vector that indicates if outcome is under treatment or control
#' @param B block ids
#' @param data alternatively is matrix of Y,Z,B
#' @param p.mat  matrix with first column B,second column prop treated in that block, p
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
block_data_sim <- function(Y, Z, B, p.mat, data = NULL) {
  if (!is.null(data)) {
    Y <- data[, 1]
    Z <- data[, 2]
    B <- data[, 3]
  }
  n <- length(Y)
  #Quick test that input is correct
  if (is.numeric(Z) == FALSE) {
    return("Treatment indicator should be vector of ones and zeros")
  }
  if ((sum(Z == 1) + sum(Z == 0)) != n) {
    return("Treatment indicator should be vector of ones and zeros")
  }
  #First convert block ids into numbers
  B <- factor(B)
  B <- as.numeric(B)
  p.mat$B <- factor(p.mat$B)
  p.mat$B <- as.numeric(p.mat$B)
  #Get number of units assigned to each treatment
  #In each block
  n_matrix <- aggregate(list(n_k = B), list(B = B), FUN = length)
  n_matrix$n_k <- n_matrix$n_k / 2
  n_ctk_matrix <- merge(n_matrix, p.mat, by = "B")
  n_ctk_matrix$n1 <- n_ctk_matrix$n_k * n_ctk_matrix$p
  n_ctk_matrix$n0 <- n_ctk_matrix$n_k - n_ctk_matrix$n1
  treated_mat <- cbind(Y[Z == 1], B[Z == 1])
  control_mat <- cbind(Y[Z == 0], B[Z == 0])
  Y1_matrix <- aggregate(list(Ybar1 = treated_mat[, 1]), list(B = treated_mat[, 2]), FUN = mean)
  Y0_matrix <- aggregate(list(Ybar0 = control_mat[, 1]), list(B = control_mat[, 2]), FUN = mean)
  Ybar_matrix <- merge(Y1_matrix, Y0_matrix, by = "B")
  var1_matrix <- aggregate(list(var1 = treated_mat[, 1]), list(B = treated_mat[, 2]), FUN = var)
  var0_matrix <- aggregate(list(var0 = control_mat[, 1]), list(B = control_mat[, 2]), FUN = var)
  var_matrix <- merge(var1_matrix, var0_matrix, by = "B")
  overall_mat <- merge(n_ctk_matrix, Ybar_matrix, by = "B")
  overall_mat <- merge(overall_mat, var_matrix, by = "B")
  overall_mat$se_ney <- sqrt(overall_mat$var1 / overall_mat$n1 + overall_mat$var0 / overall_mat$n0)
  drops <- c("n_k", "p")
  overall_mat <- overall_mat[ , !(names(overall_mat) %in% drops)]
  return(overall_mat)
}


#' Calculates variance of treatment effects.
#'
#' Function that helps calculate bias by caculating the true variance of treatment effects.
#' @param tau_vec  vector of treatment effects
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
s_tc_func <- function(tau_vec) {
  s.tc<-var(tau_vec)
  return(s.tc)
}

#' Rescaled variance
#'
#' This is the residual sum of squares
#'
#' Function that calculates scaled variance to help with bias calculation.
#' @param Y  vector of potential outcomes
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
within_blk_var <- function(Y) {
  sum((Y - mean(Y)) ^ 2)
}

