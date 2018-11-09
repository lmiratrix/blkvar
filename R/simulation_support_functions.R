##
## Functions to make fake data for testing and simulations
##


#' Summarise simulation data by block (where we know both Y0 and Y1)
#'
#' @param dat Dataframe with defined Y0, Y1, and blk variables.
#'
#' @return dataframe with summary statistics by block
#' @export
calc.summary.stats.oracle = function( data, Y0="Y0", Y1="Y1", Z="Z", blk="blk" ) {
    require( tidyverse )

            data = rename( data, Y0 = !!rlang::sym(Y0),
                       Y1 = !!rlang::sym(Y1),
                       Z = !!rlang::sym(Z),
                       blk = !!rlang::sym(blk) )

            sdat <- data %>% dplyr::group_by( blk ) %>%
        dplyr::summarise( n = n(),
                          mu0 = mean( Y0 ),
                          mu1 = mean( Y1 ),
                          tau = mu1 - mu0,
                          sd0 = sd( Y0 ),
                          sd1 = sd( Y1 ),
                          corr = cor( Y0, Y1 ) )


    as.data.frame( sdat )
}


#' Make blocks out of a continuous covariate by cutting it into pieces that are
#' ideally relatively homogenous.
#'
#' Given a dataset, create a bunch of blocks based on passed covariate and
#' return the factor of blocks
#' @param X covariate vector to block on
#' @param method How to block.
#' @param num.blocks If method is small, how many blocks to attempt to make
#'   (optional argument)
#' @return Vector with one element per element of `X`
#' @export
make.blocks = function(X,
                       method = c("small", "pair", "big", "none"),
                       num.blocks) {
    method = match.arg(method)
    X.orig = X
    X.order = rank(X, ties.method = "first")
    X = sort(X)
    if (method == "small") {
        dels = diff(X)
        dels = jitter(dels)
        N = length(X)
        if (!missing(num.blocks)) {
            ct = sort(dels, decreasing = TRUE)[num.blocks - 1]
        } else {
            ct = quantile(dels, 2 / 3 + 1 / N)
        }
        cuts = (dels >= ct) &
            (c(dels[-1], Inf) < ct) & (c(Inf, dels[-(N - 1)]) < ct)
        cuts = 1 + cumsum(cuts)
        B = paste("B", c(1, cuts), sep = "")
    } else if (method == "pair") {
        B = paste("B", rep(1:(length(X) / 2), each = 2), sep = "")
    } else if (method == "big") {
        minblk = trunc(length(X) / 4)
        B = cut(X, quantile(X, seq(0, 1, length.out = minblk + 1)), include.lowest =
                    TRUE)
    } else {
        B = paste("B", rep(1, length(X)), sep = "")
    }
    B[X.order]
}




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
make.data.linear = function(X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65),
                     a = 0,
                     b = 0,
                     ATE = 0.2,
                     d = 0) {
    X.sd = round((X - mean(X)) / sd(X), digits = 1)


    # Quadratic relationship, OLS no help
    #Y0 = a + b*X^2 + rnorm( length(X), 0, 1 )

    # Linear relationship, OLS help!
    Y0 = a + b * X.sd + rnorm(length(X), 0, 1)

    Y1 = Y0 + ATE + d * X.sd

    Y1 = Y1 / sd(Y0)
    Y0 = Y0 / sd(Y0)


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
add.obs.data = function(dat,
                        p = 0.5,
                        Y0 = "Y0",
                        Y1 = "Y1",
                        blockvar = "blk") {
    N = nrow(dat)
    blk = as.factor( dat[[blockvar]] )
    K = nlevels( blk )

    if ( length( p ) == 1 ) {
        p = rep( p, K )
    }

    # Make initial treatment assignment vector to shuffle in subsequent step
    Z = rep( NA, N )
    for ( i in 1:nlevels(blk) ) {
        nk = sum( blk == levels(blk)[[i]] )
        stopifnot( nk > 1 )
        ntx = round( nk * p[[i]] )
        if ( ntx == 0 ) {
            ntx = 1
        }
        if ( ntx == nk ) {
            ntx = nk - 1
        }
        Z[ blk == levels(blk)[[i]] ] = sample( nk ) <= ntx
    }
    #    dat$Z = randomizationInference::blockRand(Z, 1, dat[[blockvar]])[[1]]

    dat$Z = as.numeric( Z )

    dat$Yobs = ifelse(Z, dat[[Y1]], dat[[Y0]])

    dat
}


if ( FALSE ) {
    dat = make.data( c( 2, 5, 10 ) )
    debug( add.obs.data )
    add.obs.data( dat, p=0.2 )
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
generate.individuals.from.blocks <- function( n_k, alpha = 0, beta = 0, sigma_c = 1, sigma_t = 1, corr = 1, exact=FALSE){
    library(MASS)
    K = length( n_k )
    if ( length( alpha ) == 1 ) {
        alpha = rep( alpha, K )
    } else {
        stopifnot( length( alpha ) == K )
    }
    if ( length( beta ) == 1 ) {
        beta = rep( beta, K )
    } else {
        stopifnot( length( beta ) == K )
    }
    if ( length( sigma_c ) == 1 ) {
        sigma_c = rep( sigma_c, K )
    } else {
        stopifnot( length( sigma_c ) == K )
    }
    if ( length( sigma_t ) == 1 ) {
        sigma_t = rep( sigma_t, K )
    } else {
        stopifnot( length( sigma_t ) == K )
    }
    if ( length( corr ) == 1 ) {
        corr = rep( corr, K )
    } else {
        stopifnot( length( corr ) == K )
    }

    Y<-matrix(nrow=sum(n_k), ncol=2)
    j<-1
    for(i in 1:K){
        mu<-c(alpha[i], alpha[i]+beta[i])
        Sigma<-matrix(c(sigma_c[i], corr[i]*sqrt(sigma_c[i])*sqrt(sigma_t[i]), corr[i]*sqrt(sigma_c[i])*sqrt(sigma_t[i]), sigma_t[i]), ncol=2)
        Y[j:sum(n_k[1:i]),]<-mvrnorm(n_k[i], mu, Sigma, empirical=exact)
        j<-j+n_k[i]
    }
    blk =  rep( 1:K, n_k )
    blk = factor( blk, levels = 1:K, labels=paste( "B", 1:K, sep="" ) )
    data.frame( blk = blk, Y0 = Y[,1], Y1=Y[,2] )
}

if ( FALSE ) {
    dt = generate.individuals.from.blocks( c( 4, 8 ), c( 0, 10 ), c(1, 3), c( 10, 1 ), c( 3, 1 ), c( 0, 1 ), TRUE )
    dt
    levels( dt$blk )
    block.data( dt )
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
make.data = function( n_k, sigma_alpha = 1, sigma_tau = 0, tau = 5,
                      sigma_0 = 1, sigma_1 = 1, corr = 0.5, exact=FALSE ) {
    K = length( n_k )
    percents<-seq(from=(1-1/(K+1)), to=1/(K+1), by=-1/(K+1))
    alpha <-qnorm(percents, 0, sigma_alpha)
    beta <-qnorm(percents, tau, sigma_tau )
    generate.individuals.from.blocks( n_k, alpha, beta=beta, sigma_c = sigma_0, sigma_t=sigma_1, corr=corr, exact=exact )
}



#' Make a random simulated dataset from a linear model.
#'
#' Generate data, make blocks, and randomize within block and generate
#' observed potential outcomes
#'
#' @rdname make.data.linear
#' @param method How to block
#' @param p Proportion of units treated (as close as possible given block sizes)
#' @return Dataframe with original potential outcomes and observed outcome based on random assigment.
#' @export
make.obs.data.linear = function(X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65),
                        p = 0.5,
                         a = 0,
                         b = 0,
                         ATE = 0.2,
                         d = 0,
                         method = c("small", "pair", "big", "none")) {
    dat = make.data.linear(X, a, b, ATE, d)
    dat$blk = make.blocks(dat$X, method = method)
    dat = add.obs.data(dat, p = p)
    dat
}



#' Make a random simulated dataset from a list of block sizes
#'
#' Generate data, make blocks, and randomize within block and generate observed
#' potential outcomes
#'
#' @rdname make.data
#' @param n_k List of block sizes
#' @param p Proportion of units treated (as close as possible given block
#'   sizes).  This can be a vector with a probability for each block.
#' @param ... Parameters to be passed to make.data()
#' @return Dataframe with original potential outcomes and observed outcome based
#'   on random assigment.
#' @export
make.obs.data = function( n_k = c( 2, 3, 4, 8 ), p = 0.5, ... ) {
    dat = make.data( n_k = n_k, ... )
    dat = add.obs.data(dat, p = p)
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
#' @param p.mat  matrix with first column Blk_ID,second column prop treated in that block, p
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
block.data.sim<-function(Y, Z, B, p.mat, data=NULL){
  if(!is.null(data)){
    Y<-data[,1]
    Z<-data[,2]
    B<-data[,3]
  }
  n<-length(Y)
  #Quick test that input is correct
  if(is.numeric(Z)==F){
    return("Treatment indicator should be vector of ones and zeros")
  }
  if((sum(Z==1)+sum(Z==0))!=n){
    return("Treatment indicator should be vector of ones and zeros")
  }
  #First convert block ids into numbers
  blk_id<-factor(B)
  blk_id<-as.numeric(blk_id)
  p.mat$Blk_ID<-factor(p.mat$Blk_ID)
  p.mat$Blk_ID<-as.numeric(p.mat$Blk_ID)
  #Get number of units assigned to each treatment
  #In each block
  n_matrix<-aggregate(list(n_k=blk_id), list(Blk_ID=blk_id), FUN=length)
  n_matrix$n_k<-n_matrix$n_k/2
  n_ctk_matrix<-merge(n_matrix, p.mat, by="Blk_ID")
  n_ctk_matrix$Num_Trt<-n_ctk_matrix$n_k*n_ctk_matrix$p
  n_ctk_matrix$Num_Ctrl<-n_ctk_matrix$n_k-n_ctk_matrix$Num_Trt
  treated_mat<-cbind(Y[Z==1], blk_id[Z==1])
  control_mat<-cbind(Y[Z==0], blk_id[Z==0])
  Y1_matrix<-aggregate(list(Y1=treated_mat[,1]), list(Blk_ID=treated_mat[,2]), FUN=mean)
  Y0_matrix<-aggregate(list(Y0=control_mat[,1]), list(Blk_ID=control_mat[,2]), FUN=mean)
  Ybar_matrix<-merge(Y1_matrix, Y0_matrix, by="Blk_ID")
  var1_matrix<-aggregate(list(var1=treated_mat[,1]), list(Blk_ID=treated_mat[,2]), FUN=var)
  var0_matrix<-aggregate(list(var0=control_mat[,1]), list(Blk_ID=control_mat[,2]), FUN=var)
  var_matrix<-merge(var1_matrix, var0_matrix, by="Blk_ID")
  overall_mat<-merge(n_ctk_matrix, Ybar_matrix, by="Blk_ID")
  overall_mat<-merge(overall_mat, var_matrix, by="Blk_ID")
  overall_mat$se_ney<-sqrt(overall_mat$var1/overall_mat$Num_Trt + overall_mat$var0/overall_mat$Num_Ctrl)
  drops <- c("n_k","p")
  overall_mat<-overall_mat[ , !(names(overall_mat) %in% drops)]
  return(overall_mat)
}


#' Function that calculates variance of treatment effects.
#'
#' Function that helps calculate bias by caculating the tue variance of treatment effects.
#' @param tau_vec  vector of treatment effects
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
s.tc.func<-function(tau_vec){
  s.tc<-var(tau_vec)
  return(s.tc)
}

#' Variance function multiplied by n-1.
#'
#' Function that calculates scaled variance to help with bias calculation.
#' @param Y  vector of potential outcomes
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
within_blk_var<-function(Y){
  sum((Y-mean(Y))^2)
}

