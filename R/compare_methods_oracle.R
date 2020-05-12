
#' Describe an oracle (simulated) dataset.
#'
#' This utility function gives various summary statistics (at the block level)
#' for a dataset with all potential outcomes given.
#'
#' @param data Dataframe with defined Y0, Y1 (the potential outcomes), Z (sample
#'   treatment assignment vector), and sid (a categorical blocking variable).
#' @param Y0 name of Y0 column
#' @param Y1 name of Y1 column
#' @param Z example vector of treatment assignment that allows for calculation
#'   of proportion treated and so forth.
#' @param sid Name of the block ID column in the dataframe.
#' @return (invisible) the summary stats by block.  Side effect of printing to
#'   screen the description of the data.
#' @export

describe_data <- function(data, Y0 = "Y0", Y1 = "Y1", Z = "Z", sid = "sid") {
    # old param list: Y0, Y1, Z, sid, data = NULL
    data <- dplyr::rename(data,
                          Y0 = !!rlang::sym(Y0),
                          Y1 = !!rlang::sym(Y1),
                          Z = !!rlang::sym(Z),
                          sid = !!rlang::sym(sid))
    sites <- data %>% group_by(!!rlang::sym(sid)) %>%
        dplyr::summarise(Y0.bar = mean(Y0),
                         Y1.bar = mean(Y1),
                         beta = Y1.bar - Y0.bar,
                         n = n(),
                         p.Tx = mean(Z))

    site.var = sd(sites$beta)
    scat( "\n\nDescription of site distribution:" )
    scat( "\n\tsite average variation:\tsd(Y0.bar) = %.3f\tsd(Y1.bar) = %.3f",
          sd(sites$Y0.bar), sd(sites$Y1.bar))
    scat( "\n\tmarginal variation:\tsd(Y0) = %.3f\t\tsd(Y1) = %.3f\t\tcor(Y0,Y1) = %.3f",
          sd(data$Y0), sd(data$Y1), cor(data$Y0, data$Y1))
    scat( "\n\tsd( beta ) = %.3f\n\tcorr( Y0.bar, Y1.bar ) = %.3f\tcov = %.3f",
          site.var, cor(sites$Y0.bar, sites$Y1.bar), cov( sites$Y0.bar, sites$Y1.bar))
    scat( "\n\tcor( Y0.bar, beta ) = %.3f", cor(sites$Y0.bar, sites$beta))
    scat( "\n\tmean( n ) = %.3f\tsd( n ) = %.3f", mean(sites$n), sd( sites$n))
    scat( "\n\tmean( Z.bar ) = %.3f\t(%.3f)\n", mean(sites$p.Tx), sd(sites$p.Tx))
    invisible(sites)
}





#' Make table of data for simulations
#'
#' Function that returns a summary of block level true values for simulations
#' (so works with full schedule of potential outcomes).
#'
#'
#' @inheritParams compare_methods_oracle
#' @param Z Instead of passing pre-computed p_mat, one can pass an example
#'   treatment assignment vector Z of 0s and 1s. In this case method will
#'   tabulate to generate p_mat.  p_mat or Z must be null to avoid ambiguity.
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
calc_summary_stats_oracle <- function(Y0, Y1, B, data = NULL, p_mat = NULL, Z = NULL) {

    if (!is.null(data)) {
        Y0 <- eval(substitute(Y0), data)
        Y1 <- eval(substitute(Y1), data)
        B <- eval(substitute(B), data)
        Z <- eval(substitute(Z), data)
    }
    B <- as.character( B )

    stopifnot( is.null( p_mat ) || is.null(Z) )

    if ( !is.null( Z ) ) {
        stopifnot( sum( Z == 0 ) + sum( Z == 1 ) == length( Y0 ) )
        p_mat = aggregate( list( p=Z ), list( B=B ), mean )
        p_mat$B = as.character(p_mat$B)
    }

    if ( !is.null( p_mat ) ) {
        stopifnot( nrow( p_mat ) == length( unique( B ) ) )
        stopifnot( ncol( p_mat ) == 2 )
        names( p_mat ) = c("B","p" )
    }
    stopifnot( length(Y0) == length(Y1) )
    stopifnot( length(Y0) == length(B) )

    # Use aggregation to simplify and keep the B name right
    data <- data.frame( Y0 = Y0,
                        Y1 = Y1,
                        B = B )

    sdat <- data %>% dplyr::group_by( B ) %>%
        dplyr::summarise(n = n(),
                         Ybar1 = mean(Y1),
                         Ybar0 = mean(Y0),
                         tau = Ybar1 - Ybar0,
                         var1 = var(Y1),
                         var0 = var(Y0),
                         corr = cor(Y0, Y1) )

    sdat = merge( sdat, p_mat, by="B", all=TRUE )

    sdat = dplyr::mutate(sdat,
                         n0 = round( n * (1-p) ),
                         n1 = round( n * p ),
                         se_ney = sqrt( var0 / n0 + var1 / n1 ) )
    stopifnot( all( sdat$n0 + sdat$n1 == sdat$n ) )
    sdat$B = as.character( sdat$B )

    return( sdat )
}



#' Compare estimators using the true values (i.e., full schedule of potential
#' outcomes)
#'
#'
#' This method is for simulaton studies where we know all the potential
#' outcomes.
#'
#' Function that returns some variance function estimates based on all potential
#' outcomes.
#'
#' @param Y0 Vector of control potential outcomes  (or name of column in data
#'   holding same).
#' @param Y1 Vector of treatment potential outcomes  (or name of column in data
#'   holding same).
#' @param B block ids  (or name of column in data holding same).
#' @param data alternatively is matrix of Y,Z,B
#' @param p_mat  Data.frame with exactly two columns, and one row per block.
#'   First column is the categorical covariate denoting blocks, with same levels
#'   as in B, above. Second column is proportion treated in that block.
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @importFrom msm deltamethod
#' @export
compare_methods_oracle <- function(Y0, Y1, B, data = NULL, p_mat = NULL ) {

    if (!is.null(data)) {
        Y0 <- eval(substitute(Y0), data)
        Y1 <- eval(substitute(Y1), data)
        B <- eval(substitute(B), data)
    }

    stopifnot( length(Y0) == length(Y1) )
    stopifnot( length(Y0) == length(B) )

    data_table <- calc_summary_stats_oracle(Y0, Y1, B, p_mat = p_mat)

    s.tc.bk <- aggregate(list(s.tc.bk = Y1 - Y0), list(B = B), FUN = var)

    n <- sum(data_table$n1) + sum(data_table$n0)
    data_table$nk <- data_table$n1 + data_table$n0

    Y = c( Y1, Y0 )
    Z  = rep( c(1,0), each= length(Y0) )
    B = c( B, B )

    stopifnot( length(Y) == 2*n )


    # Total number of blocks
    K <- length( unique( data_table$B) )

    #Split into big and small blocks
    data.big <- data_table[data_table$n1 > 1 & data_table$n0 > 1, ]
    data.small<-data_table[data_table$n1 == 1 | data_table$n0 == 1, ]
    #Calculate variance for big blocks
    var_big <- sum(data.big$se_ney ^ 2 * (data.big$nk) ^ 2) / n ^ 2
    #I f no small blocks, just use big block variance
    if (nrow(data.small) == 0) {
        hybrid_m_est <- var_big
        hybrid_p_est <- var_big
        plug_in_big_est <- var_big
    } else{
        # Otherwise calculate variance estimates for each method considered
        mod.small <- merge(s.tc.bk, data.small, by = "B")
        var_small <- sum((mod.small$se_ney ^ 2 - mod.small$s.tc.bk / mod.small$nk) * (mod.small$nk) ^ 2) / n ^ 2
        hybrid_m_est <- hybrid_m_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big + var_small
        hybrid_p_est <- hybrid_p_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big + var_small
        plug_in_big_est <- plug_in_big(data.small, data.big) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big
    }
    #Get trt effect estimates and aggregate
    ATE_vec <- data_table$Ybar1 - data_table$Ybar0
    ATE_hat <- sum(ATE_vec * data_table$nk) / n

    #Get linear model estimates (sandwich)
    M0 <- lm(Y ~ Z + as.factor(B))
    ATE_hat_lm <- summary(M0)$coefficients[2, 1]
    var_est_lm <- sandwich::vcovHC(M0, type = "HC1")[2, 2]
    num.c.vec <- data_table[match(B,data_table$B), ]$nk / data_table[match(B, data_table$B), ]$n0
    num.t.vec <- data_table[match(B,data_table$B), ]$nk / data_table[match(B, data_table$B), ]$n1
    weight.vec <- num.c.vec * (1 - Z) + num.t.vec * (Z)
    p.overall <- sum(data_table$n1) / sum(data_table$n1 + data_table$n0)
    weight.vec <- num.c.vec * (1 - Z) * (1 - p.overall) + num.t.vec * (Z) * p.overall

    # Get linear model estimates (Weighted OLS)
    M1 <- lm(Y ~ Z + as.factor(B), weights = weight.vec)
    ATE_hat_weighted_fixed <- summary(M1)$coefficients[2, 1]
    var_est_weighted_fixed <- summary(M1)$coefficients[2, 2] ^ 2

    # Aggregate results into summary table
    methods_list <- c("hybrid_m", "hybrid_p", "plug_in_big", "fixed effects-no int", "weighted regression-fixed effects")
    var_estimates <- c(hybrid_m_est, hybrid_p_est, plug_in_big_est, var_est_lm, var_est_weighted_fixed)
    ATE_estimates <- c(rep(ATE_hat, 3), ATE_hat_lm, ATE_hat_weighted_fixed)
    summary_table <- data.frame(Methods = methods_list, tau = ATE_estimates, SE = sqrt(var_estimates))
    return(summary_table)
}
