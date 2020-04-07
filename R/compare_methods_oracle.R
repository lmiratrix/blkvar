

#'
#' #' Summarise simulation data by block (where we know both Y0 and Y1)
#' #' @param data Dataframe with defined Y0, Y1, and B variables.
#' #' @param Y0 name of Y0 column
#' #' @param Y1 name of Y1 column
#' #' @param Z vector that indicates if outcome is under treatment or control
#' #' @param B block ids
#' #' @return dataframe with summary statistics by block
#' #'
#' #' @export
# calc_summary_stats_oracle_retired <- function (data, Y0 = "Y0", Y1 = "Y1", Z = "Z", B = "B") {
#     data <- dplyr::rename(data,
#                           Y0 = !!rlang::sym(Y0),
#                           Y1 = !!rlang::sym(Y1),
#                           Z = !!rlang::sym(Z),
#                           B = !!rlang::sym(B))
#     sdat <- data %>%
#         dplyr::group_by( B ) %>%
#         dplyr::summarise(n = n(),
#                          mu0 = mean(Y0),
#                          mu1 = mean(Y1),
#                          tau = mu1 - mu0,
#                          sd0 = sd(Y0),
#                          sd1 = sd(Y1),
#                          corr = cor(Y0, Y1))
#     as.data.frame(sdat)
# }



#' Make table of data for simulations
#'
#' Function that returns a summary of block level true values for simulations
#' (so works with full schedule of potential outcomes).
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

    stopifnot( is.null( p_mat ) || is.null(Z) )

    if ( !is.null( Z ) ) {
        stopifnot( sum( Z == 0 ) + sum( Z == 1 ) == length( Y0 ) )
        p_mat = aggregate( list( p=Z ), list( B=B ), mean )
    }

    stopifnot( length(Y0) == length(Y1) )
    stopifnot( length(Y0) == length(B) )

    Y = c( Y1, Y0 )
    Z  = rep( c(1,0), each= length(Y0) )
    B = c( B, B )

    n <- length(Y)

    #First convert block ids into numbers
    B <- as.numeric( factor(B) )
    p_mat$B <- as.numeric( factor(p_mat[,1]) )
    p_mat$p = p_mat[,2]

    #Get number of units assigned to each treatment
    #In each block
    n_matrix <- aggregate(list(n_k = B), list(B = B), FUN = length)
    n_matrix$n_k <- n_matrix$n_k / 2
    n_ctk_matrix <- merge(n_matrix, p_mat, by = "B")
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
#' @noRd
s_tc_func <- function(tau_vec) {
    s.tc<-var(tau_vec)
    return(s.tc)
}



#' Compare estimators using the true values (i.e., full schedule of potential outcomes)
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
#' @param p_mat  matrix with two columns, one row per block.  First column is
#'   the categorical covariate denoting blocks, with same levels as in B, above.
#'   Second column is proportion treated in that block.
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

    s.tc.bk <- aggregate(list(s.tc.bk = Y1 - Y0), list(B = B), FUN = s_tc_func)

    n <- sum(data_table$n1) + sum(data_table$n0)
    data_table$nk <- data_table$n1 + data_table$n0

    Y = c( Y1, Y0 )
    Z  = rep( c(1,0), each= length(Y0) )
    B = c( B, B )

    stopifnot( length(Y) == 2*n )


    # Total number of blocks
    K <- max(data_table$B)

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
        combine_table <- merge(s.tc.bk, data_table, by = "B")
        mod.small <- combine_table[data_table$n1 == 1 | data_table$n0 == 1, ]
        var_small <- sum((mod.small$se_ney ^ 2 - mod.small$s.tc.bk / mod.small$nk) * (mod.small$nk) ^ 2) / n ^ 2
        hybrid_m_est <- hybrid_m_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big + var_small
        hybrid_p_est <- hybrid_p_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big + var_small
        plug_in_big_est <- plug_in_big(data.small, data.big) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big
    }
    #Get trt effect estimates and aggregate
    tau_vec <- data_table$Ybar1 - data_table$Ybar0
    tau_est <- sum(tau_vec * data_table$nk) / n

    #Get linear model estimates (sandwich)
    M0 <- lm(Y ~ Z + as.factor(B))
    tau_est_lm <- summary(M0)$coefficients[2, 1]
    var_est_lm <- sandwich::vcovHC(M0, type = "HC1")[2, 2]
    blk.id <- as.numeric(as.factor(B))
    num.c.vec <- data_table[match(blk.id,data_table$B), ]$nk / data_table[match(blk.id, data_table$B), ]$n0
    num.t.vec <- data_table[match(blk.id,data_table$B), ]$nk / data_table[match(blk.id, data_table$B), ]$n1
    weight.vec <- num.c.vec * (1 - Z) + num.t.vec * (Z)
    p.overall <- sum(data_table$n1) / sum(data_table$n1 + data_table$n0)
    weight.vec <- num.c.vec * (1 - Z) * (1 - p.overall) + num.t.vec * (Z) * p.overall

    # Get linear model estimates (Weighted OLS)
    M1 <- lm(Y ~ Z + as.factor(B), weights = weight.vec)
    tau_est_weighted_fixed <- summary(M1)$coefficients[2, 1]
    var_est_weighted_fixed <- summary(M1)$coefficients[2, 2] ^ 2

    # Aggregate into summary table
    methods_list <- c("hybrid_m", "hybrid_p", "plug_in_big", "fixed effects-no int", "weighted regression-fixed effects")
    var_estimates <- c(hybrid_m_est, hybrid_p_est, plug_in_big_est, var_est_lm, var_est_weighted_fixed)
    tau_estimates <- c(rep(tau_est, 3), tau_est_lm, tau_est_weighted_fixed)
    summary_table <- data.frame(Methods = methods_list, tau = tau_estimates, SE = sqrt(var_estimates))
    return(summary_table)
}
