


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



#' Table of data for simulations
#'
#' Function that returns a summary of block level true values for sims.
#'
#' @param Y vector of all outcomes
#' @param B block ids
#' @param data alternatively is matrix of Y,B
#' @param p.mat  matrix with first column B,second column prop treated in that block, p
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
block_data_sim <- function(Y, B, p.mat, data = NULL) {
    if (!is.null(data)) {
        Y <- data[, 1]
        B <- data[, 2]
    }
    n <- length(Y)

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
    treated_mat <- cbind(Y1, B)
    control_mat <- cbind(Y0, B)
    Y1_matrix <- aggregate(list(Ybar1 = treated_mat[, 1]),
                           list(B = treated_mat[, 2]), FUN = mean)
    Y0_matrix <- aggregate(list(Ybar0 = control_mat[, 1]),
                           list(B = control_mat[, 2]), FUN = mean)
    Ybar_matrix <- merge(Y1_matrix, Y0_matrix, by = "B")
    var1_matrix <- aggregate(list(var1 = treated_mat[, 1]),
                             list(B = treated_mat[, 2]), FUN = var)
    var0_matrix <- aggregate(list(var0 = control_mat[, 1]),
                             list(B = control_mat[, 2]), FUN = var)
    var_matrix <- merge(var1_matrix, var0_matrix, by = "B")
    overall_mat <- merge(n_ctk_matrix, Ybar_matrix, by = "B")
    overall_mat <- merge(overall_mat, var_matrix, by = "B")
    overall_mat$se_ney <- sqrt(overall_mat$var1 / overall_mat$n1 +
                                   overall_mat$var0 / overall_mat$n0)
    drops <- c("n_k", "p")
    overall_mat <- overall_mat[ , !(names(overall_mat) %in% drops)]
    return(overall_mat)
}





#' Function to compare estimators using the true values (i.e., whole science
#' table).  This is for simulaton studies where we know all the potential
#' outcomes.
#'
#' Function that returns some variance function estimates based on all potential
#' outcomes.
#'
#' @param Y0 Vector of control potential outcomes.
#' @param Y1 Vector of treatment potential outcomes
#' @param B block ids
#' @param data Dataframe holding Y0,Y1,B
#' @param p.mat  matrix with two columns, one row per block.  First column is B,
#'   second column prop treated in that block, ("p")
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @importFrom msm deltamethod
#' @export
compare_methods_oracle <- function(Y0, Y1, B, data = NULL, p.mat = NULL ) {

    if (!is.null(data)) {
            d2 <- data
            d2$Y0 <- eval(substitute(Y0), data)
            d2$Y1 <- eval(substitute(Y1), data)
            d2$B <- eval(substitute(B), data)
            data <- d2
            rm(d2)
    } else {
        data <- data.frame(Y0 = Y0, Y1=Y1, B = B)
    }
    Y0 = data$Y0
    Y1 = data$Y1
    B = data$ B

    n <- nrow( data )


    #Get data into table
    s.tc.bk <- aggregate(list(s.tc.bk = Y1 - Y0), list(B = B[1:(n / 2)]), FUN = s_tc_func)
    data_table <- block_data_sim(Y, B, p.mat)
    K <- max(data_table$B)
    n <- sum(data_table$n1) + sum(data_table$n0)
    data_table$nk <- data_table$n1 + data_table$n0
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
        hybrid_m_est <- hybrid_m_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 +
            var_big + var_small
        hybrid_p_est < -hybrid_p_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 +
            var_big + var_small
        plug_in_big_est <- plug_in_big(data.small, data.big) * sum(data.small$nk) ^ 2 / n ^ 2 +
            var_big
    }
    #Get trt effect estimates and aggregate
    tau_vec < -data_table$Ybar1 - data_table$Ybar0
    tau_est <- sum(tau_vec * data_table$nk) / n

    #Get linear model estimates (sandwich)
    d2 = bind_rows( data, data )
    d2$Z = rep( c(0,1), each=nrow(data) )
    M0 <- lm(Y ~ Z + as.factor(B), data=d2)
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
