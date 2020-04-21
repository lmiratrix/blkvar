
##
## This file contains the various estimators discussed in Pashley et al
##



# #' Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}




#' Matched pairs type variance function
#'
#' Function that calculates matched pairs type variance.
#' @param data_table data frame containing info for set of blocks of same size
#' @param weighted indicates whether variance should be weighted by number of units
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @noRd

paired_var <- function(data_table, weighted = TRUE) {
  #Vector of trt effect estimates
  ATE_vec <- data_table$Ybar1 - data_table$Ybar0
  var_est <- sum((ATE_vec - mean(ATE_vec)) ^ 2) / ((length(ATE_vec) - 1) * length(ATE_vec))
  if (weighted){
    var_est <- var_est * sum(data_table$nk) ^ 2
  }
  return(var_est)
}

# #' Aggregated matched pairs type variance estimator
# #'
# #' Function to calculate hybrid_m variance estimator.
# #' @importFrom stats aggregate lm quantile rnorm sd var

hybrid_m_small <- function(data.small) {
  #First split into groups of blocks of same size
  by_size <- split(data.small, data.small$nk, drop = FALSE)
  #Estimate variance in each of these groups
  var_est <- sum(sapply(by_size, paired_var)) / sum(data.small$nk) ^ 2
  if (is.nan(var_est)) {
    var_est <- NA
  }
  return(var_est)
}

# #' Pooled matched pairs type variance estimator
# #'
# #' Function to calculate hybrid_p variance estimator.
# #' @param data.small data frame containg info for small blocks
# #' @importFrom stats aggregate lm quantile rnorm sd var

hybrid_p_small <- function(data.small) {
  #First check we can use this estimator
  n <- sum(data.small$nk)
  if (max(data.small$nk) >= n/ 2 ){
    return(NA)
  }
  #Vector of treatment effect estimates
  ATE_vec <- data.small$Ybar1 - data.small$Ybar0
  #Constant that adjusts for varying block sizes
  constant <- data.small$nk ^ 2 / ((n - 2 * data.small$nk) * (n + sum(data.small$nk ^ 2 / (n - 2 * data.small$nk))))
  var_est <- sum(constant * (ATE_vec - sum(ATE_vec * data.small$nk) / n) ^2 )
  return(var_est)
}

# #' Plug in type variance estimator
# #'
# #'Function to calculate variance estimate using average big block variances as plug-ins.
# #'
# #' @param data.small data frame containg info for small blocks
# #' @param data.big data frame containg info for small blocks
# #' @importFrom stats aggregate lm quantile rnorm sd var

plug_in_big <- function(data.small, data.big) {
  n <- sum(data.small$nk)
  #Get average variances in big blocks
  var0_avg <- mean(data.big$var0)
  var1_avg <- mean(data.big$var1)
  #Put in average variances in missing pieces of small blocks
  var.small.c <- as.numeric(data.small$n0 == 1) * var0_avg
  var.small.t <- as.numeric(data.small$n1 == 1) * var1_avg
  #Not all small blocks are missing things
  var.big.c <- as.numeric(data.small$n0 > 1) * data.small$var0 / data.small$n0
  var.big.t <- as.numeric(data.small$n1 > 1) * data.small$var1 / data.small$n1
  var.big.c[is.na(var.big.c)] <- 0
  var.big.t[is.na(var.big.t)] <- 0
  #Combine all variance terms
  var.c <- (var.small.c + var.big.c) * data.small$nk ^ 2
  var.t <- (var.small.t + var.big.t) * data.small$nk ^ 2
  var_est <- sum(var.c + var.t) / n ^ 2
  return(var_est)
}


# #' Print method for block_estimator output
# #'
# #' Function that compares difference variance estimators for blocked designs.
# #' @param x output from block_estimator (class=var_data)
# #' @param ... further arguments to match print
# #' @importFrom stats aggregate lm quantile rnorm sd var
# #' @export
print.var_dat <- function(x, ... ) {
  scat("Randomization Inference Treatment Estimate (method = %s)\n", x$method)
  SE <- sqrt(x$var_est)
  scat("Estimate of ATE: %.2f \t(SE: %.2f)\n", x$ATE_hat, SE)
  scat("Nominal 95 Confidence Interval: %.2f - %.2f\n",  x$ATE_hat - 2 * SE, x$ATE_hat + 2 * SE)
  cat(paste(x$percent_small_blocks, "% of units are in small blocks", sep=""), "\n")
  bsize <- x$block_sizes
  scat("\nBlock Sizes (%d blocks):", nrow(bsize))
  bsize$name <- paste0(bsize$n1, "-", bsize$n0)
  tb <- table(bsize$name)
  tb <- sort(tb, decreasing = TRUE)
  print(tb)
}




#' Block variance estimation function.
#'
#' This function takes the block-level summary statistics of a dataset and
#' returns a treatment effect estimate, variance estimate and some summary info
#' about the block structure.
#'
#' @param summary_stats  Summary statistics of all the blocks in a dataset.  In
#'   particular, this is the output of calc_summary_stats method.
#' @param method The method to be used for variance estimation, defauly
#'   "hybrid_m"
#' @param throw.warnings TRUE means throw warnings if the hybrid estimators are
#'   breaking down due to violation of assumptions.
#' @export
block_estimator_tabulated <- function(summary_stats, method = c("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small"),
                                      throw.warnings = TRUE) {

    stopifnot(is.data.frame( summary_stats))
    stopifnot(all(c("Ybar0","Ybar1","n", "var1","var0", "n1","n0", "se_ney") %in% names(summary_stats)))

    method <- match.arg(method)
    K <- nrow(summary_stats)
    summary_stats$nk <- summary_stats$n1 + summary_stats$n0
    n <- sum(summary_stats$nk)

    #Split into big and small blocks
    data.big <- summary_stats[summary_stats$n1 > 1&summary_stats$n0 > 1, ]
    data.small <- summary_stats[summary_stats$n1 == 1 | summary_stats$n0 == 1, ]
    #Get variance for big blocks
    var_big <- sum(data.big$se_ney ^ 2 * (data.big$nk) ^ 2) / n ^ 2
    #If no small blocks, set small block variance to 0
    if (nrow(data.small) == 0) {
        var_small <- 0
    } else if (method == "hybrid_m" ) {
        #Estimate variance according to methof given
        var_small <- hybrid_m_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2
        if (throw.warnings && is.na(var_small)) {
            warning("Need multiple small blocks of the same size for hybrid_m")
        }
    } else if (method == "hybrid_p") {
        var_small <- hybrid_p_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2
        if (throw.warnings && is.na(var_small)) {
            warning("Largest small block is more than half the small block units in hybrid_p so no variance calc possible")
        }
    } else if (method == "plug_in_big") {
        var_small <- plug_in_big(data.small, data.big) * sum(data.small$nk) ^ 2 / n ^ 2
        if (throw.warnings && is.na(var_small)) {
            warning("Need some big blocks for plug_in_big")
        }
    } else if (method == "rct_yes_small") {
        var_small <- (estimate_ATE_design_based_from_stats(summary_stats, method="superpop.original" , weight = "individual")$SE) ^ 2 * sum(data.small$nk) ^ 2 / n ^ 2
        if (throw.warnings && is.na(var_small)) {
            warning("Need multiple small blocks for rct_yes_small")
        }
    } else if (method == "rct_yes_mod_small") {
        var_small <- (estimate_ATE_design_based_from_stats(summary_stats, method="superpop", weight = "individual")$SE) ^ 2 * sum(data.small$nk) ^ 2 / n ^ 2
        if (throw.warnings && is.na(var_small)) {
            warning("Need multiple small blocks for rct_yes_small")
        }
    }
    #Get trt effect estimates and aggregate
    ATE_vec <- summary_stats$Ybar1 - summary_stats$Ybar0
    ATE_hat <- sum(ATE_vec * summary_stats$nk) / n
    #Get overall variance estimate
    if (method == "rct_yes_all") {
        var_est <- (estimate_ATE_design_based_from_stats(summary_stats,
                                                         method="superpop.original" , weight = "individual")$SE) ^ 2
    } else if (method == "rct_yes_mod_all") {
        var_est <- (estimate_ATE_design_based_from_stats(summary_stats, method="superpop" , weight = "individual")$SE) ^ 2
    } else {
        var_est <- var_small + var_big
    }
    # Table summarizing block sizes
    size_table <- summary_stats[, c("B", "n1", "n0")]
    #Percent of blocks that are small
    perc_small <- round(sum(data.small$nk) / n * 100, 2)
    return_val<-list(ATE_hat, var_est, perc_small, size_table)
    names(return_val) <- c("ATE_hat", "var_est", "percent_small_blocks", "block_sizes")
    return_val$method <- method
    return_val$se_est <- sqrt(var_est)
    class(return_val) <- "var_dat"
    return(return_val)
}





#' Block variance estimation function.
#'
#' This function takes observed data and returns a treatment effect estimate,
#' variance estimate and some summary info about the block structure.
#'
#' @param Yobs vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data matrix of Y, Z, B, as alternative to using vectors
#' @param method is method of variance estimation, defauly "hybrid_m"
#' @param throw.warnings TRUE means throw warnings if the hybrid estimators are
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
block_estimator <- function(Yobs, Z, B, data = NULL, method = c("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all",
                                                                "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small"), throw.warnings = TRUE) {
    if (!is.null(data)) {
        if (missing( "Yobs")) {
            Yobs <- data[, 1]
            Z <- data[, 2]
            B <- data[, 3]
        } else {
            Yobs <- eval(substitute(Yobs), data)
            Z <- eval(substitute(Z), data)
            B <- eval(substitute(B), data)
        }
    } else {
        if (is.data.frame(Yobs)) {
            stopifnot(is.null(data))
            B <- Yobs$B
            Z <- Yobs$Z
            Yobs <- Yobs$Yobs
        }
    }

    n <- length(Yobs)
    # Quick test that input is correct
    if (is.numeric(Z) == FALSE) {
        return("Treatment indicator should be vector of ones and zeros")
    }
    if ((sum(Z ==1 ) + sum(Z == 0)) != n) {
        return("Treatment indicator should be vector of ones and zeros")
    }

    # Get summary info
    summary_table <- calc_summary_stats(Yobs, Z, B, add.neyman = TRUE)
    block_estimator_tabulated(summary_table, method = method, throw.warnings = throw.warnings)
}
