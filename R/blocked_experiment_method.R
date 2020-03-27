
##
## This file contains the various estimators discussed in Pashley et al
##

utils::globalVariables(c("."))


# #' Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}

convert.table.names = function( data.table ) {
  names(data.table) <- c( "B", "n1", "n0", "Ybar1", "Ybar0", "var.1", "var.0", "se_ney")
  data.table$n <- with(data.table, n0 + n1)
  data.table
}



#' Matched pairs type variance function
#'
#' Function that calculates matched pairs type variance.
#' @param data.table data frame containing info for set of blocks of same size
#' @param weighted indicates whether variance should be weighted by number of units
#' @importFrom stats aggregate lm quantile rnorm sd var

paired_var <- function(data.table, weighted = TRUE) {
  #Vector of trt effect estimates
  tau_vec <- data.table$Ybar1 - data.table$Ybar0
  var_est <- sum((tau_vec - mean(tau_vec)) ^ 2) / ((length(tau_vec) - 1) * length(tau_vec))
  if (weighted){
    var_est <- var_est * sum(data.table$nk) ^ 2
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
  tau_vec <- data.small$Ybar1 - data.small$Ybar0
  #Constant that adjusts for varying block sizes
  constant <- data.small$nk ^ 2 / ((n - 2 * data.small$nk) * (n + sum(data.small$nk ^ 2 / (n - 2 * data.small$nk))))
  var_est <- sum(constant * (tau_vec - sum(tau_vec * data.small$nk) / n) ^2 )
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
  scat("Estimate of tau: %.2f \t(SE: %.2f)\n", x$tau_est, SE)
  scat("Nominal 95 Confidence Interval: %.2f - %.2f\n",  x$tau_est - 2 * SE, x$tau_est + 2 * SE)
  cat(paste(x$percent_small_blocks, "% of units are in small blocks", sep=""), "\n")
  bsize <- x$block_sizes
  scat("\nBlock Sizes (%d blocks):", nrow(bsize))
  bsize$name <- paste0(bsize$n1, "-", bsize$n0)
  tb <- table(bsize$name)
  tb <- sort(tb, decreasing = TRUE)
  print(tb)
}