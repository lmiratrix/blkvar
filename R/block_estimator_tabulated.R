#' Block variance estimation function.
#'
#' This function takes the block-level summary statistics of a dataset and
#' returns a treatment effect estimate, variance estimate and some summary info
#' about the block structure.
#'
#' @param data.table  Summary statistics of all the blocks in a dataset.  In
#'   particular, this is the output of calc_summary_stats method.
#' @param method The method to be used for variance estimation, defauly
#'   "hybrid_m"
#' @param throw.warnings TRUE means throw warnings if the hybrid estimators are
#'   breaking down due to violation of assumptions.
#' @export
block_estimator_tabulated <- function(data.table, method = c("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small"),
  throw.warnings = TRUE) {

  stopifnot(is.data.frame( data.table))
  stopifnot(all(c("Ybar0","Ybar1","n", "var1","var0", "n1","n0", "se_ney") %in% names(data.table)))

  method <- match.arg(method)
  K <- nrow(data.table)
  data.table$nk <- data.table$n1 + data.table$n0
  n <- sum(data.table$nk)

  #Split into big and small blocks
  data.big <- data.table[data.table$n1 > 1&data.table$n0 > 1, ]
  data.small <- data.table[data.table$n1 == 1 | data.table$n0 == 1, ]
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
    dt <- convert.table.names(data.small)
    var_small <- (calc.RCT.Yes.SE(dt, "superpop.original" , weight = "individual")$SE) ^ 2 * sum(data.small$nk) ^ 2 / n ^ 2
    if (throw.warnings && is.na(var_small)) {
      warning("Need multiple small blocks for rct_yes_small")
    }
  } else if (method == "rct_yes_mod_small") {
    dt <- convert.table.names(data.small)
    var_small <- (calc.RCT.Yes.SE(dt, "superpop", weight = "individual")$SE) ^ 2 * sum(data.small$nk) ^ 2 / n ^ 2
    if (throw.warnings && is.na(var_small)) {
      warning("Need multiple small blocks for rct_yes_small")
    }
  }
  #Get trt effect estimates and aggregate
  tau_vec <- data.table$Ybar1 - data.table$Ybar0
  tau_est <- sum(tau_vec * data.table$nk) / n
  #Get overall variance estimate
  if (method == "rct_yes_all") {
    dt <- convert.table.names(data.table)
    var_est <- (calc.RCT.Yes.SE(dt, "superpop.original" , weight = "individual")$SE) ^ 2
  } else if (method == "rct_yes_mod_all") {
    dt <- convert.table.names(data.table)
    var_est <- (calc.RCT.Yes.SE(dt, "superpop" , weight = "individual")$SE) ^ 2
  } else {
    var_est <- var_small + var_big
  }
  # Table summarizing block sizes
  size_table <- data.table[, c("B", "n1", "n0")]
  #Percent of blocks that are small
  perc_small <- round(sum(data.small$nk) / n * 100, 2)
  return_val<-list(tau_est, var_est, perc_small, size_table)
  names(return_val) <- c("tau_est", "var_est", "percent_small_blocks", "block_sizes")
  return_val$method <- method
  return_val$se_est <- sqrt(var_est)
  class(return_val) <- "var_dat"
  return(return_val)
}