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
  data.table <- calc_summary_stats(Yobs, Z, B, add.neyman = TRUE)
  block_estimator_tabulated(data.table, method = method, throw.warnings = throw.warnings)
}