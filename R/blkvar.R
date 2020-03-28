#' \code{blkvar} package
#'
#' Blocked Design Variance Estimation
#'
#' See the README on
#' \href{https://github.com/lmiratrix/blkvar#readme}{GitHub}
#'
#' @docType package
#' @name blkvar
#' @importFrom dplyr %>%
#' @importFrom purrr %||%
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines (and .MULTISITE_CANONICAL is a true global variable).
globalVariables(names = c("n", "n1", "var0", "n0", "Yobs", "Z", "B", "mu1", ".weight", "MSE.C", "MSE.T",
  "X1.bar", "X2.bar", "Y0", "Y0.bar", "Y1", "Y1.bar", "Ybar0", "Ybar1", "mu0",
  "nC", "nT", "nj", "p", "prec", "sid", "tau_vec", "var1", "weight", "weight.site", ".MULTISITE_CANONICAL"),
  package = "blkvar", add = FALSE)
