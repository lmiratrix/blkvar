#' \code{blkvar} package
#'
#' Estimating ATE and Cross-Site Variation for Multisite (block) RCTs
#'
#' This code bundle is a package for estimating impacts and cross-site variation
#' in multi-site or blocked RCTs.
#'
#' It is the compainion package for Pashley & Miratrix (2020) and Miratrix,
#' Weiss & Henderson (2020).
#'
#' It also has a bunch of routines for running multisite simulation scenarios.
#' It is used in a variety of research projects connected with the Miratrix
#' CARES Lab (\url{https://cares.gse.harvard.edu/}) and other research that
#' Miratrix is part of.
#'
#' The github has development versions, with more experimental code. For working
#' (cleaner) versions see CRAN. See
#' \href{https://github.com/lmiratrix/blkvar#readme}{GitHub} for current
#' development version.
#'
#' @references
#' Pashley & Miratrix (2020) Insights on Variance Estimation for
#' Blocked and Matched Pairs Designs on \url{https://arxiv.org/abs/1710.10342}
#'
#' Miratrix, Weiss and Henderson (2020) "An Applied Researcher’s Guide to
#' Estimating Effects from Multisite Individually Randomized Trials: Estimands,
#' Estimators, and Estimates"
#'
#' @docType package
#' @name blkvar
#' @importFrom dplyr %>%
#' @importFrom purrr %||%
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines (and .MULTISITE_CANONICAL is a true global variable).
globalVariables(names = c("n", "n1", "var0", "n0", "Yobs", "Z", "B", "mu1", ".weight", "MSE.C", "MSE.T",
  "X1.bar", "X2.bar", "Y0", "Y0.bar", "Y1", "Y1.bar", "Ybar0", "Ybar1", "mu0",
  "nC", "nT", "nj", "p", "prec", "sid", "ATE_vec", "var1", "weight", "weight.site", ".MULTISITE_CANONICAL"),
  package = "blkvar", add = FALSE)
