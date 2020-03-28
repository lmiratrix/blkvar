#' \code{blkvar} package
#'
#' Estimating ATE and Cross-Site Variation for Multisite (block) RCTs
#'
#' This code bundle is a package for estimating impacts and cross-site variation
#' in multi-site or blocked RCTs.
#'
#' It is the compainion package for Pashley & Miratrix (2020) and Miratrix &
#' Weiss (2020).
#'
#' It also has a bunch of routines for running multisite simulation scenarios.
#' It is used in a variety of research projects connected with the Miratrix
#' CARES Lab (https://cares.gse.harvard.edu/) and other research that Miratrix
#' is part of.
#'
#' For working (cleaner) versions see CRAN.
#'
#' See \href{https://github.com/lmiratrix/blkvar#readme}{GitHub} for current
#' development version.
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
