#' \code{blkvar} package
#'
#' Estimating ATE and Cross-Site Variation for Multisite (block) RCTs
#'
#' This code bundle is a package for estimating impacts and cross-site variation
#' in multi-site or blocked RCTs.
#'
#' It is the companion package for three papers on blocking and multisite experiments (see references below).
#'
#'
#' It also has a bunch of routines for running multisite simulation scenarios.
#' It is used in a variety of research projects connected with the Miratrix
#' CARES Lab (\url{https://cares.gse.harvard.edu/}) and other research that
#' Miratrix is part of.
#'
#' While planned for posting on CRAN, the current version is currently hosted on GitHub. See
#' \href{https://github.com/lmiratrix/blkvar#readme}{GitHub} for current
#' development version.
#'
#' @references
#'
#' Miratrix, L. W., Weiss, M. J., & Henderson, B. (2021). An applied researcher’s guide to estimating effects from multisite individually randomized trials: Estimands, estimators, and estimates. Journal of Research on Educational Effectiveness, 14(1), 270-308.
#'
#' Pashley, N. E., & Miratrix, L. W. (2021). Insights on variance estimation for blocked and matched pairs designs. Journal of Educational and Behavioral Statistics, 46(3), 271-296.
#'
#' Pashley, N. E., & Miratrix, L. W. (In press). Block what you can, except when you shouldn't. Journal of Educational and Behavioral Statistics
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
