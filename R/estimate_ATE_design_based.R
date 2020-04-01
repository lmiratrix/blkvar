#' Implementation of the formula found in RCT-YES documentation.
#'
#' This can handle a superpopulation model, non-clustered but blocked.
#'
#' Taken from page 83 of Shochet RCT-YES paper (eq 6.25).
#'
#' The `superpop` variant is a modification of the original 'superpop.original',
#' pulling the weights from inside the squared term to outside. This method was
#' suggested in personal correspondance with Schochet.  If the weights are not
#' all 1, this can make a difference.
#' @param formula Input formula for analysis
#' @param control.formula What variables to control for, in the form of "~ X1 + X2".
#' @param data Dataframe with defined Yobs, Z, and B variables.
#' @param siteID Vector of site IDs if there are randomization blocks nested in
#'    site that should be aggregated (will change results for site weighting only).
#' @param weight Individual weight (i.e., number of individuals in each block)
#'   or site weight (average site estimates (which will be considered block
#'   estimates if siteID is null)).
#' @param method finite, superpop, or superpop2 to give SEs that either capture
#'   uncertainty due to targeting a superpopulation quantity or not.
#' @return dataframe with calculated impacts and standard errors.
#' @export

estimate_ATE_design_based <- function(formula,
                                      control.formula = NULL,
                                      siteID = NULL,
                                      data,
                                      method = c("finite", "superpop", "superpop.original"),
                                      weight = c("individual", "site")) {
 if (is.null(control.formula)) {
   data_table <- calc_summary_stats_formula(formula, data=data, siteID=siteID, add.neyman = TRUE)
   estimate_ATE_design_based_from_stats( data_table, siteID = siteID, method = method, weight = weight)
  } else {
    estimate_ATE_design_based_adjusted(formula = formula,
                                       control.formula = control.formula,
                                       siteID = siteID, data = data,
                                       method = method, weight = weight)
  }
}
