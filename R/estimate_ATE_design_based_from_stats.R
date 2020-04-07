#' Implementation of the formula found in RCT-YES documentation.
#'
#' This can handle a superpopulation model, non-clustered but blocked.
#'
#' Taken from page 83 of Shochet RCT-YES paper (eq 6.25).
#' The `superpop` variant is a modification of the original 'superpop.original',
#' pulling the weights from inside the squared term to outside. This method was
#' suggested in personal correspondance with Schochet.  If the weights are not
#' all 1, this can make a difference.
#'
#' @param sum_tab Table of summary statistics by block, from, e.g.,
#'   `block.data()`
#' @param siteID Vector of site IDs if there are randomization blocks nested in
#'    site that should be aggregated (will change results for site weighting only).
#' @param weight Individual weight (i.e., number of individuals in each block)
#'   or site weight (average site estimates (which will be considered block
#'   estimates if siteID is null)).
#' @param method finite, superpop, or superpop2 to give SEs that either capture
#'   uncertainty due to targeting a superpopulation quantity or not.
#' @return dataframe with calculated impacts and standard errors.
#' @importFrom magrittr "%>%"
#' @export

estimate_ATE_design_based_from_stats <- function(sum_tab,
                                                 siteID = NULL,
                                                 method = c( "finite", "superpop", "superpop.original"),
                                                 weight = c("individual", "site")) {
  tau_hat <- NA
  stopifnot(is.data.frame(sum_tab))
  stopifnot(all(c("Ybar0","Ybar1","n", "var1","var0", "n1","n0" ) %in% names( sum_tab)))
  method <- match.arg(method)
  weight <- match.arg(weight)
  h <- nrow(sum_tab)
  # Get our block weights depending on target estimand
  if (weight == "individual") {
    sum_tab$.weight = sum_tab$n
  } else {
    if (!is.null(siteID)) {
      sum_tab <- sum_tab %>% dplyr::group_by(!!as.name(siteID)) %>% dplyr::mutate(.weight = n / sum(n)) %>% ungroup()
    } else {
      sum_tab$.weight <- rep(1, h)
    }
  }

  # calculate individual block treatment impact estimates
  sum_tab <- mutate( sum_tab, tau_hat_b = Ybar1 - Ybar0)

  # calculate overall ATE estimate by taking a weighted average of the
  # individual
  tau_hat <- with(sum_tab, sum(tau_hat_b * .weight) / sum(.weight))
  # Now do the SEs.
  if (method == "finite") {
    # finite pop (Neyman)
    w.tot <- sum(sum_tab$.weight)
    # calculate SEs for each block by itself
    sum_tab <- mutate(sum_tab, block.vars = (var1 / n1 + var0 / n0 ))
    # and then take a weighted sum of these
    var <- sum (sum_tab$.weight ^ 2 * sum_tab$block.vars) / w.tot ^ 2
    SE <- sqrt(var)
  } else {  # superpopulation!
    # First aggregate to get sites, if needed
    if (!is.null(siteID)) {
      sum_tab <- sum_tab %>% group_by(!!as.name(siteID)) %>%
          dplyr::summarise(tau_hat_b = sum(.data$tau_hat_b * .data$.weight) / sum(.data$.weight),
                           .weight = sum(.data$.weight))
      h <- nrow( sum_tab )
    }
    # Calculate average weight across sites
    wbar <- mean(sum_tab$.weight)
    if (method == "superpop.original") {
      # This is the formula 6.25
      asyVar <- sum((sum_tab$.weight * sum_tab$tau_hat_b - wbar * tau_hat) ^ 2) / ((h - 1) * h * wbar ^ 2)
    } else if (method == "superpop") {
      # This is based on the email chain with Weiss, Pashley, etc.
      asyVar <- sum(sum_tab$.weight ^ 2 * (sum_tab$tau_hat_b - tau_hat) ^ 2) / ((h -1 ) * h * wbar ^ 2)
    }
    SE <- sqrt(asyVar)
  }
  data.frame(tau_hat = tau_hat, SE = SE, weight = weight, method = method, stringsAsFactors = FALSE)
}
