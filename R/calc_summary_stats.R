#' Summarise data by block.
#'
#' Given dataframe or list of variables, return table of stats for each
#' randomization block.
#'
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z vector that indicates if outcome is under treatment or control
#' @param B block ids
#' @param data Dataframe with defined Yobs, Z, and B variables.
#' @param siteID If not null, name of siteID that has randomization blocks
#'   nested inside.
#' @param add.neyman If TRUE, add block-specific SEs using Neyman formula of
#'   $s^2_T/n_t + s^2_C/n_c$ as column in output.
#' @param data Dataframe with defined Yobs, Z, and B variables.
#' @family calc_summary_stats
#'
#' @return dataframe with summary statistics by block
#' @importFrom rlang .data
#' @export

calc_summary_stats <- function(Yobs, Z, B, data = NULL, siteID = NULL, add.neyman = FALSE) {
  if (missing("Z") && is.null(data)) {
    data <- Yobs
  }
  if (is.null(data)) {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
  } else {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    }
    if (!is.null( siteID ) && (length(siteID) == 1)) {
      siteID <- data[[siteID]]
    }
    if (missing( "Z")) {
      stopifnot(all(c("Yobs", "Z", "B") %in% names(data)))
    } else {
      data <- data.frame(Yobs = eval(substitute(Yobs), data),
                         Z = eval(substitute(Z), data),
                         B = eval(substitute(B), data))
    }
  }
  dat <- data
  if (!is.null(siteID)) {
    dat $siteID <- siteID
  }
  if (!is.null(siteID)) {
    sdat <- dat %>% dplyr::group_by(B, Z, siteID)
  } else {
    sdat <- dat %>% dplyr::group_by(B, Z)
  }

  sdat <- sdat %>% dplyr::summarise(n = n(), Ybar = mean(Yobs), var = var(Yobs))
  sdat <- reshape(as.data.frame(sdat), direction = "wide", v.names = c("n", "Ybar", "var"),
    idvar = "B", timevar = "Z", sep= "")
  sdat$n <- sdat$n0 + sdat$n1
  if (add.neyman) {
    sdat <- mutate(sdat, se_ney = sqrt(var1 / n1 + var0 / n0))
  }
  if (any(is.na(sdat$n1) | is.na(sdat$n0))) {
    # sdat$n1[is.na(sdat$n1)] <- 0
    # sdat$n0[is.na(sdat$n0)] <- 0
    sdat <- dplyr::filter(sdat, !is.na(n1), !is.na(n0))
    warning("Some blocks have no treatment or no control units; they are being dropped")
  }
  sdat
}



#' Summarise data by block.
#'
#' Given dataframe or list of variables, return table of stats for each
#' randomization block.
#'
#' @inheritParams calc_summary_stats
#' @param  formula  Formula of form Yobs ~ Z * B (ONLY)
#' @family calc_summary_stats
#' @return dataframe with summary statistics by block
#' @export

calc_summary_stats_formula <- function(formula, data = NULL, siteID = NULL, add.neyman = FALSE) {

    data = make_canonical_data( formula=formula, data=data, siteID=siteID )
    calc_summary_stats(Yobs, Z, B, data = data, siteID = siteID, add.neyman = add.neyman )
}

