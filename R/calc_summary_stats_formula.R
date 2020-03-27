#' Summarise data by block.
#'
#' Given dataframe or list of variables, return table of stats for each
#' randomization block.
#'
#' @param data Dataframe with defined Yobs, Z, and B variables.
#' @param siteID If not null, name of siteID that has randomization blocks
#'   nested inside.
#' @param siteID if blocks B nested in sites, then pass the site indicator.
#' @param  formula  Formula of form Yobs ~ Z * B (ONLY)
#' @param  add.neyman TOADD
#'
#' @return dataframe with summary statistics by block
#' @export
calc_summary_stats_formula <- function(formula, data = NULL, siteID = NULL, add.neyman = FALSE) {
  # Figure out the covariates we are using
  if (length(formula.tools::lhs.vars(formula)) != 1 | length(formula.tools::rhs.vars(formula)) != 2) {
    stop("The formula argument must be of the form outcome ~ treatment*block_id.")
  }
  main.vars <- formula.tools::get.vars(formula, data=data)
  if (any(!(main.vars %in% colnames(data)))) {
    browser()
    stop("Some variables in formula are not present in your data.")
  }
  out.name <- formula.tools::lhs.vars(formula)[[1]]
  main.name <- formula.tools::rhs.vars(formula)

  # Copy over the variables to our names
  data$Yobs <- data[[out.name]]
  data$Z <- data[[ main.name[[1]]]]
  data$B <- data[[ main.name[[2]]]]
  if (!is.null(siteID) && (length(siteID) == 1)) {
    data$siteID <- data[[ siteID ]]
  } else {
    data$siteID <- siteID
  }
  calc_summary_stats(Yobs, Z, B, data = data, siteID = siteID, add.neyman = add.neyman )
}
