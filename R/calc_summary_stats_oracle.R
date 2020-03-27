#' Summarise simulation data by block (where we know both Y0 and Y1)
#' @param data Dataframe with defined Y0, Y1, and B variables.
#' @param Y0 name of Y0 column
#' @param Y1 name of Y1 column
#' @param Z vector that indicates if outcome is under treatment or control
#' @param B block ids
#' @return dataframe with summary statistics by block
#' @export
calc_summary_stats_oracle <- function (data, Y0 = "Y0", Y1 = "Y1", Z = "Z", B = "B") {
  data <- rename(data, Y0 = !!rlang::sym(Y0), Y1 = !!rlang::sym(Y1), Z = !!rlang::sym(Z), B = !!rlang::sym(B))
  sdat <- data %>% dplyr::group_by( B ) %>%
  dplyr::summarise(n = n(), mu0 = mean(Y0), mu1 = mean(Y1), tau = mu1 - mu0, sd0 = sd(Y0), sd1 = sd(Y1), corr = cor(Y0, Y1))
  as.data.frame(sdat)
}