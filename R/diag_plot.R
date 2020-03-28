#' Make a diagnostic plot
#'
#' Plot block-specific variance estimates versus size of treatment group.
#'
#' @param Y vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data matrix of Y, Z, B, as alternative to using vectors
#' @importFrom graphics axis par plot
make_diagnostic_plot <- function(Y, Z, B, data = NULL) {
  if (!is.null(data)) {
    Y < -data[, 1]
    Z <- data[, 2]
    B <- data[, 3]
  }
  n <- length(Y)
  #Quick test that input is correct
  if (is.numeric(Z) == FALSE) {
    return("Treatment indicator should be vector of ones and zeros")
  }
  if ((sum(Z ==1 ) + sum(Z == 0)) != n) {
    return("Treatment indicator should be vector of ones and zeros")
  }
  data_table <- calc_summary_stats(Y, Z, B)
  #data_table$nk<-data_table$n1+data_table$n0
  var1_plot <- plot(data_table$n1[!is.na(data_table$var1)], data_table$var1[!is.na(data_table$var1)], xlab = "Number treated",
    ylab = "Estimated variance for treated units", xaxt = "n")
  axis(1, at = 1:max(data_table$n1[!is.na(data_table$var1)]))
  var0_plot <- plot(data_table$n0[!is.na(data_table$var0)], data_table$var0[!is.na(data_table$var0)], xlab = "Number in control",
    ylab = "Estimated variance for control units", xaxt = "n")
  axis(1, at = 1:max(data_table$n0[!is.na(data_table$var0)]))
  par(mfrow = c(1, 2))
  var1_plot
  var0_plot
}
