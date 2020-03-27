#' Make diagnostic plot
#'
#' Function that plots variance estimates versus size of treatment group.
#'
#' @param Y vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data matrix of Y, Z, B, as alternative to using vectors
#' @importFrom graphics axis par plot
diag_plot <- function(Y, Z, B, data = NULL) {
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
  data.table <- calc_summary_stats(Y, Z, B)
  #data.table$nk<-data.table$n1+data.table$n0
  var1_plot <- plot(data.table$n1[!is.na(data.table$var1)], data.table$var1[!is.na(data.table$var1)], xlab = "Number treated",
    ylab = "Estimated variance for treated units", xaxt = "n")
  axis(1, at = 1:max(data.table$n1[!is.na(data.table$var1)]))
  var0_plot <- plot(data.table$n0[!is.na(data.table$var0)], data.table$var0[!is.na(data.table$var0)], xlab = "Number in control",
    ylab = "Estimated variance for control units", xaxt = "n")
  axis(1, at = 1:max(data.table$n0[!is.na(data.table$var0)]))
  par(mfrow = c(1, 2))
  var1_plot
  var0_plot
}