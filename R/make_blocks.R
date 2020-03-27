#' Make blocks out of a continuous covariate by cutting it into pieces that are
#' ideally relatively homogenous.
#'
#' Given a dataset, create a bunch of blocks based on passed covariate and
#' return the factor of blocks
#' @param X covariate vector to block on
#' @param method How to block.
#' @param num.blocks If method is small, how many blocks to attempt to make
#' @return Vector with one element per element of `X`
#' @export

make_blocks <- function(X, method = c("small", "pair", "big", "none"), num.blocks) {
  method <- match.arg(method)
  X.orig = X
  X.order <- rank(X, ties.method = "first")
  X <- sort(X)
  if (method == "small") {
    dels <- diff(X)
    dels <- jitter(dels)
    N <- length(X)
    if (!missing(num.blocks)) {
      ct <- sort(dels, decreasing = TRUE)[num.blocks - 1]
    } else {
      ct <- quantile(dels, 2 / 3 + 1 / N)
    }
    cuts <- (dels >= ct) & (c(dels[-1], Inf) < ct) & (c(Inf, dels[ - (N - 1)]) < ct)
    cuts <- 1 + cumsum(cuts)
    B <- paste("B", c(1, cuts), sep = "")
  } else if (method == "pair") {
    B <- paste("B", rep(1:(length(X) / 2), each = 2), sep = "")
  } else if (method == "big") {
    minB <- trunc(length(X) / 4)
    B <- cut(X, quantile(X, seq(0, 1, length.out = minB + 1)), include.lowest = TRUE)
  } else {
    B <- paste("B", rep(1, length(X)), sep = "")
  }
  B[X.order]
}