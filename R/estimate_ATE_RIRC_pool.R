#' Fit Random-intercept, random-slope model
#'
#' This is a wrapper for a simpler lmer call.
#' @inheritParams estimate_ATE_FIRC
#' @importFrom lme4 lmer VarCorr
# #' @export
estimate_ATE_RIRC_pool <- function(Yobs, Z, B, data = NULL, REML = FALSE, include_testing = TRUE, control.formula = NULL) {
  if (!is.null(control.formula)) {
    stopifnot(!is.null( data ))
    stopifnot(!missing("Yobs"))
  }
  if (!is.null(data)) {
    if (missing("Yobs")) {
      data <- data.frame( Yobs = data[[1]], Z = data[[2]], B = data[[3]])
    } else {
      if (!is.null(siteID)) {
        siteID <- data[[siteID]]
        stopifnot(!is.null(siteID))
      }
      d2 <- data
      if (!is.null(siteID)) {
        d2$siteID <- data[[siteID]]
        stopifnot(!is.null(d2$siteID))
      }
      d2$Yobs <- eval(substitute(Yobs), data)
      d2$Z <- eval(substitute(Z), data)
      d2$B <- eval(substitute(B), data)
      data <- d2
      rm(d2)
    }
  } else {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
    if (!is.null(siteID)) {
      data$siteID <- siteID
    }
  }
  stopifnot(length(unique(data$Z)) == 2)
  stopifnot(is.numeric(data$Yobs))
  M0.full <- lme4::lmer( Yobs ~ 1 + Z + (1 + Z|B), data = data, REML = REML)
  if (include_testing) {
    M0.null <- lme4::lmer( Yobs ~ 1 + Z + (1|B), data = data, REML = REML)
    # I _think_ this is what is suggested to handle the boundary by Snijders and Bosker
    td <- deviance(M0.null) - deviance(M0.full)
    pv <- 0.5 * pchisq(td, 2, lower.tail = FALSE) + 0.5 * pchisq(td, 1, lower.tail = FALSE)
  } else {
    pv <- NA
    td <- NA
  }

  # get ATE and SE
  a <- summary(M0.full)
  a <- a$coefficients[2, ]

  # Cross site variation
  tau_hat <- sqrt( VarCorr(M0.full)$B[2, 2])
  res <- list(ATE = a[[1]], SE.ATE = a[[2]], tau_hat = tau_hat, SE_tau = NA)
  if (include_testing) {
    res$p_variation <- pv
    res$deviance <- td
  }
  res
}
