#' Compare different estimates of cross site variation.
#'
#' Given a dataframe, use the different methods to pull out point estimates and
#' (if desired) pvalues and return them all.
#'
#' @inheritParams compare_methods
#' @param long.results TRUE means each estimator gets a line in a data.frame.  FALSE gives all as columns in a 1-row dataframe.
#' @param include.testing TOADD
#' @param siteID if blocks B nested in sites, then pass the site indicator.
#'
#' @export
compare_methods_variation <- function(Yobs, Z, B, siteID = NULL, data = NULL, include.testing = TRUE, long.results = FALSE) {
  if (!is.null(data)) {
    if (missing("Yobs")) {
      data <- data.frame(Yobs = data[, 1], Z = data[, 2], B = data[, 3])
      n.tx.lvls <- table(data$Z)
      stopifnot(n.tx.lvls == 2)
      stopifnot(is.numeric(data$Yobs))
    } else {
      if (!is.null(siteID)) {
        siteID <- data[[siteID]]
      }
      data <- data.frame(Yobs = eval( substitute(Yobs), data), Z = eval(substitute(Z), data), B = eval( substitute(B), data))
      data$siteID <- siteID

      # if ( !missing( siteID ) ) {
	  #     data = data.frame( Yobs = eval( substitute( Yobs ), data ),
	  #                        Z = eval( substitute( Z ), data ),
 	  #                        B = eval( substitute( B ), data),
	  #                        siteID = eval( substitute( siteID ), data ) )
	  # } else {
	  #     data = data.frame( Yobs = eval( substitute( Yobs ), data ),
	  #                        Z = eval( substitute( Z ), data ),
	  #                        B = eval( substitute( B ), data) )
	  # }
    }
  } else {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
    if (!is.null(siteID)) {
      data$siteID = siteID
    }
  }

  n <- nrow(data)

  # Quick check that input is correct
  if (is.numeric(data$Z) == FALSE) {
    stop("Treatment indicator should be vector of ones and zeros")
  }
  if ((sum(data$Z == 1) + sum(data$Z == 0)) != n) {
    stop("Treatment indicator should be vector of ones and zeros")
  }
  if (!is.null(siteID)) {
    siteID <- "siteID" #quote( siteID )
  }

  # FIRC model (separate variances)
  FIRC <- estimate_ATE_FIRC(Yobs, Z, B, siteID = siteID, data = data, include.testing = include.testing)

  # FIRC model (with pooled residual variances)
  FIRC.pool <- estimate_ATE_FIRC(Yobs, Z, B, siteID = siteID, data = data, include.testing = include.testing, pool = TRUE)

  # the random-intercept, random-coefficient (RIRC) model
  RIRC <- estimate_ATE_RIRC(Yobs, Z, B, data, include.testing = include.testing)

  # the random-intercept, random-coefficient (RIRC) model
  RIRC.pool <- estimate_ATE_RIRC(Yobs, Z, B, data, include.testing = include.testing, pool = TRUE)
  Qstat <- analysis_Qstatistic(Yobs, Z, B, siteID = siteID, data = data )

  # collect results
  res <- data.frame(tau.hat.FIRC = FIRC$tau.hat, tau.hat.RIRC = RIRC$tau.hat, tau.hat.FIRC.pool = FIRC.pool$tau.hat, tau.hat.RIRC.pool = RIRC.pool$tau.hat)
  res$tau.hat.Q <- Qstat$tau.hat
  
  if (include.testing) {
  	res$pv.FIRC <- FIRC$p.variation
    res$pv.RIRC <- RIRC$p.variation
    res$pv.FIRC.pool <- FIRC.pool$p.variation
    res$pv.RIRC.pool <- RIRC.pool$p.variation
    res$pv.Qstat <- Qstat$p.value
  }

  if (long.results) {
    if (include.testing) {
      res <- data.frame( method = c("FIRC", "RIRC", "FIRC.pool", "RIRC.pool", "Q"), tau.hat = c( as.numeric( res[1:4] ), NA), pv = as.numeric(res[5:9]))
    } else {
      res = data.frame( method = c("FIRC", "RIRC", "FIRC.pool", "RIRC.pool"), tau.hat = as.numeric(c(res[1:4])))
    }
  }
  res
}