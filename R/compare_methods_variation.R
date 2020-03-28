#' Compare different estimates of cross site variation.
#'
#' Given a dataframe, use the different methods to pull out point estimates for cross site variation and
#' (if desired) pvalues and return them all.
#'
#' @inheritParams compare_methods
#' @param long_results TRUE means each estimator gets a line in a data.frame.
#'                FALSE gives all as columns in a 1-row dataframe.
#' @param include_testing Include tests for null of no cross site variation.
#' @param include_Q_estimate Include a estimate for tau_hat from Q method (involves calculating the full confidence interval so potentially time intensive).
#' @param siteID if blocks B nested in sites, then pass the site indicator.
#'
#' @export
compare_methods_variation <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                                      include_testing = TRUE,
                                      include_Q_estimate = TRUE,
                                      long_results = FALSE ) {
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
      data <- data.frame(Yobs = eval( substitute(Yobs), data),
                         Z = eval(substitute(Z), data),
                         B = eval( substitute(B), data))
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
  FIRC <- estimate_ATE_FIRC(Yobs, Z, B, siteID = siteID, data = data, include_testing = include_testing)

  # FIRC model (with pooled residual variances)
  FIRC_pool <- estimate_ATE_FIRC(Yobs, Z, B, siteID = siteID, data = data,
                                 include_testing = include_testing, pool = TRUE)

  # the random-intercept, random-coefficient (RIRC) model
  RIRC <- estimate_ATE_RIRC(Yobs, Z, B, data,
                            include_testing = include_testing)

  # the random-intercept, random-coefficient (RIRC) model
  RIRC_pool <- estimate_ATE_RIRC(Yobs, Z, B, data,
                                 include_testing = include_testing,
                                 pool = TRUE)

  Qstat <- analysis_Qstatistic(Yobs, Z, B, siteID = siteID,
                               data = data, calc.CI = include_Q_estimate )

  # collect results
  res <- data.frame(tau_hat_FIRC = FIRC$tau_hat,
                    tau_hat_RIRC = RIRC$tau_hat,
                    tau_hat_FIRC_pool = FIRC_pool$tau_hat,
                    tau_hat_RIRC_pool = RIRC_pool$tau_hat,
                    tau_hat_Q = Qstat$tau_hat )

  if (include_testing) {
  	res$pv_FIRC <- FIRC$p_variation
    res$pv_RIRC <- RIRC$p_variation
    res$pv_FIRC_pool <- FIRC_pool$p_variation
    res$pv_RIRC_pool <- RIRC_pool$p_variation
    res$pv_Qstat <- Qstat$p.value
  }

  if (long_results) {
    if (include_testing) {
      res <- data.frame( method = c("FIRC", "RIRC", "FIRC_pool", "RIRC_pool", "Q"),
                         tau_hat = c( as.numeric( res[1:5] ) ),
                         pv = as.numeric(res[6:10]))
    } else {
      res = data.frame( method = c("FIRC", "RIRC", "FIRC_pool", "RIRC_pool", "Q"),
                        tau_hat = as.numeric(c(res[1:5])))
    }
  }
  res
}
