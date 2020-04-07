#' Block variance method comparison function
#'
#' This function calculates the point estimates and SE estimates for a variety
#' of the blocked designs.
#'
#' @param Yobs vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column name in data)
#' @param B vector of block ids (or column name in data)
#' @param siteID site ids (variable name as string if data frame passed) (if randomization blocks are nested in site).
#' @param data frame holding Y, Z, B and (possibly a column with name specified by siteID).
#' @param include_MLM Include MLM estimators
#' @param include_DB Include Design-Based estimators (taken from RCTYes documentation and prior literature).
#' @param include_LM Include Linear Model-based estimators (including
#'   Huber-White SEs, etc.)
#' @param include_DBBlended Include DB estimators applied to small block
#'   and classic Neyman to large blocks.
#' @param weight.LM.method Argument passed to weight.method of weighted_linear_estimators
#' @param weight.LM.scale.weights Argument passed to sclae.weights of weighted_linear_estimators
#' @param control.formula What variables to control for, in the form of "~ X1 + X2".
#' @param include_block Include the Pashley blocking variants.
#' @param include_method_characteristics Include details of the methods (target estimands and sampling framework assumed).
#'
#' @importFrom stats aggregate lm quantile rnorm sd var as.formula cor cov filter model.matrix na.exclude pf pnorm predict qchisq qf qnorm rbinom reshape resid runif vcov weighted.mean
#' @export
compare_methods <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                            include_block = TRUE, include_MLM = TRUE, include_DB = TRUE,
                            include_LM = TRUE, include_DBBlended = FALSE,
                            include_method_characteristics = FALSE,
                            weight.LM.method = "survey",
                            weight.LM.scale.weights = TRUE, control.formula = NULL) {

  if (!is.null(control.formula)) {
    stopifnot( !is.null(data))
    # stopifnot( !missing("Yobs"))
  }

  # This code block takes the parameters of
  # Yobs, Z, B, siteID = NULL, data=NULL, ...
  # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
  if (!is.null(data)) {
    if (missing("Yobs")) {
      data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]])
      n.tx.lvls <- length(unique(data$Z))
      stopifnot(n.tx.lvls == 2)
      stopifnot(is.numeric(data$Yobs))
    } else {
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
    if (!is.null( siteID)) {
      data$siteID <- siteID
    }
  }
  n <- nrow( data )
  #Quick check that input is correct
  if (is.numeric(data$Z) == FALSE ) {
    stop("Treatment indicator should be vector of ones and zeros")
  }
  if ((sum(data$Z == 1) + sum(data$Z == 0)) != n) {
    stop("Treatment indicator should be vector of ones and zeros")
  }

  #Get data into table
  summary_stats <- calc_summary_stats(Yobs, Z, B, data = data, siteID= siteID, add.neyman = TRUE)

  if (include_block || include_DBBlended) {
    method_list <- c()
    if (include_block) {
      methods_list <- c("hybrid_m", "hybrid_p", "plug_in_big")
     }
     if (include_DBBlended) {
       methods_list <- c( methods_list, "rct_yes_atll", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small")
     }
     fits <- sapply( methods_list, function(m) {
       dd <- block_estimator_tabulated(summary_stats = summary_stats, method = m, throw.warnings = FALSE)
       c(dd$tau_est, dd$se_est)
      })
      SE_estimates <- fits[2, ]
      tau_estimates <- fits[1, ]
      summary_table <- data.frame(method = methods_list, tau = tau_estimates, SE = SE_estimates, stringsAsFactors = FALSE)
  } else {
    summary_table <- data.frame()
  }

  # stash canonical name for site ID, if we have one.
  if (!is.null(siteID)) {
    siteID <- "siteID"
  }

  if (include_DB) {
    if (is.null(control.formula)) {
      # Design based methods
      DB.fi <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID, method="finite", weight="individual" )
      DB.fs <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID, method="finite", weight="site" )
      DB.si <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID, method="superpop", weight="individual" )
      DB.ss <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID, method="superpop", weight="site" )
      DB <- dplyr::bind_rows( DB.fi, DB.fs, DB.si, DB.ss )
      DB$method <- c( "DB-FP-Persons", "DB-FP-Sites", "DB-SP-Persons", "DB-SP-Sites")
      DB$weight <- NULL
      names(DB)[1] <- "tau"
      summary_table <- dplyr::bind_rows(summary_table, DB)
    } else {
      # Design based methods with covariate adjustment
      DB.fi <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID, method = "finite", weight = "individual", control.formula = control.formula )
      DB.fs <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID, method = "finite", weight = "site", control.formula = control.formula )
      DB.si <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID, method = "superpop", weight = "individual", control.formula = control.formula)
      DB.ss <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID, method = "superpop", weight = "site", control.formula = control.formula )
      DB <- dplyr::bind_rows(DB.fi, DB.fs, DB.si, DB.ss)
      DB$method <- c( "DB-FP-Persons-adj", "DB-FP-Sites-adj", "DB-SP-Persons-adj", "DB-SP-Sites-adj")
      DB$weight <- NULL
      names(DB)[1] <- "tau"
      summary_table <- dplyr::bind_rows(summary_table, DB)
    }
  }

  if (include_LM) {
    lms <- linear_model_estimators(Yobs, Z, B, data = data, siteID = siteID, block.stats = summary_stats, control.formula = control.formula, weight.LM.method = weight.LM.method,
      weight.LM.scale.weights = weight.LM.scale.weights)
    summary_table = dplyr::bind_rows(summary_table, lms)
  }
  if (include_MLM) {
    mlms <- compare_MLM_methods(Yobs, Z, B, data = data, siteID = siteID, control.formula = control.formula)
    summary_table <- dplyr::bind_rows( summary_table, mlms )
  }

  # Add info on the methods (e.g., what estimand they are targeting)
  if (include_method_characteristics) {
    mc <- method_characteristics()
    #mcm = mc$method
    #names(mcm) = mc$fullname
    #summary_table$method = mcm[ as.character( summary_table$method ) ]
    if (!is.null(control.formula)) {
      mc$method <- paste0(mc$method, "-adj")
    }
    summary_table <- merge( summary_table, mc, by = "method", all.x = TRUE, all.y = FALSE)
  }
  return(summary_table)
}
