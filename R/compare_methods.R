
#' Compare different MLM methods
#'
#' @param Yobs vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column name in data)
#' @param B vector of block ids (or column name in data)
#' @param siteID site ids (variable name as string if data frame passed) (if randomization blocks are nested in site).
#' @param data frame holding Y, Z, B and (possibly a column with name specified by siteID).
#' @param control_formula What variables to control for, in the form of "~ X1 + X2".
#'
#' @export
compare_MLM_methods <- function( Yobs, Z, B, siteID = NULL, data = NULL, control_formula = NULL ) {
    if (!is.null(data)) {
        if (missing( "Yobs")) {
            data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]])
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
            siteID = "siteID"
        }
    }
    stopifnot(length(unique(data$Z)) == 2)
    stopifnot(is.numeric(data$Yobs))
    RICC <- estimate_ATE_RICC(Yobs, Z, B, data = data, REML = TRUE, control_formula = control_formula )
    FIRC <- estimate_ATE_FIRC(Yobs, Z, B, data = data, siteID = siteID, REML = TRUE,
                              include_testing = FALSE, control_formula = control_formula )
    RIRC <- estimate_ATE_RIRC( Yobs, Z, B, data=data, REML = TRUE, include_testing = FALSE,
                               control_formula = control_formula )
    mlms <- data.frame(method=c("RICC", "FIRC", "RIRC"),
                       ATE_hat = c( RICC$ATE_hat, FIRC$ATE_hat, RIRC$ATE_hat ),
                       SE = c( RICC$SE_ATE, FIRC$SE_ATE, RIRC$SE_ATE ),
                       stringsAsFactors = FALSE)
    if (!is.null(control_formula)) {
        mlms$method <- paste0(mlms$method, "-adj")
    }
    mlms
}




# Testing and demo of this code #

# if  (FALSE ) {
# dat = generate_blocked_data_obs( n_k = 4:10, p = 0.2 )
# dat
# table( dat$Z, dat$B )
# compare_methods( data = dat[ c("Yobs", "Z","B" ) ] )
# }




#' Block variance method comparison function
#'
#' This function calculates the point estimates and SE estimates for a variety
#' of the blocked designs.
#'
#' @param Yobs vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column name in
#'   data)
#' @param B vector of block ids (or column name in data)
#' @param siteID site ids (variable name as string if data frame passed) (if
#'   randomization blocks are nested in site).
#' @param data frame holding Y, Z, B and (possibly a column with name specified
#'   by siteID).
#' @param include_MLM Include MLM estimators
#' @param include_DB Include Design-Based estimators (taken from RCTYes
#'   documentation and prior literature).
#' @param include_LM Include Linear Model-based estimators (including
#'   Huber-White SEs, etc.)
#' @param include_DBBlended Include DB estimators applied to small block (blocks
#'   with singleton treatment or control units) and classic Neyman to large
#'   blocks.
#' @param weight_LM_method Passed to weight.method of
#'   weighted_linear_estimators; specifies the weighting method.
#' @param weight_LM_scale_weights This argument is passed to scale.weights of
#'   weighted_linear_estimators
#' @param control_formula What variables to control for, in the form of "~ X1 +
#'   X2".
#' @param include_block Include the standard error estimation found in Pashley
#'   (in particular, these work for blocks with singleton units).
#' @param include_method_characteristics Include details of the methods (target
#'   estimands and sampling framework assumed) in the return value.
#'
#' @return Dataframe of point estimates and standard errors for each method
#'   considered. If \code{include_method_characteristics=TRUE} also add some
#'   features of the methods as additional columns.
#' @importFrom stats aggregate lm quantile rnorm sd var as.formula cor cov
#'   filter model.matrix na.exclude pf pnorm predict qchisq qf qnorm rbinom
#'   reshape resid runif vcov weighted.mean
#' @export
compare_methods <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                            include_block = TRUE, include_MLM = TRUE,
                            include_DB = TRUE,
                            include_LM = TRUE, include_DBBlended = FALSE,
                            include_method_characteristics = FALSE,
                            weight_LM_method = "survey",
                            weight_LM_scale_weights = TRUE, control_formula = NULL) {

  if (!is.null(control_formula)) {
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
       methods_list <- c( methods_list,
                          "rct_yes_atll", "rct_yes_small",
                          "rct_yes_mod_all", "rct_yes_mod_small")
     }
     fits <- sapply( methods_list, function(m) {
       dd <- block_estimator_tabulated(summary_stats = summary_stats, method = m, throw.warnings = FALSE)
       c(dd$ATE_hat, dd$se_est)
      })
      SE_estimates <- fits[2, ]
      ATE_estimates <- fits[1, ]
      summary_table <- data.frame(method = methods_list, ATE_hat = ATE_estimates,
                                  SE = SE_estimates, stringsAsFactors = FALSE)
  } else {
    summary_table <- data.frame()
  }

  # stash canonical name for site ID, if we have one.
  if (!is.null(siteID)) {
    siteID <- "siteID"
  }

  if (include_DB) {
      if (is.null(control_formula)) {
          # Design based methods
          DB.fi <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID,
                                                         method="finite", weight="individual" )
          DB.fs <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID,
                                                         method="finite", weight="site" )
          DB.si <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID,
                                                         method="superpop", weight="individual" )
          DB.ss <- estimate_ATE_design_based_from_stats( summary_stats, siteID=siteID,
                                                         method="superpop", weight="site" )
          DB <- dplyr::bind_rows( DB.fi, DB.fs, DB.si, DB.ss )
          DB$method <- c( "DB-FP-Persons", "DB-FP-Sites", "DB-SP-Persons", "DB-SP-Sites")
      } else {
          # Design based methods with covariate adjustment
          DB.fi <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID,
                                                      method = "finite", weight = "individual",
                                                      control_formula = control_formula )
          DB.fs <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID,
                                                      method = "finite", weight = "site",
                                                      control_formula = control_formula )
          DB.si <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID,
                                                      method = "superpop", weight = "individual",
                                                      control_formula = control_formula)
          DB.ss <- estimate_ATE_design_based_adjusted(Yobs ~ Z * B, data = data, siteID = siteID,
                                                      method = "superpop", weight = "site",
                                                      control_formula = control_formula )
          DB <- dplyr::bind_rows(DB.fi, DB.fs, DB.si, DB.ss)
          DB$method <- c( "DB-FP-Persons-adj", "DB-FP-Sites-adj", "DB-SP-Persons-adj", "DB-SP-Sites-adj")
      }
      DB$weight <- NULL
      names(DB)[1] <- "ATE_hat"
      summary_table <- dplyr::bind_rows(summary_table, DB)
      # ensure canonical ordering of method name first
      summary_table = dplyr::select( summary_table, method, everything() )
  }

  if (include_LM) {
    lms <- linear_model_estimators(Yobs, Z, B, data = data, siteID = siteID,
                                   block.stats = summary_stats,
                                   control_formula = control_formula,
                                   weight_LM_method = weight_LM_method,
      weight_LM_scale_weights = weight_LM_scale_weights)
    summary_table = dplyr::bind_rows(summary_table, lms)
  }
  if (include_MLM) {
    mlms <- compare_MLM_methods(Yobs, Z, B, data = data, siteID = siteID,
                                control_formula = control_formula)
    summary_table <- dplyr::bind_rows( summary_table, mlms )
  }

  # Add info on the methods (e.g., what estimand they are targeting)
  if (include_method_characteristics) {
    mc <- method_characteristics()
    #mcm = mc$method
    #names(mcm) = mc$fullname
    #summary_table$method = mcm[ as.character( summary_table$method ) ]
    if (!is.null(control_formula)) {
      mc$method <- paste0(mc$method, "-adj")
    }
    summary_table <- merge( summary_table, mc, by = "method",
                            all.x = TRUE, all.y = FALSE)
  }
  if ( nrow( summary_table ) > 0 ) {
    summary_table = tibble::remove_rownames( summary_table ) %>%
      arrange( method )
  }

  return(summary_table)
}
