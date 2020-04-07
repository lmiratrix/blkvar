#' Compare different MLM methods
#'
#' @param Yobs vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column name in data)
#' @param B vector of block ids (or column name in data)
#' @param siteID site ids (variable name as string if data frame passed) (if randomization blocks are nested in site).
#' @param data frame holding Y, Z, B and (possibly a column with name specified by siteID).
#' @param control_formula What variables to control for, in the form of "~ X1 + X2".
#'
#' @importFrom stats aggregate lm quantile rnorm sd var as.formula cor cov filter model.matrix na.exclude pf pnorm predict qchisq qf qnorm rbinom reshape resid runif vcov weighted.mean
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
  FIRC <- estimate_ATE_FIRC(Yobs, Z, B, data = data, siteID = siteID, REML = TRUE, include_testing = FALSE, control_formula = control_formula )
  RIRC <- estimate_ATE_RIRC( Yobs, Z, B, data=data, REML = TRUE, include_testing = FALSE, control_formula = control_formula )
  mlms <- data.frame(method=c("RICC", "FIRC", "RIRC"), tau = c( RICC$ATE, FIRC$ATE, RIRC$ATE ), SE = c( RICC$SE.ATE, FIRC$SE.ATE, RIRC$SE.ATE ),
    stringsAsFactors = FALSE)
  if (!is.null(control_formula)) {
    mlms$method <- paste0(mlms$method, "-adj")
  }
  mlms
}




# #### Testing and demo of this code ####

# if  (FALSE ) {
    # dat = generate_blocked_data_obs( n_k = 4:10, p = 0.2 )
    # dat
    # table( dat$Z, dat$B )
    # compare_methods( data = dat[ c("Yobs", "Z","B" ) ] )
# }
