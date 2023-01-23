
#' Estimate the ATE using Random-Intercept, Random-Coefficient (RIRC) Model
#'
#' @inheritParams estimate_ATE_FIRC
#' @rdname estimate_ATE_RIRC
#'@param data Dataframe with all needed variables.
#' @export
estimate_ATE_RIRC <- function(Yobs, Z, B, data = NULL,
                              include_testing = FALSE,  REML = !include_testing,
                              pool = FALSE, control_formula = NULL) {
  if (!is.null( control_formula)) {
    stopifnot(!is.null(data))
    stopifnot(!missing("Yobs"))
  }
    if ( include_testing && REML ) {
        stop( "Cannot do likelihood ratio test when REML=TRUE" )
    }

  if (!is.null(data)) {
    if (missing("Yobs")) {
      data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]])
    } else {
      d2 <- data
      d2$Yobs <- eval(substitute(Yobs), data)
      d2$Z <- eval(substitute(Z), data)
      d2$B <- eval(substitute(B), data)
      data <- d2
      rm(d2)
    }
  } else {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
  }
  stopifnot( length( unique( data$Z ) ) == 2 )
  stopifnot( is.numeric( data$Yobs ) )

  data$B = droplevels( as.factor(data$B) )


  #fit multilevel model and extract tau
  method <- ifelse( REML, "REML", "ML" )

  formula <- make_base_formula(control_formula = control_formula, data = data)
  if (pool) {
    re.mod <- nlme::lme(formula, data = data, random = ~ 1 + Z | B,
                        na.action = stats::na.exclude,
                        method = method,
                        control = nlme::lmeControl(opt="optim",
                                                   returnObject = TRUE))
  } else {
    re.mod <- nlme::lme(formula, data = data, random = ~ 1 + Z | B,
                        weights = nlme::varIdent(form = ~ 1 | Z),
                        na.action = stats::na.exclude,
                        method = method,
                        control = nlme::lmeControl(opt="optim",
                                                   returnObject = TRUE))
  }

  if (include_testing) {
    # Test for cross site variation (???)
    if (pool) {
      re.mod.null <- nlme::lme(formula, data = data, random = ~ 1 | B, na.action = stats::na.exclude, method = method,
        control = nlme::lmeControl(opt = "optim", returnObject = TRUE))
    } else {
      re.mod.null <- nlme::lme(formula, data = data, random = ~ 1 | B, weights = nlme::varIdent(form = ~ 1 | Z), na.action = stats::na.exclude,
        method = method, control = nlme::lmeControl(opt="optim", returnObject = TRUE))
    }
    # M0.null <- lm( Yobs ~ 0 + B + Z, data=data )
    td <- as.numeric(deviance(re.mod.null) - deviance(re.mod))
    p_variation <- 0.5 * pchisq(td, 2, lower.tail = FALSE) + 0.5 * pchisq(td, 1, lower.tail = FALSE)
  } else {
    p_variation <- NA
    td <- NA
  }

  # get ATE and SE
  ATE <- nlme::fixef(re.mod)[[2]]
  SE_ATE <- sqrt(vcov(re.mod)[2, 2])

  #cextract tau
  vc <- nlme::VarCorr(re.mod)
  suppressWarnings(storage.mode(vc) <- "numeric")
  tau_hat <- vc["Z", "StdDev"]

  res <- list(ATE_hat = ATE, SE_ATE = SE_ATE,
              tau_hat = tau_hat, SE_tau = NA,
              p_variation = p_variation, deviance = td)
  class( res ) = "multisiteresult"
  attr( res, "args" ) = list(  model = "RIRC",
                               method = method )

  return( res )

}
