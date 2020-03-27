#' Estimate the ATE using Random-Intercept, Constant-Coefficient (RICC) Model.
#'
#' This model has a single treatment coefficient, and a random intercept for the
#' site control average. So it is analogous to a fixed effect model, but with a
#' random effect.
#'
#' There is no test for cross site variation for this method, since we assume none.
#'
#' @param B Name of the block indicator.
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z vector of assignment indicators (1==treated)
#' @param REML TOADD
#' @param control.formula The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)
#' @param data Dataframe with all needed variables.
#' @rdname estimate_ATE_RIRC
#'
#' @export
estimate_ATE_RICC <- function(Yobs, Z, B, data = NULL, REML = FALSE, control.formula = NULL) {
  if (!is.null(control.formula)) {
    stopifnot(!is.null(data))
    stopifnot(!missing("Yobs"))
  }
  if (!is.null(data)) {
    if (missing( "Yobs")) {
      data <- data.frame( Yobs = data[[1]], Z = data[[2]], B = data[[3]])
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
  stopifnot(length(unique(data$Z)) == 2)
  stopifnot(is.numeric(data$Yobs))
  
  #fit multilevel model and extract tau
  method <- ifelse( REML, "REML", "ML")
  formula <- make_base_formula(control.formula = control.formula, data = data)
  
  re.mod <- nlme::lme(formula, data = data, random = ~ 1 | B, weights = nlme::varIdent(form = ~ 1 | Z), na.action = stats::na.exclude, method = method,
    control = nlme::lmeControl(opt="optim", returnObject = TRUE))
    
  # get ATE and SE
  ATE <- nlme::fixef(re.mod)[[2]]
  SE.ATE <- sqrt(vcov(re.mod)[2, 2])
  return(list(ATE = ATE, SE.ATE = SE.ATE, tau.hat = NA, SE.tau = NA, p.variation = NA, deviance = NA))
}