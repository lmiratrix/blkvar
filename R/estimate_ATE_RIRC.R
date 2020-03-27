#' Estimate the ATE using Random-Intercept, Random-Coefficient (RIRC) Models
#'
#' @inheritParams estimate_ATE_FIRC
#' @rdname estimate_ATE_RIRC
#' @param pool  TRUE means tx and co have same reBual variance. FALSE gives seperate estimates for each (recommended, default).
#' @param B Name of the block indicator.
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z vector of assignment indicators (1==treated)
#' @param REML TOADD
# '@param control.formula The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)
#'@param data Dataframe with all needed variables.
#' @export

estimate_ATE_RIRC <- function(Yobs, Z, B, data = NULL, REML = FALSE, include.testing = TRUE, pool = FALSE, control.formula = NULL) {
  stopifnot(!(include.testing && REML))
  if (!is.null( control.formula)) {
    stopifnot(!is.null(data))
    stopifnot(!missing("Yobs"))
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

  #fit multilevel model and extract tau
  method <- ifelse( REML, "REML", "ML" )

  formula <- make_base_formula(control.formula = control.formula, data = data)
  if (pool) {
    re.mod <- nlme::lme(formula, data = data, random = ~ 1 + Z | B, na.action = stats::na.exclude, method = method,
      control = nlme::lmeControl(opt="optim", returnObject = TRUE))
  } else {
    re.mod <- nlme::lme(formula, data = data, random = ~ 1 + Z | B, weights = nlme::varIdent(form = ~ 1 | Z), na.action = stats::na.exclude,
      method = method, control = nlme::lmeControl(opt="optim", returnObject = TRUE))
  }

  if (include.testing) {
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
    p.variation <- 0.5 * pchisq(td, 2, lower.tail = FALSE) + 0.5 * pchisq(td, 1, lower.tail = FALSE)
  } else {
    p.variation <- NA
    td <- NA
  }

  # get ATE and SE
  ATE <- nlme::fixef(re.mod)[[2]]
  SE.ATE <- sqrt(vcov(re.mod)[2, 2])

  #cextract tau
  vc <- nlme::VarCorr(re.mod)
  suppressWarnings(storage.mode(vc) <- "numeric")
  tau.hat <- vc["Z", "StdDev"]
  
  return(list(ATE = ATE, SE.ATE = SE.ATE, tau.hat = tau.hat, SE.tau = NA, p.variation = p.variation, deviance = td))
}
