## ---------------- unpooled FIRC: accounting for covariates -----------------
##

#' Covariate adjusted test for cross-site variation
#'
#' This method fits a FIRC model but also includes a site-level covariate, X,
#' potentially predictive of treatment impact.
#'
#' This fits unpooled FIRC models.
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data Dataframe with defined Yobs, Z, and B variables.
#' @param siteID If not null, name of siteID that has randomization blocks
# '@param control.formula The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)
#' @param X Site-level covariate ideally predictive of treatment variation.
#' @param anova Use the anova() method to do the test for significance between the models. FALSE means do the modified chi-squared test.
#' @param REML TOADD
#'
analysis_FIRC_cov <- function(Yobs, Z, B, X, siteID = NULL, data = NULL, REML = FALSE, anova = FALSE) {
  # stopifnot( !( include.testing && REML ) )
  # This code block takes the parameters of
  # Yobs, Z, B, siteID = NULL, data=NULL, ...
  # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
  if (!is.null(data)) {
    if (missing("Yobs")) {
      data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]], X = data[[4]])
      n.tx.lvls <- length(unique(data$Z))
      stopifnot(n.tx.lvls == 2)
      stopifnot(is.numeric(data$Yobs))
    } else {
      if (!is.null(siteID)) {
        siteID <- data[[siteID]]
      }
      data <- data.frame(Yobs = eval(substitute(Yobs), data), Z = eval(substitute(Z), data), B = eval(substitute(B), data), X = eval(substitute(X), data))
      data$siteID <- siteID
    }
  } else {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B, X = X)
    if (!is.null(siteID)) {
      data$siteID = siteID
    }
  }
  if (is.null(data$siteID)) {
    data$siteID <- data$B
  }
  
  stopifnot(length(unique(data$Z)) == 2)
  stopifnot(is.numeric(data$Yobs))

  # fit multilevel model and extract tau
  method <- ifelse(REML, "REML", "ML")
  # fit multilevel model and extract pvalue
  re.mod <- nlme::lme(Yobs ~ 0 + Z + Z:X + B, data = data, random = ~ 0 + Z | B, weights = nlme::varIdent(form = ~ 1 | Z), na.action = stats::na.exclude,
    method = "ML", control = nlme::lmeControl(opt = "optim", returnObject = TRUE))
  
  # Test for cross site variation
  re.mod.null <- nlme::gls(Yobs ~ 0 + Z + Z:X + B, data = data, weights = nlme::varIdent(form = ~ 1 | Z), na.action = stats::na.exclude,
    method = "ML", control = nlme::lmeControl(opt = "optim", returnObject = TRUE))

  # Generate the different tests we can have for cross site variation in this model

  # Test 1: LR test (combination)
  if (anova) {
    stopifnot(REML == FALSE)
    myanova <- anova(re.mod.null, re.mod)
    # TODO: Fix this!
    p.value.comb <- myanova[2, 9] / 2  # divide by 2 by same logic as chi-squared test.
  } else {
    td <- abs(as.numeric(deviance(re.mod) - deviance(re.mod.null)))
    p.value.comb <- 0.5 * pchisq(td, 2, lower.tail = FALSE) + 0.5 * pchisq(td, 1, lower.tail = FALSE)
  }

  # Test 2: Test the systematic coefficient
  SEs <- sqrt(vcov(re.mod )["Z:X","Z:X"])
  tstat <- nlme::fixef(re.mod )[["Z:X"]] / SEs
  p.value.sys <- 2 * pnorm( - abs( tstat ))

  data.frame(method = c("FIRC combination", "FIRC systematic.RTx" ), pvalue = c( p.value.comb, p.value.sys))
}
