## ---------- Both pooled and unpooled FIRC method ------------

#' Estimate the ATE using Random-Intercept, Random-Coefficient (RIRC)
#' Model
#'
#' This method fits the FIRC model to estimate (1) ATE across sites
#' and (2) cross site treatment variation.
#'
#' Warning: If you want to test for cross-site variation you need REML
#' = FALSE to allow for likelihood ratio testing.  If you want to
#' estimate this variation, REML = TRUE is recommended as the ML
#' estimate tends to be wildly biased towards 0 when modest amounts of
#' cross-site variation are present.
#'
#' Acknowledgement: Unpooled version taken and adapted from
#' Catherine's weiss.tau() method.
#'
#'
#' @param anova Use the anova() method to do the test for significance
#'   between the models.  FALSE means do the modified chi-squared
#'   test.
#' @param pool  TRUE means tx and co have same reBual variance. FALSE
#'   gives seperate estimates for each (recommended, default).
#' @param B Name of the block indicator.
#' @param siteID Character name of the ID variable of site (blocks are
#'   conBered nested in site).  If omitted, then blocks are considered
#'   sites (the default).
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z vector of assignment indicators (1==treated)
#' @param REML Logical, Restricted maximum likelihood or maximum
#'   likelihood estimation.  Default of TRUE if include_testing =
#'   FALSE and FALSE otherwise. Note that FIRC can fail badly in
#'   estimating cross-site variation when REML=FALSE.
#' @param siteID If not null, name of siteID that has randomization
#'   blocks
#' @param control_formula The control_formula argument must be of the
#'   form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)
#' @param include_testing Logical Include likelihood ratio test for
#'   cross-site treatment variation.
#' @param data Dataframe with all needed variables.
#' @param keep_EB_estimates TRUE means returned object has EB
#'   estimates.  FALSE means do not keep them.
#' @export
estimate_ATE_FIRC <- function(Yobs, Z, B, siteID = NULL,
                              control_formula = NULL, data = NULL,
                              include_testing = FALSE,  REML = !include_testing,
                              anova = FALSE, pool = FALSE,
                              keep_EB_estimates = TRUE) {

  if ( include_testing && REML ) {
      stop( "Cannot do likelihood ratio test when REML=TRUE" )
  }

  if (!is.null(control_formula)) {
    stopifnot(!is.null(data))
    stopifnot(!missing( "Yobs"))
  }

  # This code block takes the parameters of
  # Yobs, Z, B, siteID = NULL, data=NULL, ...
  # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
  if (!is.null(data)) {
    if (missing( "Yobs")) {
  	  data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]])
  	  n.tx.lvls <- length(unique(data$Z))
  	  stopifnot(n.tx.lvls == 2)
  	  stopifnot(is.numeric(data$Yobs))
    } else {
      Yobs.n <- as.character(substitute(Yobs))
      Z.n <- as.character(substitute(Z))
      B.n <- as.character(substitute(B))
      data <- dplyr::rename(data,
                            Yobs = dplyr::all_of(Yobs.n),
                            Z = dplyr::all_of(Z.n),
                            B = dplyr::all_of(B.n) )
      if (!is.null(siteID)) {
        data$siteID <- data[[siteID]]
        stopifnot(!is.null(data$siteID))
      }
    }
  } else {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
    if (!is.null( siteID)) {
      data$siteID <- siteID
    }
  }

  if (is.null(siteID)) {
    data$siteID <- data$B
  }

  stopifnot( length(unique(data$Z)) == 2)
  stopifnot(is.numeric(data$Yobs))

  # ensure no weird numerical
  data$B = droplevels( as.factor(data$B) )


  #fit multilevel model and extract tau
  method <- ifelse( REML, "REML", "ML" )

  # Make control variable function
  formula <- make_FE_formula("Yobs", "Z", "B", control_formula = control_formula, data = data)

  if (pool) {
    re.mod <- nlme::lme(formula, data = data, random = ~ 0 + Z | siteID, method = method,
                        control=nlme::lmeControl(opt="optim", returnObject = TRUE))
  } else {
    re.mod <- nlme::lme(formula, data = data, random = ~ 0 + Z | siteID, method = method,
                        weights = nlme::varIdent(form = ~ 1 | Z),
                        na.action=stats::na.exclude,
                        control=nlme::lmeControl(opt="optim",returnObject=TRUE))
  }

  if (include_testing) {
    # Test for cross site variation
    if (pool) {
      re.mod.null <- nlme::gls(formula, data = data, method = method, control=nlme::lmeControl(opt="optim", returnObject = TRUE))
    } else {
      re.mod.null <- nlme::gls(formula, data = data, weights = nlme::varIdent(form = ~ 1 | Z), na.action = stats::na.exclude, method = method,
        control = nlme::lmeControl(opt = "optim", returnObject = TRUE))
    }
    # M0.null = lm( Yobs ~ 0 + B + Z, data=data )
    if (anova) {
      stopifnot(REML == FALSE)
      # This tosses a warning that the linear model is being swapped into a MLM
      # "original model was of class "gls", updated model is of class "lme""
      myanova <- suppressWarnings( lmtest::lrtest( re.mod.null, re.mod ) ) #anova(re.mod.null, re.mod)
      p.value.anova <- myanova[2,5] #myanova[2,9]
      p_variation <- (p.value.anova / 2)  # divide by 2 by same logic as chi-squared test.
      td <- NA
    } else {
      td <- as.numeric(deviance(re.mod.null) - deviance(re.mod))
      p_variation <- 0.5 * pchisq(td, 1, lower.tail = FALSE )
    }
  } else {
    p_variation <- td <- NA
  }

  # get ATE and SE
  ATE_hat <- nlme::fixef(re.mod)[[1]]
  SE_ATE <- sqrt(vcov(re.mod)[1, 1])

  # extract tau (estimated cross site variation)
  vc <- nlme::VarCorr(re.mod)
  suppressWarnings(storage.mode(vc) <- "numeric")
  tau_hat <- vc["Z","StdDev"]

  # extract se of tau
  ## not source code: https://stackoverflow.com/questions/31694812/standard-error-of-variance-component-from-the-output-of-lmer/31704646#31704646
  ## soure code: https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-standard-errors-for-variance-components-from-mixed-models/

  # let logged sds first
  var <- re.mod$apVar
  if (is.character(var)) {
    SE_tau <- NA
  } else {
    par <- attr(var, "Pars")
    #transform to variances
    var.re <- exp(par) ^ 2
    #use delta method from msm package
    SE_tau <- msm::deltamethod(~ exp(x1) ^ 2, par, var)
  }
  res <- list(ATE_hat = ATE_hat, SE_ATE = SE_ATE,
              tau_hat = tau_hat, SE_tau=SE_tau,
              p_variation = p_variation, deviance = td )
  class( res ) = "multisiteresult"
  attr( res, "args" ) = list(  model = "FIRC",
                               method = method )

  if ( keep_EB_estimates ) {
      site_ates = tibble( sid = rownames( coef(re.mod) ),
                          beta_hat = coef( re.mod )$Z )

      res$eb_estimates = site_ates
  }

  return( res )
}


#' Get empirical Bayes estimates from fit FIRC model
#'
#' @param x A FIRC model multisiteresult object.
#' @return Tibble with site IDs and estimated site average impacts.
#' @export
get_EB_estimates = function( x ) {
    stopifnot( is.multisiteresult(x) )
    return( x$eb_estimates )
}







