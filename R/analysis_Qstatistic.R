


#' Calculate confidence interval for cross site variation using Q-statistic test
#' inversion.
#'
#' Two versions of this, one that takes raw data (in which case a linear model
#' is fit to it to get individual site impact estimates), and a second, '_stat',
#' that takes precalculated point estimates and standard errors.
#'
#' The first simply calculates ATE_hat and SE_hat, and then calls the second.
#'
#' Code taken from prior work of Catherine Armstrong.
#'
#' @param ATE_hat List of estimated ATEs, one for each site.
#' @param SE_hat List of associated estimated SEs, assumed to be known.
#' @param alpha The level of the test.  The CI will be a 1-2alpha confidence
#'   interval.
#' @param calc_CI Logical
#' @return List of several estimated quantities: cross-site variation, p-value
#'   for any variation, whether the p-value is less than passed alpha, the Q
#'   statistic itself, and the confidence interval range.
#' @export
#' @rdname analysis_Qstatistic
analysis_Qstatistic_stat <- function( ATE_hat, SE_hat, alpha = 0.05, calc_CI = FALSE) {

  bj = ATE_hat
  vj = SE_hat^2
  wj <- 1 / vj
  s <- length(ATE_hat)

  bbar <- sum(wj * bj) / (sum(wj))
  q <- sum((bj - bbar) ^ 2 / vj)
  pval <- pchisq(q, df = (length(bj) - 1), lower.tail = FALSE)
  reject <- (pval < alpha)

  if (calc_CI) {
  	## get confidence interval
    ## confidence interval
    #test tau values from 0 to 5 (what is a good way to set this upper bound?)
    tau_test <- seq(0,5,.01) ##come back and make it increment by 0.01
    ## test two-sided version
    lowbound <- qchisq(alpha, s - 1)
    highbound <- qchisq(1 - alpha, df = (s - 1))
    ## intialize values
    q_invert <- c()
    CI_95 <- c()
    for (i in 1:length(tau_test)) {
      #calculate new denominator that takes into account tau^2
      denom <- vj + tau_test[i] ^ 2
      #new q-stat
      q_invert[i] <- sum((bj - bbar) ^ 2 / denom)
      #compare q.stat to lower and upperbounds (Weiss et al, JREE p. 55)
      CI_95[i] <- (q_invert[i] >= lowbound & q_invert[i] <= highbound)
      # Also calculate p-values for Hodges-Lehman estimate
    }

    #extract bounds
    #special condition to explore when tau_test[CI_90 == 1]
    pval.hi <- pchisq(q_invert, df = s - 1, lower.tail = FALSE)
    pval.low <- pchisq(q_invert, df = s - 1, lower.tail = TRUE)
    pvals <- 2 * pmin(pval.hi, pval.low)
    mx <- which.max(pvals)
    tau_hat <- tau_test[[mx]]

    if (length(tau_test[CI_95 == 1]) == 0) {
      CI_low <- NA
      CI_high <- NA
    } else {
      CI_high <- max(tau_test[CI_95 == 1])
      CI_low <- min(tau_test[CI_95 == 1])
    }

  } else { # we didn't bother calculating the Conf Int.
    tau_hat = NA
    CI_low = NULL
    CI_high = NULL
  }

  return(list(tau_hat = tau_hat, p.value = pval, reject = reject, Q = q, CI_low = CI_low, CI_high = CI_high))
}



#' @param data Dataframe with outcome, treatment, and blocking factor.
#' @param B (site id)
#' @param Yobs (outcome)
#' @param Z (binary treatment 0/1)
#' @param siteID If not null, name of siteID that has randomization blocks
#' @param data frame holding Y, Z, B and (possibly a column with name specified
#'   by siteID).
#' @export
#' @rdname analysis_Qstatistic
analysis_Qstatistic <- function(Yobs, Z, B, siteID = NULL, data = NULL, alpha = 0.05, calc_CI = FALSE) {
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
        if (!is.null(siteID)) {
            data$siteID <- siteID
        }
    }

    # make sites RA blocks if there are no sites.
    if (is.null(siteID)) {
        data$siteID = data$B
    }

    n <- nrow(data)
    s <- length(unique(data$B))
    s.site <- length(unique(data$siteID))
    # calculate Q-statistic
    # run ols model with no intercept and no "treatment intercept"
    ols <- nlme::gls(Yobs ~ 0 + Z + factor(B) + factor(siteID):Z - Z,
                     data = data,
                     weights = nlme::varIdent(form = ~ 1 | Z),
                     na.action = stats::na.exclude,
                     control = nlme::lmeControl(opt = "optim", returnObject = TRUE))
    ATE_hat <- ols$coefficients[(s + 1):(s + s.site)]
    SE_hat <- (coef(summary(ols))[(s + 1):(s + s.site), 2])

    analysis_Qstatistic_stat( ATE_hat, SE_hat, alpha=alpha, calc_CI = calc_CI)
}
