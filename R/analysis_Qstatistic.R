


scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



find_p_value = function( tau, resid2, vj ) {
    if ( length( tau ) > 1 ) {
        map_dbl( tau, find_p_value, resid2=resid2, vj=vj )
    }
    s = length( resid2 )

    denom <- vj + tau^2
    q = sum( resid2 / denom )
    pval.hi <- pchisq(q, df = s - 1, lower.tail = FALSE)
    pval.low <- pchisq(q, df = s - 1, lower.tail = TRUE)
    pvals <- 2 * pmin(pval.hi, pval.low)
    pvals
}


find_CI_endpoint = function( resid2, vj, threshold, tau_start, tau,
                             tolerance = 0.01,
                             bracketed = FALSE ) {

    qmax = sum( resid2 / vj )
    if ( qmax < threshold ) {
        return( list( tau=NA, Q=qmax ) )
    }

    stopifnot( tau_start < tau )

    if ( ! bracketed ) {
        denom <- vj + tau ^ 2
        q = sum( resid2 / denom )

        while( q > threshold ) {
            #scat( "tau = %.2f\n", tau )
            tau_start = tau
            tau = tau * 2
            denom <- vj + tau ^ 2
            q = sum( resid2 / denom )
        }
        #scat( "bracketing %.2f-%.2f (q=%.2f)\n", tau_start, tau, q )
    }

    tau_avg = NA
    while( abs( tau - tau_start ) > tolerance ) {
        tau_avg = (tau + tau_start) / 2
        #scat( "tau_avg = %.2f\n", tau_avg )
        denom <- vj + tau_avg ^ 2
        q = sum( resid2 / denom )
        if ( q < threshold ) {
            tau = tau_avg
        } else {
            tau_start = tau_avg
        }
    }
    #scat( "Search done: %.2f\n", tau_avg )

    return( list( tau = tau_avg,
                  Q = q ) )
}





#' Calculate confidence interval for cross site variation using Q-statistic test
#' inversion.
#'
#' Two versions of this, one that takes raw data (in which case a linear model
#' is fit to it to get individual site impact estimates), and a second, '_stat',
#' that takes precalculated point estimates and standard errors.
#'
#' The first simply calculates ATE_hat and SE_hat, and then calls the second.
#'
#' Code based on prior work of Catherine Armstrong.  Augmented to do binary
#' search for endpoints for the test inversion procedure to generate confidence
#' intervals. Also calls optim() to find the point estimate by maximizing
#' p-value.
#'
#' @param ATE_hat List of estimated ATEs, one for each site.
#' @param SE_hat List of associated estimated SEs, assumed to be known.
#' @param alpha The level of the test.  The CI will be a 1-2alpha confidence
#'   interval.
#' @param calc_CI Logical
#' @param verbose Print out a few diagnostics as the CI is being calculated, for
#'   those curious.
#' @param tolerance For optimization and search for CI endpoints, how close to
#'   optimal do we want?
#' @param mean_method How to calculate the mean effect.  'weighted' is a
#'   precision weighted mean. 'raw' is the simple mean.  Or pass a number which
#'   will be used as the mean.
#' @return List of several estimated quantities: cross-site variation, p-value
#'   for any variation, whether the p-value is less than passed alpha, the Q
#'   statistic itself, and the confidence interval range.
#'
#' @export
#' @rdname analysis_Qstatistic
#' @examples
#' analysis_Qstatistic_stat( rnorm(30, sd=0.1), rep(0.1,30), calc_CI = TRUE )
analysis_Qstatistic_stat <- function( ATE_hat, SE_hat, alpha = 0.05,
                                      calc_CI = FALSE,
                                      verbose = FALSE,
                                      tolerance = 0.00001,
                                      mean_method = c( "weighted", "raw" ) ) {

    bj = ATE_hat
    vj = SE_hat^2
    wj <- 1 / vj
    s <- length(ATE_hat)

    if ( !is.numeric( mean_method ) ) {
        stopifnot( mean_method %in% c( "weighted", "raw" ) )

        mean_method = match.arg(mean_method)
        if ( mean_method == "weighted" ) {
            bbar <- sum(wj * bj) / (sum(wj))
        } else {
            bbar = mean( bj )
        }
    } else {
        bbar = mean_method
    }




    q <- sum((bj - bbar) ^ 2 / vj)
    pval <- pchisq(q, df = (length(bj) - 1), lower.tail = FALSE)
    reject <- (pval < alpha)
    if ( verbose ) {
        scat( "Stats: bbar = %.3f\n s = %d\n alpha=%0.2f q=%.2f\n",
              bbar, s, alpha, q )
        scat( "   -> Reject: %d\n", as.numeric( reject ) )
        scat( "Summary of weights:\n" )
        print( summary( wj ) )
    }

    stopifnot( length( ATE_hat ) == length( SE_hat ) )

    if (calc_CI) {
        lowQ <- qchisq(alpha, s - 1)
        highQ <- qchisq(1 - alpha, df = (s - 1))
        qmax = q

        resid2 = (bj - bbar) ^ 2

        low_bound = find_CI_endpoint( resid2, vj, threshold=highQ,
                                      tau_start=0, tau=1,
                                      tolerance=tolerance)
        ts = ifelse( is.na( low_bound$tau ), 0, low_bound$tau )
        high_bound = find_CI_endpoint( resid2, vj, threshold=lowQ,
                                       tau_start=ts, tau=pmax( 1, ts * 2 ),
                                       tolerance=tolerance)

        if ( is.na( high_bound$tau ) ) {
            # Too little variation
            CI_low = 0
            CI_high = 0
            tau_hat = 0

        } else {
            if ( verbose ) {
                scat( "Found CI: %.2f (%.2f) - %.2f (%.2f)\n",
                      low_bound$tau, low_bound$Q,
                      high_bound$tau, high_bound$Q )
            }
            CI_low = ts
            CI_high = high_bound$tau

            tau_hat = NA
            if ( !is.na( low_bound$tau ) ) {
                tau_start = (low_bound$tau + high_bound$tau)/2
                rs = optim( tau_start, find_p_value, resid2=resid2, vj=vj,
                            lower = low_bound$tau, upper = high_bound$tau,
                            method="L-BFGS-B",
                            control = list( fnscale = -1 ) )
                tau_hat = rs$par
            } else {
                tau_hat = 0
            }
        }
    } else { # we didn't bother calculating the Conf Int.
        tau_hat = NA
        CI_low = NULL
        CI_high = NULL
    }

    return(list(tau_hat = tau_hat, p.value = pval, reject = reject,
                Q = q, CI_low = CI_low, CI_high = CI_high))
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
analysis_Qstatistic <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                                alpha = 0.05, calc_CI = FALSE) {
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
