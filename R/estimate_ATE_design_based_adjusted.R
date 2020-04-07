#' Adjusted design based estimator for ATE
#'
#' Given two covariates, X1 and X2, calculate adjusted ATE estimates for the
#' multisite trial.  This uses the formula from the Schochet RCT Yes technical
#' document.
#'
#'@param formula  The formula argument must be of the form outcome ~ treatment:block_id.
#'@param control_formula The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)
#'@param data Dataframe with all needed variables.
#' @param siteID Vector of site IDs if there are randomization blocks nested in
#'    site that should be aggregated (will change results for site weighting only).
#'@param method What form of standard errors (finite, superpop, or adjusted superpop)
#'@param weight Site or individual weighting.
#'
#'@return Estimates of ATE along with SEs.
#'@importFrom rlang .data
#'@export
estimate_ATE_design_based_adjusted <- function(formula, control_formula, data, siteID = NULL, method = c("finite", "superpop", "superpop.adj"), weight = c("individual", "site")) {
  stopifnot(!is.null(control_formula))
  # Determine which of the 4 versions of estimator we are doing.
  method <- match.arg(method)
  weight <- match.arg(weight)
  data <- make_canonical_data(formula, control_formula, siteID, data)
  # Get control variables
  c.names <- formula.tools::rhs.vars(control_formula)
  # Count observations, etc.
  N <- nrow(data)  # number of observations
  h <- length( unique( data$B ) ) # number of sites/clusters
  v <- length(c.names) # number of covariates

  # Center the X variables around their overall grand means
  center <- function(x) {
    x - mean(x)
  }
  data <- data %>% dplyr::group_by(B) %>% dplyr::mutate_at(c.names, center) %>% dplyr::ungroup()
  new.form <- make_FE_int_formula( control_formula = control_formula, data = data)
  # Fit a linear model with indicators for each site, and each site by treatment
  # interaction.  Also have covariates entered in not interacted with treatment.
  M0 <- lm( new.form, data = data )
  data$resid <- resid( M0 )

  # Aggregate data by group, calculating various summary statistics.
  # Note: Not all of these stats are used by all the estimators.
  # Key: The MSE.T and MSE.C are from the Schochet formulas.
  sum_tab <- data %>% group_by( B, siteID ) %>% dplyr::summarise(n = n(), nT = sum(Z), nC = n - nT, Ybar.C = mean(Yobs[Z==0]), Ybar.T = mean(Yobs[Z==1]),
    MSE.T = sum(resid[Z==1] ^ 2) / ((N - v) * (nT / N) - 1), MSE.C = sum(resid[Z == 0] ^ 2) / ((N - v) * (nC / N) - 1))
  #                      X1.bar = mean( X ),
  #                      X2.bar = mean( X2 ) )
  # make sure we have one row of stats per randomization block.
  stopifnot(length(unique(data$B)) == nrow(sum_tab))

  # Copy over our block-level impact estimates
  sum_tab$tau_hat_b <- coef(M0)[(h + v + 1):(2 * h + v)]

  # Get our block weights depending on target estimand
  if (weight == "individual") {
    sum_tab$.weight <- sum_tab$n
  } else {
    if (!is.null(siteID)) {
      sum_tab <- sum_tab %>% dplyr::group_by( !!as.name( siteID ) ) %>% dplyr::mutate( .weight = n / sum( n ) ) %>% ungroup()
    } else {
      sum_tab$.weight <- rep(1, h)
    }
  }

  tau_hat <- NA
  SE <- NA

  # calculate overall ATE estimate by taking a weighted average of the
  # estimated block effects.
  tau_hat <- sum( sum_tab$tau_hat_b * sum_tab$.weight ) / sum( sum_tab$.weight )

  if (method == "finite") {
    # Finite population SE calculation
    # finite pop (Neyman)
    w.tot <- sum(sum_tab$.weight)

    # calculate SEs for each block by itself
    sum_tab = dplyr::mutate(sum_tab, block.var = MSE.T / nT + MSE.C / nC )

    # and then take a weighted sum of these
    var <- with(sum_tab, sum(.weight ^ 2 * block.var) / w.tot ^ 2)
    SE <- sqrt(var)
  } else {
    # Superpopulation SE calculation
    if (method == "superpop.adj") {
      # Center our averaged covariates around their site-weighted means.
      X.w <- weighted.mean(sum_tab$X1.bar, w = sum_tab$.weight)
      X2.w <- weighted.mean(sum_tab$X2.bar, w = sum_tab$.weight)
      sum_tab = mutate(sum_tab, X1c = X1.bar - X.w, X2c = X2.bar - X2.w)

      # Our second stage regression
      M.sp <- lm( tau_hat_b ~ 1 + X1c + X2c, weights = .weight, data = sum_tab)
      # The intercept of this model should be our tau_hat
      stopifnot(abs(coef(M.sp)[[1]] - tau_hat) < 0.00001)

      # Now calculate the Superpopoulation SE
      sum_tab$tau.pred <- predict(M.sp)

      # Equation 6.28, page 86 (using the residuals)
      w.bar <- mean(sum_tab$.weight)
      SE.SP <- sum(sum_tab$.weight ^ 2 * (sum_tab$tau_hat_b - sum_tab$tau.pred) ^ 2) / ((h - v - 1) * h * w.bar ^ 2)
      SE <- sqrt(SE.SP)
    } else {
      # First aggregate to get sites, if needed
      if (!is.null(siteID)) {
        sum_tab <- sum_tab %>% group_by(!!as.name(siteID)) %>%
            dplyr::summarise( tau_hat_b = sum(.data$tau_hat_b * .data$.weight) / sum(.data$.weight),
                              .weight = sum(.data$.weight))
        h <- nrow(sum_tab)
      }
      w.bar <- mean(sum_tab$.weight)
      SE.SP <- sum(sum_tab$.weight ^ 2 * (sum_tab$tau_hat_b - tau_hat) ^ 2) / ((h - 1) * h * w.bar ^ 2)
      SE <- sqrt(SE.SP)
    }
  }
  data.frame(tau_hat = tau_hat, SE = SE, weight = weight, method = method, stringsAsFactors = FALSE)
} # end estimator function
