#
# Using linear regression to deal with multi-site or blocked experiments
#

# source( "control_formula_utilities.R")

# Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


# See https://www.jepusto.com/handmade-clubsandwich/
clubsandwich_variance <- function(w, tau_hat_b, tau_hat) {
  W <- sum(w)
  V <- (1 / W ^ 2) * sum((w ^ 2 * (tau_hat_b - tau_hat) ^ 2) / (1 - w / W))
  df.inv <- sum(w ^ 2 / (W - w) ^ 2) - (2 / W) * sum(w ^ 3 / (W - w) ^ 2) + (1 / W ^ 2) * sum(w ^ 2 / (W - w)) ^ 2
  list(var.hat = V, df = 1 / df.inv)
}


#' Fixed effect linear regression
#'
#' Fit regressions with block fixed effects (and a single treatment parameter).
#' Calculate standard errors in a variety of ways.
#'
#' For the club sandwich estimation ("Club"), see code proposed at
#' https://www.jepusto.com/handmade-clubsandwich/
#'
#' @inheritParams linear_model_estimators
#' @return Dataframe of results for different estimators.
#' @export

fixed_effect_estimators <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                                    block.stats = NULL, control.formula = NULL) {
  if (!is.null(control.formula)) {
    stopifnot(!is.null(data))
    stopifnot(!missing("Yobs"))
  }

  # This code block takes the parameters of
  # Yobs, Z, B, siteID = NULL, data=NULL, ...
  # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
  if (!is.null(data)){
    if (missing( "Yobs")) {
      data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]])
      n.tx.lvls <- length(unique(data$Z))
      stopifnot(n.tx.lvls == 2)
      stopifnot( is.numeric(data$Yobs))
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
      data$siteID = siteID
    }
  }

  # make sites RA blocks if there are no sites.
  if (is.null(siteID)) {
    data$siteID = data$B
  }

  # Make control variable function
  formula <- make_FE_formula( "Yobs", "Z", "B", control.formula, data)

  # simple linear model
  M0 <- lm( formula, data = data)
  SE.lm <- summary(M0)$coeff["Z", 2]

  # est ATE (precision weighted)
  tau_hat <- coef(M0)[["Z"]]

  # Huber-White SEs
  vcov_sand <- sandwich::vcovHC(M0, type = "HC1")
  SE.lm.sand <-sqrt( vcov_sand[1,1] )

  # Cluster robust SEs (clustering at site level)
  vcov_clust <- sandwich::vcovCL(M0, data$siteID)
  SE.lm.clust <- sqrt(vcov_clust[1, 1])

  # Cluster robust SEs (clustering at site level using clubSandwich)
  # aggregate!
  if (is.null( block.stats)) {
    block.stats <- calc_summary_stats(Yobs, Z, B, data = data, siteID = siteID, add.neyman = FALSE)
  }
  block.stats <- mutate(block.stats, tau_hat = Ybar1 - Ybar0, prec = n * (n0 / n) * (n1 / n))
  if (!is.null(siteID)) {
    # aggregate blocks into sites and calculate site weights and tau_hats
    block.stats <- block.stats %>% group_by(siteID) %>%
        dplyr::summarise(tau_hat = sum( prec * tau_hat ) / sum(prec), prec = sum(prec))
  }
  cs.var <- clubsandwich_variance(block.stats$prec, block.stats$tau_hat, tau_hat)
  SE.lm.clust.club <- sqrt(cs.var$var.hat)
  FEmodels <- data.frame(method = c("FE", "FE-Het", "FE-CR", "FE-Club"),
                         tau = rep(tau_hat, 4),
                         SE = c(SE.lm, SE.lm.sand, SE.lm.clust, SE.lm.clust.club),
                         stringsAsFactors = FALSE)
  if (!is.null( control.formula)) {
    FEmodels$method <- paste0( FEmodels$method, "-adj")
  }
  FEmodels
}

# #Commentary: some wrong ways of doing things.  The following code contains the
# #wrong way of using weights (i.e. passing directly to lm).

# if ( FALSE ) {

    # M0w = lm( Yobs ~ Z + B, weights=dat$weight, data=dat )
    # SE.w = summary( M0w )$coeff["Z",2]

    # M0w2 = lm( Yobs ~ Z + B, weights=dat$weight, data=dat )
    # SE.w2 = summary( M0w2 )$coeff["Z",2]

    # # Site weighted regression models
    # M0w.site = lm( Yobs ~ Z + B, weights=dat$weight.site, data=dat )
    # tau.w.site = coef( M0w.site )[[2]]
    # SE.w.site = summary( M0w.site )$coeff["Z",2]

# }


# # For the debugging code to get the DGP files
# localsource = function( filename ) {
    # source( file.path( dirname( rstudioapi::getActiveDocumentContext()$path ), filename ) )
# }


# #### Some overall testing code #####
# if ( FALSE ) {
    # dat = make_obs_data(p = 0.2)
    # head( dat )
    # localsource("control_formula_utilities.R" )

    # fixed_effect_estimators( Yobs, Z, B, data=dat )

    # #debug( weighted_linear_estimators.naive )
    # weighted_linear_estimators.naive( Yobs, Z, blk, data=dat )

    # dat2 = rename( dat, B= blk )
    # head( dat2 )
    # weighted_linear_estimators.naive( dat2 )

    # weighted_linear_estimators( Yobs, Z, blk, data=dat )

    # debug( linear_model_estimators)
    # linear_model_estimators( dat$Yobs, dat$Z, dat$blk )
# }



# #### Some overall testing code comparing adustment vs not #####
# if ( FALSE ) {
    # dat = make_obs_data(p = 0.2)
    # head( dat )
    # localsource("control_formula_utilities.R" )

    # set.seed( 1019 )
    # dat = gen_dat( n.bar = 30, J = 30 )
    # nrow( dat )
    # head( dat )
    # dat$X1 = dat$W + rnorm( nrow(dat) )
    # dat$X2 = dat$Y0 + rnorm( nrow( dat ) )

    # wt = weighted_linear_estimators( Yobs ~ Z*sid, data=dat )
    # wt
    # weighted_linear_estimators( Yobs ~ Z*sid, data=dat )
    # weighted_linear_estimators( Yobs ~ Z*sid, data=dat, scaled.weights = FALSE )
    # weighted_linear_estimators( Yobs ~ Z*sid, data=dat, weight.method = "precision" )
    # weighted_linear_estimators( Yobs ~ Z*sid, data=dat, scaled.weights = FALSE,
                                # weight.method = "precision" )

    # rs = linear_model_estimators( Yobs, Z, sid, data=dat )
    # rs.adj = linear_model_estimators( Yobs, Z, sid, data=dat,
                                      # control.formula = ~X1 + X2)
    # rs
    # rs.adj


    # # problem child
    # rs = interacted_linear_estimators( Yobs, Z, sid, data=dat )
    # debug( interacted_linear_estimators )
    # rs.adj = interacted_linear_estimators( Yobs, Z, sid, data=dat,
                                      # control.formula = ~X1 + X2)
    # rs.adj$tau.adj = rs$tau
    # rs.adj = mutate( rs.adj, delta = tau - tau.adj )
    # rs.adj


# }
