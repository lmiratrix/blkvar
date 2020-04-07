#' Interacted linear regression models
#'
#' These linear models have block by treatment interaction terms.  The final ATE
#' estimates are then weighted average of the block (site) specific ATE
#' estimates.
#'
#'#' If siteID passed, it will weight the RA blocks within site and then average
#' these site estimates.
#'
#' SEs come from the overall variance-covariance matrix.
#'
#' @inheritParams linear_model_estimators
#'
#' @return Dataframe of the different versions of this estimator (person and
#'   site weighted)
#' @export
interacted_linear_estimators <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                                         control_formula = NULL) {
  if (!is.null(control_formula)) {
    stopifnot(!is.null(data))
    stopifnot(!missing( "Yobs"))
  }

  # This code block takes the parameters of
  # Yobs, Z, B, siteID = NULL, data=NULL, ...
  # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
  if (!is.null(data)) {
    if (missing( "Yobs")) {
      data <- data.frame( Yobs = data[[1]], Z = data[[2]], B = data[[3]])
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
    if (!is.null( siteID)) {
      data$siteID = siteID
    }
  }

  data$B <- droplevels(as.factor(data$B))
  J <- length(unique(data$B))
  nj <- table(data$B)
  n <- nrow(data)

  formula <- make_FE_int_formula("Yobs", "Z", "B", control_formula, data)
  M0.int <- lm(formula, data = data)
  ids <- grep( "Z:", names(coef(M0.int)))
  stopifnot(length(ids) == J)
  VC <- vcov(M0.int)
  tau_hats <- coef(M0.int)[ids]

  # site weighting
  if (!is.null( siteID)) {
    # aggregate!
    wts <- data %>% group_by(B, siteID) %>%
        dplyr::summarise(n = n()) %>% group_by(siteID) %>% mutate(wts = n / sum(n))

    # some checks to make sure we are matching RA blocks and sites to the
    # right things
    stopifnot( nrow( wts ) == J )
    nms <- gsub( "Z:B", "", names(tau_hats))
    stopifnot(all(nms == wts$B))
    wts <- wts$wts / sum(wts$wts)
  } else {
    wts <- rep(1 / J, J)
  }

  # the block SEs from our linear model
  SE.hat <- diag(VC)[ids]
  tau.site <- weighted.mean(tau_hats, wts)
  # Calculate SE for tau.site
  SE.site <- sqrt(sum(wts ^ 2 * SE.hat))
  wts.indiv <- nj / n
  tau.indiv <- weighted.mean(tau_hats, wts.indiv)
  SE.indiv <- sqrt(sum(wts.indiv ^ 2 * SE.hat))

  # This is the cautious way we don't need since we have 0s in the off diagonal
  # SE.site = sqrt( t(wts) %*% VC %*% wts )
  # faster way---this should work easily.
  # sqrt( t(wts.indiv) %*% VC %*% wts.indiv )

  interactModels <- data.frame(method = c("FE-Int-Sites", "FE-Int-Persons"), tau = c(tau.site, tau.indiv), SE = c(SE.site, SE.indiv), stringsAsFactors = FALSE)
  if (!is.null(control_formula)) {
    interactModels$method <- paste0(interactModels$method, "-adj")
  }
  interactModels
}

