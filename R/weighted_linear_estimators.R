


# Utility function
grab_SE <- function(MOD, coef="Z") {
    vc <- vcov(MOD)
    sqrt(vc[coef, coef])
}




#' Survey-weighted adjusted linear regression
#'
#' Use survey weight regression to reweight blocks to target unbiased ATE estimators.
#'
#' @inheritParams linear_model_estimators
#' @param formula Specification of outcome, treatment, and block ID variables as a formula of form "Yobs ~ Z*B" (order of variables is important).
#' @param scaled.weights Logical Scale the weights by overall proportions of treated and control.
#' @param weight.method Use survey package (svgglm) or classic OLS (lm) for model fitting.
#' @return Dataframe of results for different estimators.
#' @importFrom survey svydesign svyglm
#' @importFrom stats gaussian
#' @export

weighted_linear_estimators <- function(formula, control_formula = NULL, siteID = NULL, data, scaled.weights = TRUE,
                                       weight.method = c("survey", "precision")) {
  weight.method <- match.arg(weight.method)
  data <- make_canonical_data(formula, control_formula, siteID, data)
  data$B <- as.factor(data$B)
  n <- nrow(data)
  J <- length(unique(data$B))
  n.bar <- n / J
  data <- data %>% dplyr::group_by(B) %>%
        dplyr::mutate( p = mean(Z),
                       nj = n(),
                       weight = ifelse( Z, 1/p, 1/(1-p) ),
                       weight.site = weight * n.bar / nj ) %>%
        dplyr::ungroup()

  if (!is.null(siteID)) {
    n.site <- length(unique(data$siteID))
    n.bar <- n / n.site
    # adjust weights to account for blocks nested within sites
    data <- data %>% group_by( siteID ) %>%
            mutate( weight.site = weight * n.bar / n() )
  }
  if (scaled.weights) {
    Z.bar <- mean(data$Z)
    data <- mutate( data,
                      weight = weight * ifelse( Z, Z.bar, (1-Z.bar) ),
                      weight.site = weight.site * ifelse( Z, Z.bar, (1-Z.bar) ) )
  }
  formula <- make_FE_formula("Yobs", "Z", "B", control_formula, data )
  if (weight.method == "survey") {
  	M0w2 <- svyglm( formula,
                       design=svydesign(ids = ~1, weights = ~weight, data=data ),
                       family = gaussian() )

    M0w.site <- survey::svyglm(formula,
                           design = svydesign(ids = ~1, weights = ~weight.site, data = data),
                           family = gaussian() )
  } else {
    M0w2 <- lm( formula, data=data, weights = data$weight)
    M0w.site <- lm( formula, data=data, weights = data$weight.site)
  }
  tau <- coef( M0w2 )[["Z"]]
  SE.w2 <- grab_SE( M0w2 )
  tau.w.site <- coef( M0w.site )[["Z"]]
  SE.w.site <- grab_SE( M0w.site )
  weightModels <- data.frame( method=c("FE-IPTW", "FE-IPTW-Sites"),
                               tau = c( coef( M0w2 )[["Z"]], coef( M0w.site )[["Z"]] ),
                               SE = c( SE.w2, SE.w.site ),
                               stringsAsFactors = FALSE )
  if (!scaled.weights) {
    weightModels$method <- paste0(weightModels$method, "(n)")
  }
  if (weight.method != "survey") {
    weightModels$method <- paste0( weightModels$method, "(pr)" )
  }
  if (!is.null( control_formula)) {
    weightModels$method <- paste0( weightModels$method, "-adj" )
  }
  weightModels
}
