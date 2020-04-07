#' Estimate a series of linear models using different weighting schemes and
#' standard errors.
#'
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param siteID If not null, name of siteID that has randomization blocks
#' @param control_formula The control_formula argument must be of the form ~ X1
#'   + X2 + ... + XN. (nothing on left hand side of ~)
#' @param weight_LM_method Argument passed to weight.method of
#'   weighted_linear_estimators
#' @param weight_LM_scale_weights Argument passed to sclae.weights of
#'   weighted_linear_estimators
#' @param block.stats Table of precomputed block-level statistics (optional, for
#'   speed concerns; this gets precomputed in compare_methods).
#' @param data Dataframe of the data to analyse.
#' @importFrom dplyr group_by ungroup mutate
#' @importFrom sandwich vcovHC vcovCL
#' @importFrom stats coef
#' @return Data frame of the various results.
#' @export
linear_model_estimators <- function(Yobs, Z, B, siteID = NULL, data = NULL,
                                    block.stats = NULL,
                                    control_formula = NULL,
                                    weight_LM_method = "survey",
                                    weight_LM_scale_weights = TRUE) {

  if (!is.null(control_formula)) {
    stopifnot( !is.null(data))
    stopifnot( !missing("Yobs"))
  }
  if (missing("Z") && is.null(data)) {
    data <- Yobs
  }
  if (is.null(data)) {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
  } else {
    if (missing( "Z")) {
      stopifnot(all(c("Yobs", "Z", "B") %in% names(data)))
    } else {
      d2 <- data
      if (!is.null( siteID)) {
        d2$siteID <- data[[siteID]]
        stopifnot(!is.null(d2$siteID))
      }
      d2$Yobs <- eval(substitute(Yobs), data)
      d2$Z <- eval(substitute(Z), data)
      d2$B <- eval(substitute(B), data)
      data <- d2
      rm(d2)
      #data = rename_( data,
      #                Yobs = as.character( substitute( Yobs ) ),
      #                Z = as.character( substitute( Z ) ),
      #                B = as.character( substitute( B ) ) )
    }
  }
  dat <- data
  dat$B <- as.factor(dat$B)
  FEmodels <- fixed_effect_estimators(Yobs, Z, B, siteID = siteID, data=dat,
                                      block.stats = block.stats, control_formula = control_formula )
  weightModels <- weighted_linear_estimators(Yobs ~ Z*B, siteID = siteID, data=dat, control_formula = control_formula,
    weight.method = weight_LM_method, scaled.weights = weight_LM_scale_weights)
  interactModels <- interacted_linear_estimators(Yobs, Z, B, siteID = siteID, data = dat, control_formula = control_formula )
  # combine and return results
  dplyr::bind_rows(FEmodels, weightModels, interactModels)
}
