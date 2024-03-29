% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_ATE_RICC.R, R/estimate_ATE_RIRC.R
\name{estimate_ATE_RICC}
\alias{estimate_ATE_RICC}
\alias{estimate_ATE_RIRC}
\title{Estimate the ATE using Random-Intercept, Constant-Coefficient (RICC) Model.}
\usage{
estimate_ATE_RICC(Yobs, Z, B, data = NULL, REML = TRUE, control_formula = NULL)

estimate_ATE_RIRC(
  Yobs,
  Z,
  B,
  data = NULL,
  include_testing = FALSE,
  REML = !include_testing,
  keep_EB_estimates = TRUE,
  pool = FALSE,
  control_formula = NULL
)
}
\arguments{
\item{Yobs}{Name of outcome variable (assumed to exist in data)}

\item{Z}{vector of assignment indicators (1==treated)}

\item{B}{Name of the block indicator.}

\item{data}{Dataframe with all needed variables.}

\item{REML}{Logical, Restricted maximum likelihood or maximum
likelihood estimation.  Default of TRUE if include_testing =
FALSE and FALSE otherwise. Note that FIRC can fail badly in
estimating cross-site variation when REML=FALSE.}

\item{control_formula}{The control_formula argument must be of the
form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)}

\item{include_testing}{Logical Include likelihood ratio test for
cross-site treatment variation.}

\item{keep_EB_estimates}{TRUE means returned object has EB
estimates.  FALSE means do not keep them.}

\item{pool}{TRUE means tx and co have same reBual variance. FALSE
gives seperate estimates for each (recommended, default).}
}
\description{
This model has a single treatment coefficient, and a random intercept for the
site control average. So it is analogous to a fixed effect model, but with a
random effect.
}
\details{
There is no test for cross site variation for this method, since we assume none.
}
