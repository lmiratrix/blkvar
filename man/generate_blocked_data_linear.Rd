% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_blocked_data.R
\name{generate_blocked_data_linear}
\alias{generate_blocked_data_linear}
\alias{generate_blocked_data_obs_linear}
\title{Make data from a linear model specification}
\usage{
generate_blocked_data_linear(
  X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65),
  a = 0,
  b = 0,
  ATE = 0.2,
  d = 0
)

generate_blocked_data_obs_linear(
  X = c(0, 2, 3, 19, 20, 21, 24, 31, 32, 40, 41, 43, 45, 55, 60, 65),
  p = 0.5,
  a = 0,
  b = 0,
  ATE = 0.2,
  d = 0,
  method = c("small", "pair", "big", "none")
)
}
\arguments{
\item{X}{vector of indiviual level covariates}

\item{a}{Intercept of Y0}

\item{b}{Main effect of X}

\item{ATE}{Average Tx}

\item{d}{Interaction effect term.}

\item{p}{Proportion of units treated (as close as possible given block sizes)}

\item{method}{How to block}
}
\value{
Dataframe with original potential outcomes and observed outcome based on random assigment.
}
\description{
Function that, given a covariate vector X, returns a dataframe of potential
outcomes and blocks.

Generate data, form_blocks_from_continuous, and randomize within block and generate
observed potential outcomes
}
\details{
This uses the model: Y = a + bX + ATE Z + d X ATE + epsilon
(It will standardize X for this model, and then standardize Y0, Y1 so sd(Y0) = 1.)
}
