% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_multilevel_data.R
\name{make_covariate_set}
\alias{make_covariate_set}
\title{Generate block of noisy covariates}
\usage{
make_covariate_set(X1, k, var_prefix = "X")
}
\arguments{
\item{k}{Number of covariates}

\item{Single}{covariate.  This will be expanded to desired number of covariates.}
}
\value{
tibble with covariates.
}
\description{
Generate block of noisy covariates
}
