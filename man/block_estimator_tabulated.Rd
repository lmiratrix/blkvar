% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blocked_experiment_method.R
\name{block_estimator_tabulated}
\alias{block_estimator_tabulated}
\title{Block variance estimation function.}
\usage{
block_estimator_tabulated(
  summary_stats,
  method = c("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all", "rct_yes_small",
    "rct_yes_mod_all", "rct_yes_mod_small"),
  throw.warnings = TRUE
)
}
\arguments{
\item{summary_stats}{Summary statistics of all the blocks in a dataset.  In
particular, this is the output of calc_summary_stats method.}

\item{method}{The method to be used for variance estimation, defauly
"hybrid_m"}

\item{throw.warnings}{TRUE means throw warnings if the hybrid estimators are
breaking down due to violation of assumptions.}
}
\description{
This function takes the block-level summary statistics of a dataset and
returns a treatment effect estimate, variance estimate and some summary info
about the block structure.
}
