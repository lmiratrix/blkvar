## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library( blkvar )
library( dplyr )
library( knitr )

## -----------------------------------------------------------------------------
df <- gen_dat_no_cov(n.bar = 50, J = 30, gamma.10 = 0, tau.11.star = 0.2 ^ 2, ICC = 0.65, variable.n = TRUE)
head( df)

