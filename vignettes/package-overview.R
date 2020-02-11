## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library( blkvar )
library( dplyr )
library( ggplot2 )
library( knitr )

## -----------------------------------------------------------------------------
df = gen.dat.no.cov( n.bar=50, J=30,
                         gamma.10 = 0,
                         tau.11.star = 0.2^2,
                         ICC = 0.65,
                         variable.n = TRUE )
head( df )

## -----------------------------------------------------------------------------
compare_methods(Yobs, Z, sid, data=df )

## -----------------------------------------------------------------------------
compare_methods(Yobs, Z, sid, data=df, include.block = FALSE, include.MLM = FALSE,
                include.DB = FALSE)

## -----------------------------------------------------------------------------
df$X1 = df$Y0 + rnorm( nrow( df ) ) # add a fake covariate
df$X2 = df$Y1 + rnorm( nrow( df ) ) # add another fake covariate

rs = compare_methods( Yobs, Z, sid, data=df, control.formula = ~X1 + X2 )
knitr::kable( rs )

## -----------------------------------------------------------------------------
ests <- estimate.ATE.FIRC( Yobs, Z, sid, data=df )
t( ests )

## -----------------------------------------------------------------------------
ests <- estimate.ATE.FIRC( df$Yobs, df$Z, df$sid )

## -----------------------------------------------------------------------------
ests.RIRC = estimate.ATE.RIRC( Yobs, Z, sid, data=df )
t( ests.RIRC )

## -----------------------------------------------------------------------------
t( analysis.Qstatistic( Yobs, Z, sid, data=df, alpha=0.05, calc.CI=TRUE ) )

## -----------------------------------------------------------------------------
compare_methods_variation( Yobs, Z, sid, data=df, long.results=TRUE )

