## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)
library(blkvar)
library(nlme)
library(ggplot2)
set.seed( 102030 )

## -----------------------------------------------------------------------------
df = gen_dat_no_cov( n.bar=40, J=20,
                     gamma.10 = 0,
                     tau.11.star = 0.2^2,
                     ICC = 0.55,
                     variable.n = TRUE )
head( df )

## -----------------------------------------------------------------------------
re.mod <- nlme::lme(Yobs ~ 0 + Z + sid,
                    data = df,
                    random = ~ 0 + Z | sid,
                    weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                    method = "ML",
                    control=nlme::lmeControl(opt="optim",returnObject=TRUE))

# Estimated average treatment impact (with SEs)
ATE = nlme::fixef( re.mod )[[1]]
SE.ATE = sqrt( vcov( re.mod )[1,1] )

# Estimated cross site variation.
vc <- nlme::VarCorr(re.mod)
tau.hat <- as.numeric( vc["Z","StdDev"] )
tau.hat

## -----------------------------------------------------------------------------
re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                         data=df,
                         method = "ML",
                         control=nlme::lmeControl(opt="optim",returnObject=TRUE))

myanova = anova(re.mod.null, re.mod)
myanova
p.value.anova = myanova[2,9]

p.variation = ( p.value.anova / 2 )  # divide by 2 by same logic as chi-squared test.
p.variation

## -----------------------------------------------------------------------------
intervals( re.mod, which="var-cov" )

## -----------------------------------------------------------------------------
s = length( unique( df$sid ) )
ols <- nlme::gls(Yobs ~ 0 + Z*sid - Z,
                 data=df,
                 method = "ML",
                 control=nlme::lmeControl(opt="optim",returnObject=TRUE))

bj <- ols$coefficients[(s+1):(2*s)]
vj <- (coef(summary(ols))[(s+1):(2*s),2])^2
wj <- 1/vj

# Look at our point estimates, just for kicks
summary( bj )

## -----------------------------------------------------------------------------
# ATE
bbar <- sum(wj*bj)/(sum(wj))
bbar

# Dispersion measure
q <- sum((bj - bbar)^2/vj)
q

# P-value
pval <- pchisq(q,df=(s-1),lower.tail=FALSE)
pval

## ---- fig.height=3------------------------------------------------------------
dstn = rchisq( 100000, df= (length(bj)-1))
ggplot2::qplot( dstn, binwidth= 1 ) + geom_vline(xintercept=q )

## -----------------------------------------------------------------------------
tau_test <- seq(0, 2, 0.005) 

alpha = 0.05

lowbound <- qchisq(alpha,s-1)
highbound <- qchisq(1 - alpha,df=(s-1))
lowbound
highbound

## -----------------------------------------------------------------------------
q_invert <- c()
CI_95 <- c()
for (i in 1:length(tau_test)){
    #calculate new denominator that takes into account tau^2
    denom <- vj + tau_test[i]^2
    #new q-stat
    q_invert[i] <- sum((bj - bbar)^2/denom)
    #Compare q.stat to lower and upperbounds (Weiss et al, JREE p. 55)
    CI_95[i] <- (q_invert[i]>=lowbound & q_invert[i]<=highbound)
}

## -----------------------------------------------------------------------------

df = data.frame( tau = tau_test, q = q_invert, hit = CI_95 )

# For illustration
lowlowbound <- qchisq(alpha/2,s-1)
highhighbound <- qchisq(1 - alpha/2,df=(s-1))

ggplot( df, aes( tau, q, col=hit ) ) +
    geom_point() +
    geom_hline(yintercept=c(lowbound,highbound, lowlowbound, highhighbound), 
               col=c("orange","blue", "grey", "grey")  )


## -----------------------------------------------------------------------------
if ( length(tau_test[CI_95 == 1]) == 0 ) {
    CI_low <- NA
    CI_high <- NA
} else {
    CI_high <- max(tau_test[CI_95 == 1])
    CI_low <- min(tau_test[CI_95 == 1])
}

CI_low
CI_high

## -----------------------------------------------------------------------------
pval.low <- pchisq(q_invert,df=(s-1),lower.tail=FALSE)
pval.high <- pchisq(q_invert,df=(s-1),lower.tail=TRUE)
df$pval = pmin( pval.low, pval.high )
ggplot( df, aes( tau, pval ) ) +
    geom_line() +
    geom_hline( yintercept = alpha )

pos = which.max( df$pval )

tau.hat.2 = df$tau[[ pos ]]
tau.hat.2

