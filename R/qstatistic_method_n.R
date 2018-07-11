##
## Methods of detection and estimation of cross site variation based on the
## Q-statistic/meta analysis approach
##
##


# This function doesn't incorporate covariates and is used for small sample simulations
analysis.Qstatistic.extended = function( df ){
  s = length( unique( df$sid ) ) # max(as.numeric(as.character(df$sid)))
  df$sid = as.factor(df$sid)
  
  # get the model matrix we need to estimate the relevant betas
  Mtest = lm( Yobs ~ sid*Z - 1 - Z, data=df)
  
  # Tidy this up!
  coef_table = summary(Mtest)$coef
  hetero_vars = grep("Z",rownames(coef_table))
  
  ### pool treatment effect estimates and compute their variance estimates
  bj = coef_table[hetero_vars,"Estimate"]
  vj = coef_table[hetero_vars,"Std. Error"]^2
  
  ## estimate q and pvalue
  wj = 1/vj
  bbar <- sum(wj*bj)/(sum(wj))
  q <- sum((bj - bbar)^2/vj)
  
  pval <- pchisq(q, df=(s-1),lower.tail=FALSE)
  
  data.frame( p.value.Q = pval,
              Q = q, b.bar = bbar, mean.vj = mean( vj ),
              tau.hat.raw = sd( bj ) )
}

# Utility wrapper function to just get the p-value without all the other stuff.
analysis.Qstatistic = function( df ){
  analysis.Qstatistic.extended(df)$p.value.Q
}
