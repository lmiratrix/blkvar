analysis.FIRC <- function( Yobs, Z, sid, data=df, REML = FALSE ) {
  
  # get variables
  if ( is.null( df ) ) {
    data = data.frame( Yobs = Yobs, Z = Z, sid= factor(sid) )
  } else {
    sid.name = as.character( quote( sid ) )
    data[ sid.name ] = factor( eval( substitute( sid ), df ) )
  }
  
  #fit multilevel model and extract pvalue
  method = ifelse( REML, "REML", "ML" )
  
  re.mod <- nlme::lme(Yobs ~ 0 + Z + sid,
                      data = df,
                      random = ~ 0 + Z | sid,
                      weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                      method = method,
                      control=nlme::lmeControl(opt="optim",returnObject=TRUE))
  
  # Test for cross site variation
  re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                           data=df,
                           weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                           method = method,
                           control=nlme::lmeControl(opt="optim", returnObject=TRUE))
  
  td = abs(as.numeric( deviance( re.mod.null ) - deviance( re.mod )))
  0.5 * pchisq(td, 1, lower.tail = FALSE )
}

# FIRC function with correct pval
analysis.FIRC <- function( df ) {
  
  #fit multilevel model and extract pvalue
  re.mod <- nlme::lme(Yobs ~ 0 + Z + sid,
                      data = df,
                      random = ~ 0 + Z | sid,
                      weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                      method="ML",
                      control=nlme::lmeControl(opt="optim",returnObject=TRUE))
  
  # Test for cross site variation
  re.mod.null <- nlme::gls(Yobs ~ 0 + Z + sid,
                           data=df,
                           weights = nlme::varIdent(form = ~ 1 | Z), na.action=na.exclude,
                           method = "ML",
                           control=nlme::lmeControl(opt="optim", returnObject=TRUE))
  
  

  # prtAnova = anova(re.mod, re.mod.null)
  # prtAnova$'p-value'
  
  # myanova = anova(re.mod, re.mod.null)
  # myanova[["p-value"]]
 
  myanova = anova(re.mod.null, re.mod)
  myanova[2,9]
  
  # deviance.rand = (-2)*logLik(re.mod) # keeping here so we have it if needed
  # deviance.null = (-2)*logLik(re.mod.null)
  # td = abs(as.numeric( deviance.rand - deviance.null))
  # 0.5 * pchisq(td, 1, lower.tail = FALSE )
  # 
  # pval = summary(myanova)[2,9]
  # pval
}

# Testing
if ( FALSE ) {
  
  df = gen.dat.no.cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)
  
  # head( df )
  # describe.data( df )

  analysis.FIRC( df )
  anova(re.mod.null, re.mod)
  
}

