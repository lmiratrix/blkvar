#'
#' This file contains functions implementing the various methods for detecting
#' and estimating treatment effect variation across sites in a multi-site trial
#' using methods of moment (MoM) approach.
#'

#' We have:
#'   * Unweighted MoM for both simulations (with and without covariate)
#'   * Weighted MoM for both simulations (with and without covariate)
#'


<<<<<<< HEAD
## MoM Unweighted 
=======
#' MoM Unweighted
#'
#' @param df  The dataframe to analyze, with sid, Z, X and Yobs
>>>>>>> 9f54f7c7d24902d3d0a6f8090a369ee1405b6432
analysis.MoM.unweighted = function( df ){
  s = max(as.numeric(as.character(df$sid)))
  df$sid = as.factor(df$sid)

  # get the model matrix we need to estimate the relevant betas
  M0 = lm( Yobs ~ sid*(-1+Z)+X, data=df)
  mmat = as.data.frame(model.matrix(M0))
  mmat$`sid1:Z` = with(mmat,sid1*Z)
  mmat$Z = NULL
  mmat$Yobs = df$Yobs

  # Estimate our model
  Mtest = lm(Yobs~.,data=mmat)

  # Tidy this up!
  coef_table = summary(Mtest)$coef
  hetero_vars = grep("Z",rownames(coef_table))

  ### pool treatment effect estimates and compute their variance estimates
  bj = coef_table[hetero_vars,"Estimate"]
  vj = coef_table[hetero_vars,"Std. Error"]^2

  ## Var bj
  var_bj = var(bj)
  # Average variance bj
  mean_vj = mean(vj)

  MoM_unweighted = var_bj - mean_vj

  return(MoM_unweighted)
}


<<<<<<< HEAD

## MoM weighted 
=======
## MoM weighted
>>>>>>> 9f54f7c7d24902d3d0a6f8090a369ee1405b6432
analysis.MoM = function( df ){
  s = max(as.numeric(as.character(df$sid)))
  df$sid = as.factor(df$sid)

  # get the model matrix we need to estimate the relevant betas
  M0 = lm( Yobs ~ sid*(-1+Z)+X, data=df)
  mmat = as.data.frame(model.matrix(M0))
  mmat$`sid1:Z` = with(mmat,sid1*Z)
  mmat$Z = NULL
  mmat$Yobs = df$Yobs

  # Estimate our model
  Mtest = lm(Yobs~.,data=mmat)

  # Tidy this up!
  coef_table = summary(Mtest)$coef
  hetero_vars = grep("Z",rownames(coef_table))

  ### pool treatment effect estimates and compute their variance estimates
  bj = coef_table[hetero_vars,"Estimate"]
  vj = coef_table[hetero_vars,"Std. Error"]^2
<<<<<<< HEAD
  
  #browser()
  
=======

  browser()

>>>>>>> 9f54f7c7d24902d3d0a6f8090a369ee1405b6432
  ## estimate q
  wj = 1/vj
  wjsq = wj^2

  ## weight*site-specific treatment effect
  wjT = wj*bj
  wjTsq = wj*(bj^2)

  ## calculate Q
  sumwj = sum(wj) #total weights
  sumwjsq = sum(wjsq) # total weight squared
  sumwjT = sum(wjT) #sum of weights interecations with site-specific effect
  sumwjTsq = sum(wjTsq) #same but for squared
  sumwjT_sq = sumwjT^2 #and we square again
  Q = sumwjTsq - sumwjT_sq/sumwj

  ## calculate c
  c = sumwj - (sumwjsq)/sumwj

  ## calculate tau^2 hat
  tausq = (Q - (s - 1))/c

<<<<<<< HEAD
  
  return(tausq)
}



## MoM Unweighted, no covarite 
analysis.MoM.unweighted = function( df ){
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
  
  ## Var bj
  var_bj = var(bj)
  # Average variance bj
  mean_vj = mean(vj)
  
  MoM_unweighted = var_bj - mean_vj
  
  return(MoM_unweighted)
}

## MoM weighted, no covariate 
analysis.MoM.weighted = function( df ){
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
  
  #browser() = command to test whether everything within the function works
  
  ## estimate q and pvalue
  wj = 1/vj
  wjsq = wj^2
  
  ## weight*site-specific treatment effect
  wjT = wj*bj
  wjTsq = wj*(bj^2)
  
  ## calculate Q
  sumwj = sum(wj) #total weights
  sumwjsq = sum(wjsq) # total weight squared
  sumwjT = sum(wjT) #sum of weights interecations with site-specific effect
  sumwjTsq = sum(wjTsq) #same but for squared
  sumwjT_sq = sumwjT^2 #and we square again
  Q = sumwjTsq - sumwjT_sq/sumwj
  
  ## calculate c
  c = sumwj - (sumwjsq)/sumwj
  
  ## calculate tau^2 hat
  tausq = (Q - (s - 1))/c
  
  return(tausq)
  
}
=======
  return(tausq)
}
>>>>>>> 9f54f7c7d24902d3d0a6f8090a369ee1405b6432
