####                     Covariate inclusion attempt                   ####

#' Weiss et al. Q-statistic test when there is an individual level covariate X
#' @param data input
#' TO DO: Is the above description correct?
analysis_Qstatistic_X <- function(data) {
  s <- length(unique(data$B)) # max(as.numeric(as.character(data$B)))
  data$B <- as.factor(data$B)
  
  # get the model matrix we need to estimate the relevant betas
  M0 <- lm(Yobs ~ B * (-1 + Z) + X, data = data)
  mmat <- as.data.frame(model.matrix(M0))
  mmat$`B1:Z` <- with(mmat, B1 * Z)
  mmat$Z <- NULL
  mmat$Yobs <- data$Yobs
  
  # Estimate our model
  Mtest <- lm(Yobs~., data = mmat)

  # Tidy this up!
  coef_table <- summary(Mtest)$coef
  hetero_vars <- grep("Z", rownames(coef_table))
  
  ### pool treatment effect estimates and compute their variance estimates
  bj <- coef_table[hetero_vars, "Estimate"]
  vj <- coef_table[hetero_vars, "Std. Error"] ^ 2

  ## estimate q and pvalue
  wj <- 1 / vj
  bbar <- sum(wj * bj) / (sum(wj))
  q <- sum((bj - bbar) ^ 2 / vj)
  pval <- pchisq(q, df = (s - 1),lower.tail = FALSE)
  return(pval)
}