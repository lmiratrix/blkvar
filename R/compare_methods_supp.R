#
# Utility function to compare a variety of estimation methods on the same dataset.
#

# Calculate estimates for the Multilevel modeling methods
compare_MLM_methods <- function( Yobs, Z, B, siteID = NULL, data = NULL, control.formula = NULL ) {
  if (!is.null(data)) {
    if (missing( "Yobs")) {
      data <- data.frame(Yobs = data[[1]], Z = data[[2]], B = data[[3]])
    } else {
      d2 <- data
      if (!is.null(siteID)) {
        d2$siteID <- data[[siteID]]
        stopifnot(!is.null(d2$siteID))
      }
      d2$Yobs <- eval(substitute(Yobs), data)
      d2$Z <- eval(substitute(Z), data)
      d2$B <- eval(substitute(B), data)
      data <- d2
      rm(d2)
    }
  } else {
    data <- data.frame(Yobs = Yobs, Z = Z, B = B)
    if (!is.null( siteID)) {
      data$siteID <- siteID
      siteID = "siteID"
      }
  }
  stopifnot(length(unique(data$Z)) == 2)
  stopifnot(is.numeric(data$Yobs))
  RICC <- estimate_ATE_RICC(Yobs, Z, B, data = data, REML = TRUE, control.formula = control.formula )
  FIRC <- estimate_ATE_FIRC(Yobs, Z, B, data = data, siteID = siteID, REML = TRUE, include.testing = FALSE, control.formula = control.formula )
  RIRC <- estimate_ATE_RIRC( Yobs, Z, B, data=data, REML = TRUE, include.testing = FALSE, control.formula = control.formula )
  mlms <- data.frame(method=c("RICC", "FIRC", "RIRC"), tau = c( RICC$ATE, FIRC$ATE, RIRC$ATE ), SE = c( RICC$SE.ATE, FIRC$SE.ATE, RIRC$SE.ATE ),
    stringsAsFactors = FALSE)
  if (!is.null(control.formula)) {
    mlms$method <- paste0(mlms$method, "-adj")
  }
  mlms
}


### This function is not called elswhere - is it needed?

# #' Function to compare estimators using the true values (i.e., whole science
# #' table).  This is for simulaton studies where we know all the potential
# #' outcomes.
# #'
# #' Function that returns some variance function estimates based on all potential
# #' outcomes.
# #'
# #' @param Y vector of all outcomes
# #' @param Z vector that indicates if outcome is under treatment or control
# #' @param B block ids
# #' @param data alternatively is matrix of Y,Z,B
# #' @param p.mat  matrix with first column B,second column prop treated in
# #'   that block, p
# #' @importFrom stats aggregate lm quantile rnorm sd var
# #' @importFrom msm deltamethod
# #'

compare_methods_oracle <- function(Y, Z, B, p.mat, data = NULL) {
  if (!is.null(data)) {
    Y <- data[ ,1]
    Z <- data[ ,2]
    B <- data[, 3]
  }
  n <- length(Y)
  Y1 <- Y[1:(n / 2)]
  Y0 <- Y[(n / 2 + 1):n]
  #Quick test that input is correct
  if (is.numeric(Z) == FALSE) {
    stop("Treatment indicator should be vector of ones and zeros")
  }
  if ((sum(Z == 1) + sum(Z == 0)) != n){
    stop("Treatment indicator should be vector of ones and zeros")
  }
  #Get data into table
  s.tc.bk <- aggregate(list(s.tc.bk = Y1 - Y0), list(B = B[1:(n / 2)]), FUN = s_tc_func)
  data.table <- block_data_sim(Y, Z, B, p.mat)
  K <- max(data.table$B)
  n <- sum(data.table$n1) + sum(data.table$n0)
  data.table$nk <- data.table$n1 + data.table$n0
  #Split into big and small blocks
  data.big <- data.table[data.table$n1 > 1 & data.table$n0 > 1, ]
  data.small<-data.table[data.table$n1 == 1 | data.table$n0 == 1, ]
  #Calculate variance for big blocks
  var_big <- sum(data.big$se_ney ^ 2 * (data.big$nk) ^ 2) / n ^ 2
  #I f no small blocks, just use big block variance
  if (nrow(data.small) == 0) {
    hybrid_m_est <- var_big
    hybrid_p_est <- var_big
    plug_in_big_est <- var_big
  } else{
    # Otherwise calculate variance estimates for each method considered
    combine_table <- merge(s.tc.bk, data.table, by = "B")
    mod.small <- combine_table[data.table$n1 == 1 | data.table$n0 == 1, ]
    var_small <- sum((mod.small$se_ney ^ 2 - mod.small$s.tc.bk / mod.small$nk) * (mod.small$nk) ^ 2) / n ^ 2
    hybrid_m_est <- hybrid_m_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big + var_small
    hybrid_p_est < -hybrid_p_small(data.small) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big + var_small
    plug_in_big_est <- plug_in_big(data.small, data.big) * sum(data.small$nk) ^ 2 / n ^ 2 + var_big
  }
  #Get trt effect estimates and aggregate
  tau_vec < -data.table$Ybar1 - data.table$Ybar0
  tau_est <- sum(tau_vec * data.table$nk) / n

  #Get linear model estimates (sandwich)
  M0 <- lm(Y ~ Z + as.factor(B))
  tau_est_lm <- summary(M0)$coefficients[2, 1]
  var_est_lm <- sandwich::vcovHC(M0, type = "HC1")[2, 2]
  blk.id <- as.numeric(as.factor(B))
  num.c.vec <- data.table[match(blk.id,data.table$B), ]$nk / data.table[match(blk.id, data.table$B), ]$n0
  num.t.vec <- data.table[match(blk.id,data.table$B), ]$nk / data.table[match(blk.id, data.table$B), ]$n1
  weight.vec <- num.c.vec * (1 - Z) + num.t.vec * (Z)
  p.overall <- sum(data.table$n1) / sum(data.table$n1 + data.table$n0)
  weight.vec <- num.c.vec * (1 - Z) * (1 - p.overall) + num.t.vec * (Z) * p.overall

  # Get linear model estimates (Weighted OLS)
  M1 <- lm(Y ~ Z + as.factor(B), weights = weight.vec)
  tau_est_weighted_fixed <- summary(M1)$coefficients[2, 1]
  var_est_weighted_fixed <- summary(M1)$coefficients[2, 2] ^ 2

  # Aggregate into summary table
  methods_list <- c("hybrid_m", "hybrid_p", "plug_in_big", "fixed effects-no int", "weighted regression-fixed effects")
  var_estimates <- c(hybrid_m_est, hybrid_p_est, plug_in_big_est, var_est_lm, var_est_weighted_fixed)
  tau_estimates <- c(rep(tau_est, 3), tau_est_lm, tau_est_weighted_fixed)
  summary_table <- data.frame(Methods = methods_list, tau = tau_estimates, SE = sqrt(var_estimates))
  return(summary_table)
}

# #### Testing and demo of this code ####

# if  (FALSE ) {
    # dat = make_obs_data( n_k = 4:10, p = 0.2 )
    # dat
    # table( dat$Z, dat$B )
    # compare_methods( data = dat[ c("Yobs", "Z","B" ) ] )
# }
