
##
## Methods of detection and estimation of cross site variation based on the
## Q-statistic/meta analysis approach
##
## The conf int method taken from Catherine's code (but renamed).



# Taken from Catherine's 'weiss.reject' method
# Dataframe has variables of
# - sid (site id)
# - Yobs (outcome)
# - Z (binary treatment 0/1)
estimate.Q.confint <- function(data=data){
    s = length( unique( data$sid ) )

    #calculate Q-statistic
    #run ols model with no intercept and no "treatment intercept"
    ols <- lm(Yobs ~ 0 + Z*factor(sid) - Z, data=data, na.action=na.exclude)
    #model ^^ does not estimate sigma^2 separate for T & C --
    #generated data have same variance across conditions, so temporarily okay
    bj <- ols$coefficients[(s+1):(2*s)]
    vj <- (coef(summary(ols))[(s+1):(2*s),2])^2
    wj <- 1/vj

    bbar <- sum(wj*bj)/(sum(wj))
    q <- sum((bj - bbar)^2/vj)
    pval <- pchisq(q,df=(length(bj)-1),lower.tail=FALSE)

    reject <- (pval < 0.05)

    ## get confidence interval

    ## confidence interval
    #test tau values from 0 to 5 (what is a good way to set this upper bound?)
    tau_test <- seq(0,5,.01) ##come back and make it increment by 0.01

    ## test two-sided version
    lowbound <- qchisq(0.025,s-1)
    low90 <- qchisq(0.05,df=(s-1))
    highbound <- qchisq(0.975,df=(s-1))
    high90 <- qchisq(0.95,df=(s-1))

    ## intialize values
    q_invert <- c()
    CI_95 <- c()
    CI_90 <- c()
    for (i in 1:length(tau_test)){
        #calculate new denominator that takes into account tau^2
        denom <- vj + tau_test[i]^2
        #new q-stat
        q_invert[i] <- sum((bj - bbar)^2/denom)
        #other way: compare q.stat to lower and upperbounds (Weiss et al, JREE p. 55)
        CI_95[i] <- (q_invert[i]>=lowbound & q_invert[i]<=highbound)
        CI_90[i] <- (q_invert[i]>=low90 & q_invert[i]<=high90)
    }

    #extract bounds
    #special condition to explore when tau_test[CI_90 == 1]
    if ( length(tau_test[CI_90 == 1]) == 0 ) {
        CI90_low <- NA
        CI90_high <- NA
        q90_low <- min(q_invert)
        q90_high <- max(q_invert)
    } else {
        CI90_high <- max(tau_test[CI_90 == 1])
        CI90_low <- min(tau_test[CI_90 == 1])
    }

    if ( length(tau_test[CI_95 == 1]) == 0 ) {
        CI_low <- NA
        CI_high <- NA
    } else {
        CI_high <- max(tau_test[CI_95 == 1])
        CI_low <- min(tau_test[CI_95 == 1])
    }


    return(list(reject=reject,
                p.value = pval,
                Q = q,
                CI_low=CI_low,CI_high=CI_high,
                CI90_low=CI90_low,CI90_high=CI90_high))
}




#' Weiss et al. Q-statistic test
#' Dataframe has variables of
#' - sid (site id)
#' - Yobs (outcome)
#' - Z (binary treatment 0/1)
#' @rdname analysis.Qstatistic
#' @export
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


#' Utility wrapper function to just get the p-value without all the other stuff.
#' @export
analysis.Qstatistic = function( df ){

    analysis.Qstatistic.extended(df)$p.value.Q
}







#' Weiss et al. Q-statistic test when there is an individual level covariate X
#'
#' TO DO: Is the above description correct?
analysis.Qstatistic.X = function( df ){
    s = length( unique( df$sid ) ) # max(as.numeric(as.character(df$sid)))
    df$sid = as.factor(df$sid)

    # get the model matrix we need to estimate the relevant betas
    M0 = lm( Yobs ~ sid*(-1+Z) + X, data=df)
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

    ## estimate q and pvalue
    wj = 1/vj
    bbar <- sum(wj*bj)/(sum(wj))
    q <- sum((bj - bbar)^2/vj)

    pval <- pchisq(q, df=(s-1),lower.tail=FALSE)
    return(pval)
}







# Testing
if ( FALSE ) {
    dat = catherine.gen.dat( 0.2, 0.2, 30, 50 )
    head( dat )

    fit.FIRC( dat )

    estimate.Q.confint( dat )


    # sanity check: same Q-statistic and same p-value?
    cTst = estimate.Q.confint( dat )
    cTst

    source( "detection_methods.R")
    qTst = analysis.Qstatistic( dat )
    qTst

    cTst$p.value

    qTst - cTst$p.value
}


## --------- Blending the MLM and the Q-statistic approach ---------



# This is unfinished code trying to do the Q-statistic with the random effects
# from a MLM.
analysis.MLM = function( df ) {
    # Fit MLM
    M0 = lmer( Yobs ~ 1 + Z * X + (1|sid), data=df )
    tstat = fixef( M0 ) / se.coef(M0)$fixef
    2 * pnorm( - abs( tstat[["Z:X"]] ) )

    # Get EB estimates of sites
    beta.j = stuff

    # Get nominal SEs

    # Calculate Q
    q <- sum((bj - bbar)^2/vj)

    pval <- pchisq(q, df=(s-1),lower.tail=FALSE)

    pval

}
