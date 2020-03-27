
# # 
# # # Testing
# if ( FALSE ) {
    # dat = catherine_gen_dat( 0.2, 0.2, 30, 50 )
    # head( dat )

    # fit.FIRC( dat )

    # estimate.Q.confint( dat )


    # # sanity check: same Q-statistic and same p-value?
    # cTst = estimate.Q.confint( dat )
    # cTst

    # source( "detection_methods.R")
    # qTst = analysis_Qstatistic( dat )
    # qTst

    # cTst$p.value

    # qTst - cTst$p.value
# }


### Not run currently

## --------- Blending the MLM and the Q-statistic approach ---------


# 
# # This is unfinished code trying to do the Q-statistic with the random effects
# # from a MLM.
# analysis_Q_MLM = function( data ) {
    # # Fit MLM
    # M0 = lme4::lmer( Yobs ~ 1 + Z * X + (1|B), data=data )
    # tstat = lme4::fixef( M0 ) / arm::se.coef(M0)$fixef
    # 2 * pnorm( - abs( tstat[["Z:X"]] ) )

    # # Get EB estimates of sites
    # beta.j = stuff

    # # Get nominal SEs

    # # Calculate Q
    # q <- sum((bj - bbar)^2/vj)

    # pval <- pchisq(q, df=(s-1),lower.tail=FALSE)

    # pval


# }
