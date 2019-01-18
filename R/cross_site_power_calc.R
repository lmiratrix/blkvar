

# Calculate power for detecting cross site variation
# (With a ratio test)
# Using work and results in Appendix A of Bloom and Spybrook


calc.omega = function( n.bar, tau, p.tx, ICC, R2.X ) {
    tau * p.tx * (1-p.tx) * n.bar / ( (1 - ICC) * (1-R2.X ) )
}



#' @title Calculate power for cross-site variation
#' @description Calculate power for detecting cross site variation with a ratio
#'   test.  This uses work and results in Appendix A of Bloom and Spybrook
#'
#' @param J Number of sites
#' @param n.bar average people per site
#' @param tau Cross site VARIANCE of site-level ATEs.
#' @param p.tx proportion treated, Default: 0.5
#' @param ICC Control-side ICC, Default: 0.5
#' @param K Number of individual covariates used to explain control-side variation, Default: 0
#' @param R2.X Explanatory power of these individual level covariates, Default: 0
#' @param alpha Significance level, Default: 0.05
#' @return Power
#'
#' @examples
#' # Example in the text
#' calc.power( 80, 60, 0.02, p.tx=0.6, ICC=0.2, K = 0, R2.X=0.25 )
#'
#' @rdname calc.power
#' @export
calc.power = function( J, n.bar, tau, p.tx= 0.5, ICC = 0.5, K = 0, R2.X = 0, alpha=0.05 ) {
    omega = calc.omega( n.bar, tau, p.tx, ICC, R2.X )
    df.num = J-1
    df.denom = J*(n.bar-2) - K
    F.crit = qf( 1-alpha, df.num, df.denom )

    1 - pf( F.crit / (omega+1), df.num, df.denom )
}


#' @title Calculate MDESSD for cross-site variation
#' @description Calculate MDESSD (minimum detectable effect size for cross site standard deviation of effects
#' assuming a ratio test.  This uses work and results found in Appendix A of Bloom and Spybrook
#'
#' @inheritParams calc.power
#' @return MDESSD for given scenario.
#' @seealso calc.power
calc.MDESSD = function( J, n.bar, p.tx= 0.5, ICC = 0.5, K = 0, R2.X = 0, alpha=0.05 ) {
    df.num = J-1
    df.denom = J*(n.bar-2) - K
    F.crit = qf( 1-alpha, df.num, df.denom )
    F.crit.20 = qf( 0.20, df.num, df.denom )
    rat = (1-ICC)*(1-R2.X)/ (n.bar * p.tx * (1-p.tx) )
    sqrt( rat * (F.crit / F.crit.20 - 1 ) )
}

if ( FALSE ) {

    debugonce( calc.power )
    calc.power( 80, 60, 0.02, p.tx=0.6, ICC=0.2, K = 0, R2.X=0.25 )

    # Example in the text
    calc.power( 80, 60, 0.02, p.tx=0.6, ICC=0.2, K = 0, R2.X=0.25 )

    calc.power( J=200, n=200, tau=0.05^2, p.tx=0.5, ICC=0.15, R2.X=0.4, alpha=0.05 )
    grid = expand.grid( J = c(5,10,20,50,100,200),
                        n = c(5,50,500) )
    head( grid )
    grid$tau = c( 1.65, 1.07, 0.78,0.57,0.45,0.37,
                  0.45, 0.30, 0.22, 0.16, 0.13, 0.11,
                  0.14, 0.09, 0.07, 0.05, 0.03, 0.03 )^2


    grid$power = pmap_dbl( grid, calc.power, p.tx=0.5, ICC=0.15, R2.X=0.4, alpha=0.05)
    grid

    calc.power( J=5, n=5, tau=1.65^2, p.tx=0.5, ICC=0.15, R2.X=0.4, alpha=0.05 )
    calc.power( J=5, n=5, tau=1.65^2, p.tx=0.5, ICC=0.15, R2.X=0.4, alpha=0.05 )

}
