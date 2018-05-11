#
# Utility function to compare a variety of estimation methods on the same dataset.
#


# Multilevel modeling methods
compare.MLM.methods = function( Y, Z, B) {
    FIRC = estimate.ATE.FIRC( Y, Z, B, REML = TRUE )
    RIRC = estimate.ATE.RIRC( Y, Z, B, REML = TRUE )
    mlms = data.frame( method=c("FIRC", "RIRC"),
                       tau = c( FIRC$ATE, RIRC$ATE ),
                       SE = c( FIRC$SE.ATE, RIRC$SE.ATE ),
                       stringsAsFactors=FALSE )

    mlms
}

#' Block variance method comparison function
#'
#' This function calculates the point estimates and SE estimates for a variety
#' of the blocked designs.
#'
#' @param Y vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data matrix of Y, Z, B, as alternative to using vectors.  The columns
#'   of the data matrix have to be in the order of Y, Z, B as column 1, 2, 3.
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
compare.methods<-function(Y, Z, B, data=NULL){
    if(!is.null(data)){
        Y<-data[,1]
        Z<-data[,2]
        B<-data[,3]
    }
    n<-length(Y)

    #Quick check that input is correct
    if(is.numeric(Z)==FALSE){
        stop("Treatment indicator should be vector of ones and zeros")
    }
    if((sum(Z==1)+sum(Z==0))!=n){
        stop("Treatment indicator should be vector of ones and zeros")
    }
    #Get data into table
    data.table<-block.data(Y,Z,B)

    methods_list<-c("hybrid_m", "hybrid_p", "plug_in_big")
    fits = sapply( methods_list, function( m ) {
        dd = fitdata.sumtable( data.table=data.table, method=m, throw.warnings=FALSE )
        c( dd$tau_est, dd$se_est )
    } )

    # Aggregate into summary table
    SE_estimates<-fits[2,]
    tau_estimates<-fits[1,]
    summary_table<-data.frame( method = methods_list,
                               tau = tau_estimates,
                               SE = SE_estimates, stringsAsFactors = FALSE )

    lms = linear.model.estimators( Y, Z, B )

    dt = convert.table.names( data.table )

    # Design based methods
    RCT.yes.fi = calc.RCT.Yes.SE( dt, method="finite", weight="individual" )
    RCT.yes.fs = calc.RCT.Yes.SE( dt, method="finite", weight="site" )
    RCT.yes.si = calc.RCT.Yes.SE( dt, method="superpop", weight="individual" )
    RCT.yes.ss = calc.RCT.Yes.SE( dt, method="superpop", weight="site" )
    rctyes = dplyr::bind_rows( RCT.yes.fi, RCT.yes.fs, RCT.yes.si, RCT.yes.ss )
    rctyes$method = with( rctyes, paste( "RCT.yes (", weight, "-", method, ")", sep="" ) )
    rctyes$weight = NULL
    names(rctyes)[1] = "tau"

    mlms = compare.MLM.methods( Y, Z, B )


    summary_table = dplyr::bind_rows( summary_table, lms, rctyes, mlms )
    return(summary_table)
}

if  (FALSE ) {
    dat = make.obs.data( n_k = 4:10, p = 0.2 )
    dat
    table( dat$Z, dat$blk )
    compare.methods( data = dat[ c("Yobs", "Z","blk" ) ] )

    debug( fitdata )
    fitdata( dat$Yobs, dat$Z, dat$blk, method="hybrid_p")
}
