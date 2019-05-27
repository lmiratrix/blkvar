#
# Utility function to compare a variety of estimation methods on the same dataset.
#


# Calculate estimates for the Multilevel modeling methods
compare.MLM.methods = function( Yobs, Z, B, siteID = NULL, data = NULL ) {

    if( !is.null(data) ){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
        } else {
            if ( !is.null( siteID ) ) {
                siteIDv = data[[siteID]]
                stopifnot( !is.null( siteIDv ) )
            }
            data = data.frame( Yobs = eval( substitute( Yobs ), data ),
                               Z = eval( substitute( Z ), data ),
                               B = eval( substitute( B ), data) )
            if ( is.null( siteID ) ) {
                data$siteID = data$B
            } else {
                data$siteID = siteIDv
                siteID = "siteID"
            }
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
            siteID = "siteID"
        }
    }
    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )

    RICC = estimate.ATE.RICC( Yobs, Z, B, data=data, REML = TRUE )
    FIRC = estimate.ATE.FIRC( Yobs, Z, B, data=data, siteID = siteID, REML = TRUE, include.testing = FALSE )
    RIRC = estimate.ATE.RIRC( Yobs, Z, B, data=data, REML = TRUE, include.testing = FALSE )
    mlms = data.frame( method=c("RICC", "FIRC", "RIRC"),
                       tau = c( RICC$ATE, FIRC$ATE, RIRC$ATE ),
                       SE = c( RICC$SE.ATE, FIRC$SE.ATE, RIRC$SE.ATE ),
                       stringsAsFactors=FALSE )

    mlms
}

#' Get list of methods in package along with characteristics of those methods
#'
#' Each method name comes with whether it targets average impact for individuals
#' or average impact of sites, and also whether it appears to be a finite-sample
#' method (assuming fixed sites, even if individuals are sampled within site) or
#' superpopulation method (assuming sites/blocks are themselves sampled).
#'
#' @return A tibble of characteristics.
#'
#' @export
method.characteristics = function() {
    # Code to make the hard-coded list of characteristics
    if ( FALSE ) {
        dat = make.obs.data( n_k = 4:10, p = 0.2 )
        a = compare_methods( data = dat[ c("Yobs", "Z","B" ) ] )
        a = a[1]
        print( a, row.names = FALSE )
        a$fullname = a$method
        a$method = c("hybrid_m","hybrid_p","plug_in_big", "DB-FP-Persons",
                   "DB-FP-Sites", "DB-SP-Persons", "DB-SP-Sites", "FE",
                   "FE-Het", "FE-CR", "FE-IPTW(n)", "FE-IPTW", "FE-IPTW-Sites",
                   "FE-Int-Sites", "FE-Int-Persons", "RICC", "FIRC", "RIRC" )
        a$finite = c(1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,0,0)
        a$site = c(0,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,1,1)
        print( a, row.names = FALSE )
        a$weight = ifelse( a$site == 1, "site", "person" )
        a$population = ifelse( a$finite, "finite", "superpop" )
        a$site = a$finite = NULL
        print( a, row.names = FALSE )
        dput( a )
        a$biased = c(0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1)
        datapasta::tribble_paste( a )

}
    tibble::tribble(
        ~fullname,            ~method,  ~weight, ~population, ~biased,
        "hybrid_m",       "hybrid_m", "person",    "finite",       0,
        "hybrid_p",       "hybrid_p", "person",    "finite",       0,
        "plug_in_big",    "plug_in_big", "person",    "finite",       0,
        "DB (individual-finite)",  "DB-FP-Persons", "person",    "finite",       0,
        "DB (site-finite)",    "DB-FP-Sites",   "site",    "finite",       0,
        "DB (individual-superpop)",  "DB-SP-Persons", "person",  "superpop",       0,
        "DB (site-superpop)",    "DB-SP-Sites",   "site",  "superpop",       0,
        "FE",             "FE", "person",    "finite",       1,
        "FE (sand)",         "FE-Het", "person",    "finite",       1,
        "FE (cluster)",          "FE-CR", "person",  "superpop",       1,
        "FE (club)",          "FE-Club", "person",  "superpop",       1,
        "IPTW weighted regression (naive)",     "FE-IPTW(n)", "person",    "finite",       0,
        "IPTW weighted regression",        "FE-IPTW", "person",    "finite",       0,
        "IPTW weighted regression (site)",  "FE-IPTW-Sites",   "site",    "finite",       0,
        "FE interact (site)",   "FE-Int-Sites",   "site",    "finite",       0,
        "FE interact (indiv)", "FE-Int-Persons", "person",    "finite",       0,
        "RICC",           "RICC", "person",    "finite",       1,
        "FIRC",           "FIRC",   "site",  "superpop",       1,
        "RIRC",           "RIRC",   "site",  "superpop",       1
    )
}

#' Block variance method comparison function
#'
#' This function calculates the point estimates and SE estimates for a variety
#' of the blocked designs.
#'
#' @param Y vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column name in data)
#' @param B vector of block ids (or column name in data)
#' @param siteID site ids (variable name as string if data frame passed) (if randomization blocks are nested in site).
#' @param data frame holding Y, Z, B and (possibly a column with name specified by siteID).
#' @param include.MLM Include MLM estimators
#' @param include.DB Include Design-Based estimators (taken from RCTYes documentation and prior literature).
#' @param include.LM Include Linear Model-based estimators (including
#'   Huber-White SEs, etc.)
#' @param include.DBBlended Include DB estimators applied to small block
#'   and classic Neyman to large blocks.
#' @param include.block Include the Pashley blocking variants.
#' @param include.method.characteristics Include details of the methods (target estimands and sampling framework assumed).
#'
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
compare_methods<-function(Yobs, Z, B, siteID = NULL, data=NULL, include.block = TRUE, include.MLM = TRUE,
                          include.DB = TRUE, include.LM = TRUE, include.DBBlended = FALSE,
                          include.method.characteristics = FALSE ){

    # This code block takes the parameters of
    # Yobs, Z, B, siteID = NULL, data=NULL, ...
    # and makes a dataframe with canonical Yobs, Z, B, and siteID columns.
    if(!is.null(data)){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
            n.tx.lvls = length( unique( data$Z ) )
            stopifnot( n.tx.lvls == 2 )
            stopifnot( is.numeric( data$Yobs ) )
        } else {
            if ( !is.null( siteID ) ) {
                siteID = data[[siteID]]
            }
            data = data.frame( Yobs = eval( substitute( Yobs ), data ),
                               Z = eval( substitute( Z ), data ),
                               B = eval( substitute( B ), data) )
            data$siteID = siteID
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
        }
    }
    n = nrow( data )

    #Quick check that input is correct
    if(is.numeric(data$Z)==FALSE){
        stop("Treatment indicator should be vector of ones and zeros")
    }
    if((sum(data$Z==1)+sum(data$Z==0))!=n){
        stop("Treatment indicator should be vector of ones and zeros")
    }

    #Get data into table
    data.table<-calc.summary.stats(Yobs, Z, B, data=data, siteID=siteID, add.neyman = TRUE )

    if ( include.block || include.DBBlended ) {
        method_list = c()
        if ( include.block ) {
            methods_list<-c("hybrid_m", "hybrid_p", "plug_in_big")
        }
        if ( include.DBBlended ) {
            methods_list<-c( methods_list, "rct_yes_all", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small")
        }

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

    } else {
        summary_table = data.frame()
    }

    # stash canonical name for site ID, if we have one.
    if ( ! is.null( siteID ) ) {
        siteID = "siteID"
    }

    if ( include.DB ) {

        # Design based methods
        DB.fi = estimate.ATE.design.based( data.table, siteID=siteID, method="finite", weight="individual" )
        DB.fs = estimate.ATE.design.based( data.table, siteID=siteID, method="finite", weight="site" )
        DB.si = estimate.ATE.design.based( data.table, siteID=siteID, method="superpop", weight="individual" )
        DB.ss = estimate.ATE.design.based( data.table, siteID=siteID, method="superpop", weight="site" )
        DB = dplyr::bind_rows( DB.fi, DB.fs, DB.si, DB.ss )
        DB$method = c( "DB-FP-Persons", "DB-FP-Sites", "DB-SP-Persons", "DB-SP-Sites" ) #with( DB, paste( "DB (", weight, "-", method, ")", sep="" ) )
        DB$weight = NULL
        names(DB)[1] = "tau"

        summary_table = dplyr::bind_rows( summary_table, DB )
    }

    if ( include.LM ) {
        lms = linear.model.estimators( Yobs, Z, B, data=data, siteID = siteID, block.stats = data.table )
        summary_table = dplyr::bind_rows( summary_table, lms )
    }

    if ( include.MLM ) {
        mlms = compare.MLM.methods( Yobs, Z, B, data=data, siteID = siteID )
        summary_table = dplyr::bind_rows( summary_table, mlms )
    }

    # Add info on the methods (e.g., what estimand they are targeting)
    if ( include.method.characteristics ) {
        mc = method.characteristics()
        #mcm = mc$method
        #names(mcm) = mc$fullname
        #summary_table$method = mcm[ as.character( summary_table$method ) ]
        summary_table = merge( summary_table, mc, by="method", all.x=TRUE, all.y=FALSE )
    }

    return(summary_table)
}










#' Function to compare estimators using the true values (i.e., whole science
#' table).  This is for simulaton studies where we know all the potential
#' outcomes.
#'
#' Function that returns some variance function estimates based on all potential
#' outcomes.
#'
#' @param Y vector of all outcomes
#' @param Z vector that indicates if outcome is under treatment or control
#' @param B block ids
#' @param data alternatively is matrix of Y,Z,B
#' @param p.mat  matrix with first column B,second column prop treated in
#'   that block, p
#' @importFrom stats aggregate lm quantile rnorm sd var
#'
#' @export
compare_methods_oracle <-function(Y, Z, B, p.mat, data=NULL){
    if(!is.null(data)){
        Y<-data[,1]
        Z<-data[,2]
        B<-data[,3]
    }
    n<-length(Y)
    Y1<-Y[1:(n/2)]
    Y0<-Y[(n/2+1):n]
    #Quick test that input is correct
    if(is.numeric(Z)==FALSE){
        stop("Treatment indicator should be vector of ones and zeros")
    }
    if((sum(Z==1)+sum(Z==0))!=n){
        stop("Treatment indicator should be vector of ones and zeros")
    }
    #Get data into table
    s.tc.bk<-aggregate(list(s.tc.bk=Y1-Y0), list(B=B[1:(n/2)]), FUN=s.tc.func)
    data.table<-block.data.sim(Y,Z,B, p.mat)
    K<-max(data.table$B)
    n<-sum(data.table$n1)+sum(data.table$n0)
    data.table$nk<-data.table$n1+data.table$n0
    #Split into big and small blocks
    data.big<-data.table[data.table$n1>1&data.table$n0>1,]
    data.small<-data.table[data.table$n1==1|data.table$n0==1,]
    #Calculate variance for big blocks
    var_big<-sum(data.big$se_ney^2*(data.big$nk)^2)/n^2
    #If no small blocks, just use big block variance
    if(nrow(data.small)==0){
        hybrid_m_est<-var_big
        hybrid_p_est<-var_big
        plug_in_big_est<-var_big
    }
    #Otherwise calculate variance estimates for each method considered
    else{
        combine_table<-merge(s.tc.bk, data.table, by="B")
        mod.small<-combine_table[data.table$n1==1|data.table$n0==1,]
        var_small<-sum((mod.small$se_ney^2-mod.small$s.tc.bk/mod.small$nk)*(mod.small$nk)^2)/n^2
        hybrid_m_est<-hybrid_m(data.small)*sum(data.small$nk)^2/n^2+var_big+var_small
        hybrid_p_est<-hybrid_p(data.small)*sum(data.small$nk)^2/n^2+var_big+var_small
        plug_in_big_est<-plug_in_big(data.small, data.big)*sum(data.small$nk)^2/n^2+var_big
    }
    #Get trt effect estimates and aggregate
    tau_vec<-data.table$Y1-data.table$Y0
    tau_est<-sum(tau_vec*data.table$nk)/n

    #Get linear model estimates (sandwich)
    M0 = lm(Y ~ Z + as.factor(B))
    tau_est_lm<-summary(M0)$coefficients[2,1]
    var_est_lm<-sandwich::vcovHC(M0, type = "HC1")[2,2]
    blk.id<-as.numeric(as.factor(B))
    num.c.vec<-data.table[match(blk.id,data.table$B),]$nk/data.table[match(blk.id,data.table$B),]$n0
    num.t.vec<-data.table[match(blk.id,data.table$B),]$nk/data.table[match(blk.id,data.table$B),]$n1
    weight.vec<-num.c.vec*(1-Z) + num.t.vec*(Z)
    p.overall<-sum(data.table$n1)/sum(data.table$n1+data.table$n0)
    weight.vec<-num.c.vec*(1-Z)*(1-p.overall) + num.t.vec*(Z)*p.overall

    # Get linear model estimates (Weighted OLS)
    M1 = lm(Y ~ Z + as.factor(B), weights=weight.vec)
    tau_est_weighted_fixed<-summary(M1)$coefficients[2,1]
    var_est_weighted_fixed<-summary(M1)$coefficients[2,2]^2

    # Aggregate into summary table
    methods_list<-c("hybrid_m", "hybrid_p", "plug_in_big", "fixed effects-no int", "weighted regression-fixed effects")
    var_estimates<-c(hybrid_m_est, hybrid_p_est, plug_in_big_est, var_est_lm, var_est_weighted_fixed)
    tau_estimates<-c(rep(tau_est, 3), tau_est_lm, tau_est_weighted_fixed)
    summary_table<-data.frame( Methods = methods_list,
                               tau = tau_estimates,
                               SE = sqrt( var_estimates ) )
    return(summary_table)
}






#' Compare different estimates of cross site variation.
#'
#' Given a dataframe, use the different methods to pull out point estimates and
#' (if desired) pvalues and return them all.
#'
#' @inheritParams compare_methods
#' @param long.results TRUE means each estimator gets a line in a data.frame.  FALSE gives all as columns in a 1-row dataframe.
#' @param siteID if blocks B nested in sites, then pass the site indicator.
#'
#' @export
compare_methods_variation = function( Yobs, Z, B, siteID = NULL, data = NULL, include.testing = TRUE, long.results = FALSE) {
    if(!is.null(data)){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs<-data[,1],
                               Z<-data[,2],
                               B<-data[,3] )
            n.tx.lvls = table( data$Z )
            stopifnot( n.tx.lvls == 2 )
            stopifnot( is.numeric( data$Yobs ) )
        } else {
            if ( !is.null( siteID ) ) {
                siteID = data[[siteID]]
            }
            data = data.frame( Yobs = eval( substitute( Yobs ), data ),
                               Z = eval( substitute( Z ), data ),
                               B = eval( substitute( B ), data) )
            data$siteID = siteID

            # if ( !missing( siteID ) ) {
            #     data = data.frame( Yobs = eval( substitute( Yobs ), data ),
            #                        Z = eval( substitute( Z ), data ),
            #                        B = eval( substitute( B ), data),
            #                        siteID = eval( substitute( siteID ), data ) )
            # } else {
            #     data = data.frame( Yobs = eval( substitute( Yobs ), data ),
            #                        Z = eval( substitute( Z ), data ),
            #                        B = eval( substitute( B ), data) )
            # }
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
        }
    }

    n<-nrow( data )

    #Quick check that input is correct
    if(is.numeric(data$Z)==FALSE){
        stop("Treatment indicator should be vector of ones and zeros")
    }
    if((sum(data$Z==1)+sum(data$Z==0))!=n){
        stop("Treatment indicator should be vector of ones and zeros")
    }

    if ( !is.null( siteID ) ) {
        siteID = "siteID" #quote( siteID )
    }

    # FIRC model (separate variances)
    FIRC = estimate.ATE.FIRC( Yobs, Z, B, siteID=siteID, data=data, include.testing=include.testing )

    # FIRC model (with pooled residual variances)
    FIRC.pool = estimate.ATE.FIRC( Yobs, Z, B, siteID=siteID, data=data, include.testing=include.testing, pool=TRUE )

    # the random-intercept, random-coefficient (RIRC) model
    RIRC = estimate.ATE.RIRC( Yobs, Z, B, data, include.testing=include.testing )

    # the random-intercept, random-coefficient (RIRC) model
    RIRC.pool = estimate.ATE.RIRC( Yobs, Z, B, data, include.testing=include.testing, pool=TRUE)

    # collect results
    res = data.frame( tau.hat.FIRC = FIRC$tau.hat,
                      tau.hat.RIRC = RIRC$tau.hat,
                      tau.hat.FIRC.pool = FIRC.pool$tau.hat,
                      tau.hat.RIRC.pool = RIRC.pool$tau.hat )

    if ( include.testing ) {
        res$pv.FIRC = FIRC$p.variation
        res$pv.RIRC = RIRC$p.variation
        res$pv.FIRC.pool = FIRC.pool$p.variation
        res$pv.RIRC.pool = RIRC.pool$p.variation
        res$pv.Qstat = analysis.Qstatistic( Yobs, Z, B, data = data )$p.value.Q
    }


    if ( long.results ) {
        if ( include.testing ) {
            res = data.frame( method = c("FIRC", "RIRC", "FIRC.pool", "RIRC.pool", "Q" ),
                              tau.hat = c( as.numeric( res[1:4] ), NA ),
                              pv = as.numeric( res[5:9] ) )
        } else {
            res = data.frame( method = c("FIRC", "RIRC", "FIRC.pool", "RIRC.pool" ),
                              tau.hat = as.numeric( c( res[1:4]) ) )
        }

    }
    res
}



#### Testing and demo of this code ####

if  (FALSE ) {
    dat = make.obs.data( n_k = 4:10, p = 0.2 )
    dat
    table( dat$Z, dat$B )
    compare_methods( data = dat[ c("Yobs", "Z","B" ) ] )

    debug( fitdata )
    fitdata( dat$Yobs, dat$Z, dat$B, method="hybrid_p")
}
