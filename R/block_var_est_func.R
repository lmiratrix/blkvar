


#' Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#' Block table function
#'
#' Function that returns a summary of block level estimates.
#' See also `calc.summary.stats` (which does the same thing with slightly different interface).
#'
#' @param Y vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data alternatively is matrix of Y,Z,B
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
block.data<-function(Y, Z, B, data=NULL){
  if(!is.null(data)){
    Y<-data[,1]
    Z<-data[,2]
    B<-data[,3]
  }
  n<-length(Y)
  #Quick test that input is correct
  if(is.numeric(Z)==F){
    return("Treatment indicator should be vector of ones and zeros")
  }
  if((sum(Z==1)+sum(Z==0))!=n){
    return("Treatment indicator should be vector of ones and zeros")
  }

  #First convert block ids into numbers
  blk_id<-factor(B)
  blk_id<-as.numeric(blk_id)

  #Get number of units assigned to each treatment
  #In each block
  n_tk_matrix<-aggregate(list(Num_Trt=Z), list(Blk_ID=blk_id), FUN=sum)
  n_ck_matrix<-aggregate(list(Num_Ctrl=(1-Z)), list(Blk_ID=blk_id), FUN=sum)
  n_ctk_matrix<-merge(n_tk_matrix, n_ck_matrix, by="Blk_ID")
  treated_mat<-cbind(Y[Z==1], blk_id[Z==1])
  control_mat<-cbind(Y[Z==0], blk_id[Z==0])
  Y1_matrix<-aggregate(list(Y1=treated_mat[,1]), list(Blk_ID=treated_mat[,2]), FUN=mean)
  Y0_matrix<-aggregate(list(Y0=control_mat[,1]), list(Blk_ID=control_mat[,2]), FUN=mean)
  Ybar_matrix<-merge(Y1_matrix, Y0_matrix, by="Blk_ID")
  var1_matrix<-aggregate(list(var1=treated_mat[,1]), list(Blk_ID=treated_mat[,2]), FUN=var)
  var0_matrix<-aggregate(list(var0=control_mat[,1]), list(Blk_ID=control_mat[,2]), FUN=var)
  var_matrix<-merge(var1_matrix, var0_matrix, by="Blk_ID")
  overall_mat<-merge(n_ctk_matrix, Ybar_matrix, by="Blk_ID")
  overall_mat<-merge(overall_mat, var_matrix, by="Blk_ID")
  overall_mat$se_ney<-sqrt(overall_mat$var1/overall_mat$Num_Trt + overall_mat$var0/overall_mat$Num_Ctrl)

  return(overall_mat)
}




#' Summarise data by block.
#'
#' @param dat Dataframe with defined Yobs, Z, and blk variables.
#'
#' @return dataframe with summary statistics by block
#' @export
calc.summary.stats = function( dat ) {
    require( tidyverse )
    sdat <- dat %>% dplyr::group_by( blk, Z ) %>%
        dplyr::summarise( n = n(),
                   Y.bar = mean( Yobs ),
                   var = var( Yobs ) )
    sdat <- reshape( as.data.frame(sdat), direction="wide",
                     v.names = c("n","Y.bar", "var"),
                     idvar = "blk",
                     timevar = "Z" )
    sdat$n = sdat$n.0 + sdat$n.1

    sdat
}

convert.table.names = function( data.table ) {
    names( data.table ) = c( "blk", "n.1", "n.0", "Y.bar.1", "Y.bar.0", "var.1", "var.0", "se_ney" )
    data.table$n = with( data.table, n.0 + n.1 )
    data.table
}

if ( FALSE ) {
    dat = make.data( c( 3,4,5 ) )
    dat = add.obs.data(dat)

    calc.summary.stats(dat)
    block.data( dat$Yobs, dat$Z, dat$blk )
}


#' Plot diagnostic function
#'
#' Function that plots variance estimates versus size of treatment group.
#' @param Y vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data matrix of Y, Z, B, as alternative to using vectors
#' @importFrom graphics axis par plot
#' @export
diag.plot<-function(Y, Z, B, data=NULL){
  if(!is.null(data)){
    Y<-data[,1]
    Z<-data[,2]
    B<-data[,3]
  }
  n<-length(Y)
  #Quick test that input is correct
  if(is.numeric(Z)==F){
    return("Treatment indicator should be vector of ones and zeros")
  }
  if((sum(Z==1)+sum(Z==0))!=n){
    return("Treatment indicator should be vector of ones and zeros")
  }
  data.table<-block.data(Y,Z,B)
  #data.table$nk<-data.table$Num_Trt+data.table$Num_Ctrl
  var1_plot<-plot(data.table$Num_Trt[!is.na(data.table$var1)], data.table$var1[!is.na(data.table$var1)], xlab="Number treated", ylab="Estimated variance for treated units", xaxt = "n")
  axis(1, at = 1:max(data.table$Num_Trt[!is.na(data.table$var1)]))
  var0_plot<-plot(data.table$Num_Ctrl[!is.na(data.table$var0)], data.table$var0[!is.na(data.table$var0)], xlab="Number in control", ylab="Estimated variance for control units", xaxt = "n")
  axis(1, at = 1:max(data.table$Num_Ctrl[!is.na(data.table$var0)]))
  par(mfrow = c(1,2))
  var1_plot
  var0_plot
}

#' Matched pairs type variance function
#'
#' Function that calculates matched pairs type variance.
#' @param data.table data frame containing info for set of blocks of same size
#' @param weighted indicates whether variance should be weighted by number of units
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
paired_var<-function(data.table, weighted=TRUE){
  #Vector of trt effect estimates
  tau_vec<-data.table$Y1-data.table$Y0
  var_est<-sum((tau_vec-mean(tau_vec))^2)/((length(tau_vec)-1)*length(tau_vec))
  if(weighted){
    var_est<-var_est*sum(data.table$nk)^2
  }
  return(var_est)
}

#' Aggregated matched pairs type variance estimator
#'
#' Function to calculate hybrid_m variance estimator.
#' @param data.small data frame containg info for small blocks
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
hybrid_m<-function(data.small){
  #First split into groups of blocks of same size
  by_size<-split(data.small, data.small$nk, drop = FALSE)
  #Estimate variance in each of these groups
  var_est<-sum(sapply(by_size, paired_var))/sum(data.small$nk)^2
  if ( is.nan( var_est ) ) {
      var_est = NA
  }
  return(var_est)
}

#' Pooled matched pairs type variance estimator
#'
#' Function to calculate hybrid_p variance estimator.
#' @param data.small data frame containg info for small blocks
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
hybrid_p<-function(data.small){
  #First check we can use this estimator
  n<-sum(data.small$nk)
  if(max(data.small$nk)>=n/2){
    return( NA )
  }
  #Vector of treatment effect estimates
  tau_vec<-data.small$Y1-data.small$Y0
  #Constant that adjusts for varying block sizes
  constant<-data.small$nk^2/((n-2*data.small$nk)*(n+sum(data.small$nk^2/(n-2*data.small$nk))))
  var_est<-sum(constant*(tau_vec-sum(tau_vec*data.small$nk)/n)^2)
  return(var_est)
}

#' Plug in type variance estimator
#'
#'Function to calculate variance estimate using average big block variances as plug-ins.
#' @param data.small data frame containg info for small blocks
#' @param data.big data frame containg info for small blocks
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
plug_in_big<-function(data.small, data.big){
  n<-sum(data.small$nk)
  #Get average variances in big blocks
  var0_avg<-mean(data.big$var0)
  var1_avg<-mean(data.big$var1)
  #Put in average variances in missing pieces of small blocks
  var.small.c<-as.numeric(data.small$Num_Ctrl==1)*var0_avg
  var.small.t<-as.numeric(data.small$Num_Trt==1)*var1_avg
  #Not all small blocks are missing things
  var.big.c<-as.numeric(data.small$Num_Ctrl>1)*data.small$var0/data.small$Num_Ctrl
  var.big.t<-as.numeric(data.small$Num_Trt>1)*data.small$var1/data.small$Num_Trt
  var.big.c[is.na(var.big.c)]<-0
  var.big.t[is.na(var.big.t)]<-0
  #Combine all variance terms
  var.c<-(var.small.c+var.big.c)*data.small$nk^2
  var.t<-(var.small.t+var.big.t)*data.small$nk^2
  var_est<-sum(var.c+var.t)/n^2
  return(var_est)
}

#' Block variance estimation function.
#'
#' This function takes the block-level summary statistics of a dataset and returns a treatment effect estimate,
#' variance estimate and some summary info about the block structure.
#'
#' @param data.table  Summary statistics of blocks in a dataset.  In particular, the output of block.data method.
#' @param method is method of variance estimation, defauly "hybrid_m"
#' @param throw.warnings TRUE means throw warnings if the hybrid estimators are breaking down due to violation of assumptions.
#' @export
fitdata.sumtable = function( data.table, method=c("hybrid_m", "hybrid_p", "plug_in_big"), throw.warnings=TRUE ) {
    K<-nrow(data.table)
    data.table$nk<-data.table$Num_Trt+data.table$Num_Ctrl
    n<-sum(data.table$nk)

    #Split into big and small blocks
    data.big<-data.table[data.table$Num_Trt>1&data.table$Num_Ctrl>1,]
    data.small<-data.table[data.table$Num_Trt==1|data.table$Num_Ctrl==1,]
    #Get variance for big blocks
    var_big<-sum(data.big$se_ney^2*(data.big$nk)^2)/n^2
    #If no small blocks, set small block variance to 0
    if(nrow(data.small)==0){
        var_small=0
    }
    #Estimate variance according to methof given
    else if(method=="hybrid_m"){
        var_small<-hybrid_m(data.small)*sum(data.small$nk)^2/n^2
        if(throw.warnings && is.na(var_small)){
            warning("Need multiple small blocks of the same size for hybrid_m")
        }
    }
    else if(method=="hybrid_p"){
        var_small<-hybrid_p(data.small)*sum(data.small$nk)^2/n^2
        if(throw.warnings && is.na(var_small)){
            warning("Largest small block is more than half the small block units in hybrid_p so no variance calc possible")
        }
    }
    else if(method=="plug_in_big"){
        var_small<-plug_in_big(data.small, data.big)*sum(data.small$nk)^2/n^2
        if(throw.warnings && is.na(var_small)){
            warning("Need some big blocks for plug_in_big")
        }
    }
    #Get trt effect estimates and aggregate
    tau_vec<-data.table$Y1-data.table$Y0
    tau_est<-sum(tau_vec*data.table$nk)/n
    #Get overall variance estimate
    var_est<-var_small+var_big
    #Table summarizing block sizes
    size_table<-data.table[,c("Blk_ID", "Num_Trt", "Num_Ctrl")]
    #Percent of blocks that are small
    perc_small<-round(sum(data.small$nk)/n*100,2)
    return_val<-list(tau_est, var_est, perc_small, size_table)
    names(return_val)<-c("tau_est", "var_est", "percent_small_blocks", "block_sizes")
    return_val$method = method
    return_val$se_est = sqrt( var_est )
    class(return_val)<-"var_dat"
    return(return_val)
}

#' Block variance estimation function.
#'
#' This function takes observed data and returns a treatment effect estimate,
#' variance estimate and some summary info about the block structure.
#'
#' @param Y vector observed outcomes
#' @param Z vector of assignment indicators (1==treated)
#' @param B block ids
#' @param data matrix of Y, Z, B, as alternative to using vectors
#' @param method is method of variance estimation, defauly "hybrid_m"
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
fitdata<-function(Y, Z, B, data=NULL, method=c("hybrid_m", "hybrid_p", "plug_in_big"), throw.warnings=TRUE ){
  if(!is.null(data)){
    Y<-data[,1]
    Z<-data[,2]
    B<-data[,3]
  }
  n<-length(Y)
  #Quick test that input is correct
  if(is.numeric(Z)==F){
    return("Treatment indicator should be vector of ones and zeros")
  }
  if((sum(Z==1)+sum(Z==0))!=n){
    return("Treatment indicator should be vector of ones and zeros")
  }
  #Get summary info
  data.table<-block.data(Y,Z,B)

  fitdata.sumtable( data.table, method=method,throw.warnings=throw.warnings)
}



#' Print method for fitdata output
#'
#' Function that compares difference variance estimators for blocked designs.
#' @param x output from fitdata (class=var_data)
#' @param ... further arguments to match print
#' @importFrom stats aggregate lm quantile rnorm sd var
#' @export
print.var_dat<-function(x, ...){
    scat("Randomization Inference Treatment Estimate (method = %s)\n", x$method )
    SE = sqrt(x$var_est)
    scat("Estimate of tau: %.2f \t(SE: %.2f)\n", x$tau_est, SE )
    scat("Nominal 95 Confidence Interval: %.2f - %.2f\n",  x$tau_est - 2*SE, x$tau_est + 2*SE )
    cat(paste(x$percent_small_blocks, "% of units are in small blocks", sep=""), "\n")
    cat("\nBlock Sizes:\n")
    print(x$block_sizes)
}


#' Estimate a series of linear models using different weighting schemes and
#' standard errors.
#'
#' @importFrom dplyr group_by ungroup mutate
#' @importFrom sandwich vcovHC vcovCL
#' @importFrom stats coef
#' @return Data frame of the various results.
linear.model.estimators = function( Y, Z, B, data=NULL ) {
    require( multiwayvcov )
    require( tidyverse )

    if ( !is.null( data ) ) {
        dat = data
        names( dat )[1:3] = c( "Yobs", "Z", "blk" )
    }  else {
        dat = data.frame( Yobs = Y, Z =Z, blk= B )
    }
    dat$blk<-as.factor(dat$blk)
    nj = table( dat$blk )

    # simple linear model
    M0 = lm( Yobs ~ Z + blk, data=dat )
    SE.lm = summary( M0 )$coeff["Z",2]

    # Huber-White SEs
    vcov_sand = sandwich::vcovHC(M0, type = "HC1")
    SE.lm.sand <-sqrt( vcov_sand[2,2] )

    # Cluster robust SEs
    #vcov_id <- multiwayvcov::cluster.vcov(M0, B)
    vcov_clust = sandwich::vcovCL( M0, B )
    SE.lm.clust = sqrt( vcov_clust[2,2] )

    FEmodels = data.frame( method=c("fixed effects (no inter)", "fixed effects (sand)", "fixed effects (cluster)" ),
                tau = c( coef(M0)[["Z"]], coef(M0)[["Z"]], coef(M0)[["Z"]] ),
                SE = c( SE.lm, SE.lm.sand, SE.lm.clust ), stringsAsFactors = FALSE )


    # Individual Weighted regression models
    Z.bar = mean( dat$Z )
    n = nrow( dat )
    J = length( unique( dat$blk ) )
    n.bar = n / J
    dat = dat %>% dplyr::group_by( blk ) %>%
        dplyr::mutate( p = mean( Z ),
                nj = n(),
                w.orig = ifelse( Z, 1/p, 1/(1-p) ),
                weight = ifelse( Z, Z.bar / p, (1-Z.bar)/(1-p) ),
                weight.site = weight * n.bar / nj ) %>%
        dplyr::ungroup()

    M0w = lm( Yobs ~ Z + blk, weights=dat$w.orig, data=dat )
    SE.w = summary( M0w )$coeff["Z",2]

    M0w2 = lm( Yobs ~ Z + blk, weights=dat$weight, data=dat )
    SE.w2 = summary( M0w2 )$coeff["Z",2]

    # Site weighted regression models
    M0w.site = lm( Yobs ~ Z + blk, weights=dat$weight.site, data=dat )
    tau.w.site = coef( M0w.site )[[2]]
    SE.w.site = summary( M0w.site )$coeff["Z",2]

    weightModels = data.frame( method=c("IPTW weighted regression (naive)", "IPTW-weighted regression", "IPTW-weighted regression (site)"),
                         tau = c( coef( M0w )[["Z"]], coef( M0w2 )[["Z"]], coef( M0w.site )[["Z"]] ),
                         SE = c( SE.w, SE.w2, SE.w.site ), stringsAsFactors = FALSE )


    # Interacted linear regression models
    M0.int = lm( Yobs ~ 0 + Z * blk - Z, data = dat )
    VC = vcov( M0.int )

    tau.hats = coef(M0.int)[J + 1:J]

    wts = c( rep( 0, J ), rep( 1/J, J ) )
    tau.site = mean( tau.hats )
    SE.site = sqrt( t(wts) %*% VC %*% wts )

    wts.indiv = c( rep( 0, J ), nj/n )
    tau.indiv = weighted.mean( tau.hats, nj )
    SE.indiv = sqrt( t(wts.indiv) %*% VC %*% wts.indiv )

    interactModels = data.frame( method=c("fixed effects interact (site)", "fixed effects-interact (indiv)"),
                                 tau = c( tau.site, tau.indiv ),
                                 SE = c( SE.site, SE.indiv ),
                                 stringsAsFactors = FALSE)
    # combine and return results
    bind_rows( FEmodels, weightModels, interactModels )
}

if ( FALSE ) {
    dat = make.obs.data(p = 0.2)
    linear.model.estimators( dat$Yobs, dat$Z, dat$blk )
}



