
##
## This file contains the various estimators discussed in Pashley et al
##


#' Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}





#' Summarise data by block.
#'
#' Given dataframe or list of variables, return table of stats for each
#' randomization block.
#'
#' @param dat Dataframe with defined Yobs, Z, and B variables.
#' @param siteID If not null, name of siteID that has randomization blocks
#'   nested inside.
#'
#' @return dataframe with summary statistics by block
#' @export
calc.summary.stats = function( Yobs, Z, B, data = NULL, siteID = NULL, add.neyman = FALSE ) {
    require( tidyverse )
    if ( missing( "Z" ) && is.null( data ) ) {
        data = Yobs
    }
    if ( is.null( data ) ) {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
    } else {
        if ( is.matrix( data ) ) {
            data = as.data.frame( data )
        }
        if ( !is.null( siteID ) && (length( siteID ) == 1 ) ) {
            siteID = data[, siteID ]
        }

        if ( missing( "Z" ) ) {
            stopifnot( all( c( "Yobs", "Z", "B" ) %in% names(data) ) )
        } else {
            data = rename_( data,
                            Yobs = as.character( substitute( Yobs ) ),
                            Z = as.character( substitute( Z ) ),
                            B = as.character( substitute( B ) ) )
        }
    }
    dat = data
    if ( !is.null( siteID ) ) {
        dat$siteID = siteID
    }

    if ( !is.null( siteID ) ) {
        sdat <- dat %>% dplyr::group_by( B, Z, siteID )
    } else {
        sdat <- dat %>% dplyr::group_by( B, Z )
    }

    sdat <- sdat %>%
        dplyr::summarise( n = n(),
                   Ybar = mean( Yobs ),
                   var = var( Yobs ) )
    sdat <- reshape( as.data.frame(sdat), direction="wide",
                     v.names = c("n","Ybar", "var"),
                     idvar = "B",
                     timevar = "Z",
                     sep="" )
    sdat$n = sdat$n0 + sdat$n1

    if ( add.neyman ) {
        sdat = mutate( sdat,
                       se_ney = sqrt( var1 / n1 + var0 / n0 ) )
    }

    if ( any( is.na( sdat$n1 ) | is.na( sdat$n0 ) ) ) {
        #sdat$n1[ is.na(sdat$n1) ] = 0
        #sdat$n0[ is.na(sdat$n0) ] = 0
        sdat = filter( sdat, !is.na( n1 ), !is.na( n0 ) )
        warning( "Some blocks have no treatment or no control units; they are being dropped" )
    }

    sdat

}




# convert.table.names = function( data.table ) {
#     names( data.table ) = c( "B", "n1", "n0", "Ybar1", "Ybar0", "var.1", "var.0", "se_ney" )
#     data.table$n = with( data.table, n0 + n1 )
#     data.table
# }

if ( FALSE ) {
    dat = make.data( c( 3,4,5 ) )
    dat = add.obs.data(dat)

    calc.summary.stats(dat)
}


#' Plot diagnostic function
#'
#' Function that plots variance estimates versus size of treatment group.
#'
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
  data.table<-calc.summary.stats(Y,Z,B)
  #data.table$nk<-data.table$n1+data.table$n0
  var1_plot<-plot(data.table$n1[!is.na(data.table$var1)], data.table$var1[!is.na(data.table$var1)], xlab="Number treated", ylab="Estimated variance for treated units", xaxt = "n")
  axis(1, at = 1:max(data.table$n1[!is.na(data.table$var1)]))
  var0_plot<-plot(data.table$n0[!is.na(data.table$var0)], data.table$var0[!is.na(data.table$var0)], xlab="Number in control", ylab="Estimated variance for control units", xaxt = "n")
  axis(1, at = 1:max(data.table$n0[!is.na(data.table$var0)]))
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
  tau_vec<-data.table$Ybar1-data.table$Ybar0
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
  tau_vec<-data.small$Ybar1-data.small$Ybar0
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
  var.small.c<-as.numeric(data.small$n0==1)*var0_avg
  var.small.t<-as.numeric(data.small$n1==1)*var1_avg
  #Not all small blocks are missing things
  var.big.c<-as.numeric(data.small$n0>1)*data.small$var0/data.small$n0
  var.big.t<-as.numeric(data.small$n1>1)*data.small$var1/data.small$n1
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
#' @param data.table  Summary statistics of all the blocks in a dataset.  In particular, this is the output of calc.summary.stats method.
#' @param method The method to be used for variance estimation, defauly "hybrid_m"
#' @param throw.warnings TRUE means throw warnings if the hybrid estimators are breaking down due to violation of assumptions.
#' @export
fitdata.sumtable = function( data.table,
                             method=c("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small"), throw.warnings=TRUE ) {
    method = match.arg( method )
    K<-nrow(data.table)
    data.table$nk<-data.table$n1+data.table$n0
    n<-sum(data.table$nk)

    #Split into big and small blocks
    data.big<-data.table[data.table$n1>1&data.table$n0>1,]
    data.small<-data.table[data.table$n1==1|data.table$n0==1,]
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
    }else if(method=="rct_yes_small"){
      dt<-convert.table.names(data.small)
      var_small<-(calc.RCT.Yes.SE(dt, "superpop.original" , weight = "individual")$SE)^2*sum(data.small$nk)^2/n^2
      if(throw.warnings && is.na(var_small)){
        warning("Need multiple small blocks for rct_yes_small")
      }
    }else if(method=="rct_yes_mod_small"){
      dt<-convert.table.names(data.small)
      var_small<-(calc.RCT.Yes.SE(dt, "superpop" , weight = "individual")$SE)^2*sum(data.small$nk)^2/n^2
      if(throw.warnings && is.na(var_small)){
        warning("Need multiple small blocks for rct_yes_small")
      }
    }
    #Get trt effect estimates and aggregate
    tau_vec<-data.table$Ybar1-data.table$Ybar0
    tau_est<-sum(tau_vec*data.table$nk)/n
    #Get overall variance estimate
    if(method=="rct_yes_all"){
      dt<-convert.table.names(data.table)
      var_est<-(calc.RCT.Yes.SE(dt, "superpop.original" , weight = "individual")$SE)^2
    }else if(method=="rct_yes_mod_all"){
      dt<-convert.table.names(data.table)
      var_est<-(calc.RCT.Yes.SE(dt, "superpop" , weight = "individual")$SE)^2
    }else{
      var_est<-var_small+var_big
    }
    #Table summarizing block sizes
    size_table<-data.table[,c("B", "n1", "n0")]
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
fitdata<-function(Y, Z, B, data=NULL,
                  method=c("hybrid_m", "hybrid_p", "plug_in_big", "rct_yes_all", "rct_yes_small", "rct_yes_mod_all", "rct_yes_mod_small"),
                  throw.warnings=TRUE ) {
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
  data.table<-calc.summary.stats(Y,Z,B)

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



