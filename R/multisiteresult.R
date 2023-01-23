



scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


#' Pretty print result from Q-statistic call
#'
#' @export
#' @param x A multisiteresult object.
#' @param ... No extra options passed.
#' @family multisiteresult
print.multisiteresult = function( x, ... ) {
    args = attr( x, "args" )

    anal_type = args$model

    suff = ""
    if ( !is.null( args$method ) ) {
        suff = paste0( " (", args$method, ")" )
    }
    scat( "Cross-site %s analysis%s\n", anal_type, suff )
    scat( "  Grand ATE: %s", format( x$ATE_hat ) )
    if ( !is.null( x$SE_ATE ) ) {
        scat( " (%s)\n", format( x$SE_ATE ) )
    } else {
        scat( "\n" )
    }

    scat( "  Cross-site SD: %s",
          format( x$tau_hat ) )
    if ( is.null( x$CI_low ) && is.null( x$CI_high ) ) {
        if ( !is.null( x$SE_tau ) ) {
            scat( " (%s)\n", format( x$SE_tau ) )
        } else {
            cat("\n" )
        }
    } else {
        scat( " (%s%% CI %s--%s)\n",
              format( 100*(1-args$alpha*2) ),
              format( x$CI_low ), format( x$CI_high ) )
    }

    stat_name = "deviance"
    stat_val = ""
    if ( args$model == "Q-statistic" ) {
        stat_name = "Q"
        stat_val = format( x$Q )
    } else {
        stat_val = format( x$deviance )
    }

    scat( "    %s = %s, p=%.4f\n", stat_name, stat_val, x$p_variation )

    if ( args$model == "Q-statistic" ) {
        est_method = ifelse( args$calc_optim_pvalue, "max p-value", "MoM" )
        scat( "  [ tolerance = %s; mean method = %s; est method = %s ]\n",
          format( args$tolerance, digits=2 ), args$mean_method,
          est_method )
    }

    invisible( x )
}



#' Is object a multisiteresult object?
#'
#' @export
#' @aliases multisiteresult
#' @param x the object to check.
#' @family multisiteresult
is.multisiteresult = function( x ) {
    inherits(x, "multisiteresult")
}




#' Cast qstat result to data.frame
#'
#' @export
#' @aliases multisiteresult
#' @param x the multisiteresult object to covert
#' @family multisiteresult
as.data.frame.multisiteresult = function( x ) {
    class(x) = "list"
    as.data.frame( x )
}

