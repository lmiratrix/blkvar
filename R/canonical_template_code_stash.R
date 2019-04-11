
# Canonical code to setup the weird parameter passing
#
# This code needs to be pasted at the top of most of the methods in this
# package.  There should be a better way.
#
# It allows for unquoted naming of variables in the method calls.


if ( FALSE ) {


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
                stopifnot( !is.null( siteID ) )
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


    #
    # To just make the variables as vectors
    #
    if(!is.null(data)){
        if ( missing( "Yobs" ) ) {
            Yobs<-data[,1]
            Z<-data[,2]
            B<-data[,3]
        } else {
            Yobs = eval( substitute( Yobs ), data )
            Z = eval( substitute( Z ), data )
            B = eval( substitute( B ), data)
        }
    } else {
        if ( is.data.frame(Yobs) ) {
            stopifnot( is.null( data ) )
            B = Yobs$B
            Z = Yobs$Z
            Yobs = Yobs$Yobs
        }
    }


    #
    # To make a dataframe with canonical variable names
    #

    if( !is.null(data) ){
        if ( missing( "Yobs" ) ) {
            data = data.frame( Yobs = data[[1]],
                               Z = data[[2]],
                               B = data[[3]] )
        } else {
            if ( !is.null( siteID ) ) {
                siteID = data[[siteID]]
                stopifnot( !is.null( siteID ) )
            }
            data = data.frame( Yobs = eval( substitute( Yobs ), data ),
                               Z = eval( substitute( Z ), data ),
                               B = eval( substitute( B ), data) )
            if ( is.null( siteID ) ) {
                data$siteID = data$B
            } else {
                data$siteID = siteID
            }
        }
    } else {
        data = data.frame( Yobs = Yobs,
                           Z = Z,
                           B = B )
        if ( !is.null( siteID ) ) {
            data$siteID = siteID
        }
    }
    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )


    #
    # Checks on dataframe integrity
    #
    stopifnot( length( unique( data$Z ) ) == 2 )
    stopifnot( is.numeric( data$Yobs ) )





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
