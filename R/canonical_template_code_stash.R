
# Canonical code to setup the weird parameter passing


if ( FALSE ) {

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



}
