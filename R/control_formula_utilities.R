
##
## Functions to help process control functions
##

library( formula.tools )

#' Make a canonical fixed effect formula, possibly with control variables.
#'
#' @param Yobs
#' @param Z
#' @param B
#' @param control.formula
#' @param data
#'
#' @return Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X"
make.base.formula = function( Yobs = "Yobs", Z = "Z", control.formula = NULL, data = NULL) {

    if ( is.null( control.formula ) ) {
        new.form = sprintf( "%s ~ 1 + %s ",
                            Yobs, Z )
        return( as.formula( new.form ) )
    }


    if(length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1){
        stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
    }

    if (!is.null( data ) ) {
        control.vars <- formula.tools::get.vars(control.formula, data=data)
        if(any(!(control.vars %in% colnames(data)))){
            stop("Some variables in control.formula are not present in your data.")
        }
    }

    c.names = formula.tools::rhs.vars(control.formula)
    new.form = sprintf( "%s ~ 1 + %s + %s",
                        Yobs, Z,
                        paste( c.names, collapse =" + " ) )
    return( as.formula( new.form ) )

}



#' Make a canonical fixed effect formula, possibly with control variables.
#'
#' @param Yobs
#' @param Z
#' @param B
#' @param control.formula
#' @param data
#'
#' @return Something like "Yobs ~ 0 + Z + B" or "Yobs ~ 0 + Z + B + X"
make.FE.formula = function( Yobs = "Yobs", Z = "Z", B = "B", control.formula = NULL, data = NULL) {

    if ( is.null( control.formula ) ) {
        new.form = sprintf( "%s ~ 0 + %s + %s",
                            Yobs, Z, B )
        return( as.formula( new.form ) )
    }

    if(length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1){
        stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
    }

    if (!is.null( data ) ) {
         control.vars <- formula.tools::get.vars(control.formula, data=data)
        if(any(!(control.vars %in% colnames(data)))){
            stop("Some variables in control.formula are not present in your data.")
        }
    }

    c.names = formula.tools::rhs.vars(control.formula)
    new.form = sprintf( "%s ~ 0 + %s + %s + %s",
                        Yobs, Z, B,
                        paste( c.names, collapse =" + " ) )
    return( as.formula( new.form ) )

}

#' Make a canonical fixed effect formula, possibly with control variables.
#'
#' @param Yobs
#' @param Z
#' @param B
#' @param control.formula
#' @param data
#'
#' @return Something like "Yobs ~ 0 + Z * B - Z" or "Yobs ~ 0 + Z * B - Z + X"
make.FE.int.formula = function( Yobs = "Yobs", Z = "Z", B = "B", control.formula = NULL, data = NULL) {

    if ( is.null( control.formula ) ) {
        new.form = sprintf( "%s ~ 0 + %s * %s - %s",
                            Yobs, Z, B, Z )
        return( as.formula( new.form ) )
    }

    if(length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1){
        stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
    }

    if (!is.null( data ) ) {
        control.vars <- formula.tools::get.vars(control.formula, data=data)
        if(any(!(control.vars %in% colnames(data)))){
            stop("Some variables in control.formula are not present in your data.")
        }
    }

    c.names = formula.tools::rhs.vars(control.formula)
    new.form = sprintf( "%s ~ 0 + %s * %s - %s + %s",
                        Yobs, Z, B, Z,
                        paste( c.names, collapse =" + " ) )
    return( as.formula( new.form ) )

}



if ( FALSE ) {

    make.base.formula()
    make.base.formula( control.formula = ~ X1 + X2 + X3 )
    make.FE.formula(  )
    make.FE.formula( control.formula = ~ X1 + X2 + X3 )

    make.FE.int.formula( )
    make.FE.int.formula( control.formula = ~ X1 + X2 + X3 )

}
