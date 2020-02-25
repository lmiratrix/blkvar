
##
## Functions to help process control functions
##

library( formula.tools )

#' Make a canonical fixed effect formula, possibly with control variables.
#'
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z Name of treatment variable (assumed to exist in data)
#' @param B Name of blocking variable (assumed to exist in data)
#' @param control.formula What variables to control for, in the form of "~ X1 + X2".
#' @param data Dataframe holding all variables to be used in formula.
#'
#' @return Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X"
#' @keywords internal
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
#' @inheritParams make.base.formula
#'
#' @return Something like "Yobs ~ 0 + Z + B" or "Yobs ~ 0 + Z + B + X"
#' @keywords internal
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
#' @inheritParams make.base.formula
#'
#' @return Something like "Yobs ~ 0 + Z * B - Z" or "Yobs ~ 0 + Z * B - Z + X"
#' @keywords internal
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
