##
## Utility functions for the rest of the package
##


scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


##### Functions to help process control functions  ####


#' Make a canonical fixed effect formula, possibly with control variables.
#'
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z Name of treatment variable (assumed to exist in data)
#' @param B Name of blocking variable (assumed to exist in data)
#' @param control_formula What variables to control for, in the form of "~ X1 + X2".
#' @param data Dataframe holding all variables to be used in formula.
#'
#' @return Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X1 + X2"
#' @noRd
make_base_formula = function( Yobs = "Yobs", Z = "Z", control_formula = NULL, data = NULL) {
  if (is.null(control_formula)) {
    new.form <- sprintf( "%s ~ 1 + %s ", Yobs, Z )
    return(as.formula(new.form ))
  }

  if (length(formula.tools::lhs.vars(control_formula)) != 0 | length(formula.tools::rhs.vars(control_formula)) < 1) {
    stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
  }

  if (!is.null(data)) {
    control.vars <- formula.tools::get.vars(control_formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
      stop("Some variables in control_formula are not present in your data.")
    }
  }

  c.names <- formula.tools::rhs.vars(control_formula)
  new.form <- sprintf( "%s ~ 1 + %s + %s", Yobs, Z, paste( c.names, collapse =" + " ))
  return(as.formula( new.form ))
}

#' Make a canonical fixed effect formula, possibly with control variables.
#'
#' @inheritParams make_base_formula
#'
#' @return Something like "Yobs ~ 0 + Z + B" or "Yobs ~ 0 + Z + B + X"
#' @noRd
make_FE_formula <- function(Yobs = "Yobs", Z = "Z", B = "B", control_formula = NULL, data = NULL) {
  if ( is.null( control_formula)) {
    new.form <- sprintf( "%s ~ 0 + %s + %s", Yobs, Z, B)
    return( as.formula(new.form))
  }

  if (length(formula.tools::lhs.vars(control_formula)) != 0 | length(formula.tools::rhs.vars(control_formula)) < 1) {
    stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
  }

  if (!is.null(data)) {
    control.vars <- formula.tools::get.vars(control_formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
      stop("Some variables in control_formula are not present in your data.")
    }
  }
  c.names <- formula.tools::rhs.vars(control_formula)
  new.form <- sprintf( "%s ~ 0 + %s + %s + %s", Yobs, Z, B, paste(c.names, collapse =" + " ))
  return(as.formula(new.form))
}

# #' Make a canonical fixed effect formula, possibly with control variables.
# #'
# #' @inheritParams make_base_formula
# #'
# #' @return Something like "Yobs ~ 0 + Z * B - Z" or "Yobs ~ 0 + Z * B - Z + X"

make_FE_int_formula <- function(Yobs = "Yobs", Z = "Z", B = "B", control_formula = NULL, data = NULL) {
  if (is.null(control_formula)) {
    new.form <- sprintf( "%s ~ 0 + %s * %s - %s", Yobs, Z, B, Z)
    return(as.formula(new.form))
  }
  if (length(formula.tools::lhs.vars(control_formula)) != 0 | length(formula.tools::rhs.vars(control_formula)) < 1) {
    stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
  }
  if (!is.null(data)) {
    control.vars <- formula.tools::get.vars(control_formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
      stop("Some variables in control_formula are not present in your data.")
    }
  }
  c.names <- formula.tools::rhs.vars(control_formula)
  new.form = sprintf( "%s ~ 0 + %s * %s - %s + %s", Yobs, Z, B, Z,
                      paste( c.names, collapse =" + " ))
  return(as.formula(new.form))
}

#' Check formula and control_formula for syntax
#'
#' Checks syntax and also presence of variables in the dataframe.
#'
#' Then makes the outcome, treatment and block variables have canonical "Yobs",
#' "Z", and "B" names.
#'
#' Expand any factors, etc., in the control formula and make the corresponding
#' variables.
#'
#' @return Dataset with now-canonical variable names, and no extraneous
#'   variables.
#'
#' @noRd
make_canonical_data <- function(formula, control_formula = NULL, siteID = NULL, data ) {
    # Figure out the covariates we are using
    if (length(formula.tools::lhs.vars(formula)) != 1 | length(formula.tools::rhs.vars(formula)) != 2) {
      stop("The formula argument must be of the form outcome ~ treatment:block_id.")
    }
    main.vars <- formula.tools::get.vars(formula, data = data)
    if (any(!(main.vars %in% colnames(data)))) {
      stop("Some variables in formula are not present in your data.")
    }

    # Control variables?
    control.vars = c()
    if (!is.null(control_formula)) {
      if(length(formula.tools::lhs.vars(control_formula)) != 0 |
         length(formula.tools::rhs.vars(control_formula)) < 1) {
        stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
      }
      control.vars <- formula.tools::get.vars(control_formula, data = data)
      if (any(!(control.vars %in% colnames(data)))) {
        stop("Some variables in control_formula are not present in your data.")
      }
    }


    # Make canonical names for outcome, treatment, etc.
    out.name <- formula.tools::lhs.vars(formula)[[1]]
    main.name <- formula.tools::rhs.vars(formula)

    # Copy over the variables to our names
    data$Yobs <- data[[out.name]]
    data$Z <- data[[ main.name[[1]] ]]

    # Check for valid treatment variable
    if ( length( unique( data$Z ) ) != 2 ) {
        stop( sprintf( "Identified treatment variable '%s' has more than two values. Did you swap treatment and block?",
                       main.name[[1]] ) )
    }

    # Make blocking variable with canonical name
    data$B <- data[[ main.name[[2]] ]]

    # Add site variable (same as block if no RA blocks in site).
    if (is.null(siteID)) {
      data$siteID <- data$B
    } else {
      data$siteID <- data[[siteID]]
    }

   # drop extra variables
   data = data[ c( c("Yobs", "Z", "siteID", "B"), control.vars ) ]

   return( data )
}

# # #
# if ( FALSE ) {

    # make_base_formula()
    # make_base_formula( control_formula = ~ X1 + X2 + X3 )
    # make_FE_formula(  )
    # make_FE_formula( control_formula = ~ X1 + X2 + X3 )

    # make_FE_int_formula( )
    # make_FE_int_formula( control_formula = ~ X1 + X2 + X3 )

# }
