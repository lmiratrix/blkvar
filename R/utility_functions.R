##
## Utility functions for the rest of the package
##


scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


##### Functions to help process control functions  ####


# #' Make a canonical fixed effect formula, possibly with control variables.
# #'
# #' @param Yobs Name of outcome variable (assumed to exist in data)
# #' @param Z Name of treatment variable (assumed to exist in data)
# #' @param B Name of blocking variable (assumed to exist in data)
# #' @param control.formula What variables to control for, in the form of "~ X1 + X2".
# #' @param data Dataframe holding all variables to be used in formula.
# #'
# #' @return Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X1 + X2"

make_base_formula = function( Yobs = "Yobs", Z = "Z", control.formula = NULL, data = NULL) {
  if (is.null(control.formula)) {
    new.form <- sprintf( "%s ~ 1 + %s ", Yobs, Z )
    return(as.formula(new.form ))
  }

  if (length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1) {
    stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
  }

  if (!is.null(data)) {
    control.vars <- formula.tools::get.vars(control.formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
      stop("Some variables in control.formula are not present in your data.")
    }
  }

  c.names <- formula.tools::rhs.vars(control.formula)
  new.form <- sprintf( "%s ~ 1 + %s + %s", Yobs, Z, paste( c.names, collapse =" + " ))
  return(as.formula( new.form ))
}

# #' Make a canonical fixed effect formula, possibly with control variables.
# #'
# #' @inheritParams make_base_formula
# #'
# #' @return Something like "Yobs ~ 0 + Z + B" or "Yobs ~ 0 + Z + B + X"

make_FE_formula <- function(Yobs = "Yobs", Z = "Z", B = "B", control.formula = NULL, data = NULL) {
  if ( is.null( control.formula)) {
    new.form <- sprintf( "%s ~ 0 + %s + %s", Yobs, Z, B)
    return( as.formula(new.form))
  }

  if (length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1) {
    stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
  }

  if (!is.null(data)) {
    control.vars <- formula.tools::get.vars(control.formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
      stop("Some variables in control.formula are not present in your data.")
    }
  }
  c.names <- formula.tools::rhs.vars(control.formula)
  new.form <- sprintf( "%s ~ 0 + %s + %s + %s", Yobs, Z, B, paste(c.names, collapse =" + " ))
  return(as.formula(new.form))
}

# #' Make a canonical fixed effect formula, possibly with control variables.
# #'
# #' @inheritParams make_base_formula
# #'
# #' @return Something like "Yobs ~ 0 + Z * B - Z" or "Yobs ~ 0 + Z * B - Z + X"

make_FE_int_formula <- function(Yobs = "Yobs", Z = "Z", B = "B", control.formula = NULL, data = NULL) {
  if (is.null(control.formula)) {
    new.form <- sprintf( "%s ~ 0 + %s * %s - %s", Yobs, Z, B, Z)
    return(as.formula(new.form))
  }
  if (length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1) {
    stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
  }
  if (!is.null(data)) {
    control.vars <- formula.tools::get.vars(control.formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
      stop("Some variables in control.formula are not present in your data.")
    }
  }
  c.names <- formula.tools::rhs.vars(control.formula)
  new.form = sprintf( "%s ~ 0 + %s * %s - %s + %s", Yobs, Z, B, Z, paste( c.names, collapse =" + " ))
  return(as.formula(new.form))
}

# #' Check formula and control.formula for syntax
# #'
# #' Checks syntax and also presence of variables in the dataframe.
# #'
# #' Then makes the outcome, treatment and block variables have canonical "Yobs", "Z", and "B" names.
# #'
# #' @return Dataset with now-canonical variable names
make_canonical_data <- function(formula, control.formula = NULL, siteID = NULL, data) {
   # Figure out the covariates we are using
    if (length(formula.tools::lhs.vars(formula)) != 1 | length(formula.tools::rhs.vars(formula)) != 2) {
      stop("The formula argument must be of the form outcome ~ treatment:block_id.")
    }
    main.vars <- formula.tools::get.vars(formula, data = data)
    if (any(!(main.vars %in% colnames(data)))) {
      stop("Some variables in formula are not present in your data.")
    }
    if (!is.null(control.formula)) {
      if(length(formula.tools::lhs.vars(control.formula)) != 0 | length(formula.tools::rhs.vars(control.formula)) < 1) {
        stop("The control formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
      }
      control.vars <- formula.tools::get.vars(control.formula, data = data)
      if (any(!(control.vars %in% colnames(data)))) {
        stop("Some variables in control.formula are not present in your data.")
      }
    }
    out.name <- formula.tools::lhs.vars(formula)[[1]]
    main.name <- formula.tools::rhs.vars(formula)

    # Copy over the variables to our names
    data$Yobs <- data[[out.name]]
    data$Z <- data[[ main.name[[1]] ]]

    if ( length( unique( data$Z ) ) != 2 ) {
        stop( sprintf( "Identified treatment variable '%s' has more than two values. Did you swap treatment and block?",
                       main.name[[1]] ) )
    }

    data$B <- data[[ main.name[[2]] ]]

    if (is.null(siteID)) {
      data$siteID <- data$B
    } else {
      data$siteID <- data[[siteID]]
    }
   return( data )
}

# # #
# if ( FALSE ) {

    # make_base_formula()
    # make_base_formula( control.formula = ~ X1 + X2 + X3 )
    # make_FE_formula(  )
    # make_FE_formula( control.formula = ~ X1 + X2 + X3 )

    # make_FE_int_formula( )
    # make_FE_int_formula( control.formula = ~ X1 + X2 + X3 )

# }
