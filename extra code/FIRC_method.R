##
## FIRC model methods
##
## This file has implementation of the FIRC model (both with pooled variances
## across Tx and Co units and unpooled where Tx and Co get their own variance
## terms)
##
## The unpooled functions are taken from Catherine's implementation of Weiss et
## al methods.
##

# # For the debugging code to get the DGP files
# localsource = function( filename ) {
    # source( file.path( dirname( rstudioapi::getActiveDocumentContext()$path ), filename ) )
# }


# # Testing
# if ( FALSE ) {

    # localsource( "multisite_data_generators.R")

    # dat = catherine_gen_dat( 0.2, 1.0, 30, 50 )
    # head( dat )
    # describe_data( dat )


    # #debug( estimate_ATE_FIRC )
    # estimate_ATE_FIRC( Yobs, Z, sid, data=dat )

    # dat$X = dat$Y0 + rnorm( nrow(dat) )
    # estimate_ATE_FIRC( Yobs, Z, sid, data=dat, control.formula = ~ X )

    # dat = catherine_gen_dat( 0.2, 0, 30, 50 )
    # head( dat )
    # describe_data( dat )

    # #debug( estimate_ATE_FIRC )
    # estimate_ATE_FIRC( Yobs, Z, B, dat )
# }


# if ( FALSE ) {

    # localsource( "multisite_data_generators.R")

    # dat = catherine_gen_dat( 0.2, 0.0, 30, 50 )
    # describe_data( dat )
    # head( dat )

    # estimate_ATE_FIRC( Yobs, Z, sid, data=dat, pool=TRUE )

    # dat = catherine_gen_dat( 0.2, 0.5, 30, 50 )
    # describe_data( dat )

    # estimate_ATE_FIRC( Yobs, Z, B, data=dat, pool=TRUE )

# }


# # Testing
# if ( FALSE ) {

    # df = gen_dat_no_cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

    # # head( df )
    # # describe_data( df )

    # analysis.FIRC( Yobs, Z, B, df )

# }

# # Testing
# # This testing is based on the DGP for small sample simulations
# if ( FALSE ) {

  # df = gen_dat_no_cov.n( n = 600,n.small = 6, J = 30, small.percentage = 0.7,tau.11.star = 0.2)

  # # head( df )
  # # describe_data( df )

  # analysis.FIRC( df )

# }