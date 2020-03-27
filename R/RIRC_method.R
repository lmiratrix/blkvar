#
# This file contains functions implementing the various methods for detecting
#  and estimating treatment effect variation across sites in a multi-site trial
# using random-slope, random-intercept approaches.
#

# We have:
#   * Ideosyncratic test
#   * Testing for systematic covariate (two versions, one with random coeff and one not)
#   * Hybrid test doing a likelihood ratio of the hybrid model vs. no-variation model
#









# localsource <- function( filename ) {
    # source( file.path( dirname( rstudioapi::getActiveDocumentContext()$path ), filename ) )
# }


# if ( FALSE ) {
    # localsource( "MLM_method.R" )
    # localsource( "multisite_data_generators.R" )

    # dat <- catherine_gen_dat( 0.2, 5, 30, 50 )
    # head( dat )
    # describe_data( dat, Y0 <- "Y0", Y1="Y1" )


    # estimate_ATE_RIRC( Yobs, Z, B, data=dat )

    # dat <- catherine_gen_dat( 0.2, 0, 30, 50 )
    # describe_data( dat, Y0 <- "Y0", Y1="Y1" )

    # estimate_ATE_RIRC( Yobs, Z, B, data=dat )
# }






##
## Testing code
##

# if ( FALSE ) {
    # source( "R/multisite_data_generators.R")
    # dat <- catherine_gen_dat( 0.2, 0.2, 30, 50 )
    # head( dat )
    # #analysis_idio.RIRC.pool( Yobs, Z, B, data=dat )
    # #fit.RIRC.pool( Yobs, Z, B, data=dat )

    # estimate_ATE_RIRC( Yobs, Z, sid, data=dat )
    # estimate_ATE_RIRC_pool( Yobs, Z, sid, data=dat )

    # #    fit.RIRC( dat$Yobs, dat$Z, dat$sid )
# }

