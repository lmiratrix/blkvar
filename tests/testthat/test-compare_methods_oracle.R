
test_that( "calc_summary_stats_oracle works", {

    dat = generate_blocked_data_linear()
    head( dat )
    nrow( dat )
    dat$blk = form_blocks_from_continuous( dat$X, method="small" )
    dat$Tx = rep( c(0,1), nrow(dat)/2 )

    p_mat<- dat %>% group_by( blk ) %>%
        dplyr::summarise( p = mean( Tx ) ) %>%
        as.data.frame()
    p_mat

    dt = calc_summary_stats_oracle( Y0, Y1, blk, data=dat, p_mat = p_mat )
    dt
    expect_equal( nrow( dt ), 5 )
    expect_equal( sum( dt$n1 + dt$n0 ), nrow( dat ) )
    expect_equal( sum( dt$n ), nrow( dat ) )
    expect_equal( dt$B, sort( unique( as.character( dat$blk ) ) ) )

    dt2 = calc_summary_stats_oracle( dat$Y0, dat$Y1, dat$blk, p_mat = p_mat )
    dt2
    expect_equal( dt, dt2 )

    dt3 = calc_summary_stats_oracle( Y0, Y1, blk, data=dat, Z=Tx )
    expect_equal( dt, dt3 )

})




test_that("compare_methods_oracle works", {

    dat = generate_blocked_data_linear()
    head( dat )
    dat$blk = form_blocks_from_continuous( dat$X, method="small" )
    p_mat<- dat %>% group_by( blk ) %>%
        dplyr::summarise( p = round( 0.5 * n() ) / n() ) %>%
        as.data.frame()
    p_mat


    datBig = bind_rows( dat, dat )
    head( datBig )
    table( datBig$blk )

    stats = calc_summary_stats_oracle( datBig$Y0, datBig$Y1, datBig$blk, p_mat=p_mat )
    stats

    v1 <- compare_methods_oracle( Y0, Y1, blk, data=datBig, p_mat = p_mat )
    v1
    expect_equal( nrow(v1), 5 )
    expect_true( all( !is.na( v1$SE ) ) )

    # Now with tiny blocks
    v2 <- compare_methods_oracle( Y0, Y1, blk, data=dat, p_mat = p_mat )
    v2
    expect_equal( nrow(v2), 5 )
 })
