test_that("compare_methods_oracle works", {

    dat = make_data_linear()
    head( dat )
    dat$blk = form_blocks_from_continuous( dat$X, method="small" )
    data = dat
    head( data )
    table( dat$blk )
    p_mat<- dat %>% group_by( blk ) %>%
        summarise( p = round( 0.5 * n() ) / n() ) %>%
        as.data.frame()
    p_mat = rename( p_mat, B = blk )
    p_mat

    v1 <- compare_methods_oracle( Y0, Y1, blk, data=data, p_mat = p_mat )

    v1
    expect_equal( nrow(v1), 5 )
 })
