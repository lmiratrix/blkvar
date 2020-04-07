#' Describe oracle dataset (for simulation studies)
#'
#' Utility function to describe different characteristics of the distribution of
#' sites (given all the potential outcomes).
#' @param data Dataframe with defined Y0, Y1, and B variables.
#' @param Y0 name of Y0 column
#' @param Y1 name of Y1 column
#' @param Z vector that indicates if outcome is under treatment or control
#' @param sid block ids
#' @export

describe_data <- function(data, Y0 = "Y0", Y1 = "Y1", Z = "Z", sid = "sid") {
  # old param list: Y0, Y1, Z, sid, data = NULL
  data <- dplyr::rename(data,
                 Y0 = !!rlang::sym(Y0),
                 Y1 = !!rlang::sym(Y1),
                 Z = !!rlang::sym(Z),
                 sid = !!rlang::sym(sid))
  sites <- data %>% group_by(!!rlang::sym(sid)) %>%
      dplyr::summarise(Y0.bar = mean(Y0),
                       Y1.bar = mean(Y1),
                       beta = Y1.bar - Y0.bar,
                       n = n(),
                       p.Tx = mean(Z))
  site.var = sd(sites$beta)
  scat( "\n\nDescription of site distribution:" )
  scat( "\n\tsite average variation:\tsd(Y0.bar) = %.3f\tsd(Y1.bar) = %.3f",
        sd(sites$Y0.bar), sd(sites$Y1.bar))
  scat( "\n\tmarginal variation:\tsd(Y0) = %.3f\t\tsd(Y1) = %.3f\t\tcor(Y0,Y1) = %.3f",
        sd(data$Y0), sd(data$Y1), cor(data$Y0, data$Y1))
  scat( "\n\tsd( beta ) = %.3f\n\tcorr( Y0.bar, Y1.bar ) = %.3f\tcov = %.3f",
        site.var, cor(sites$Y0.bar, sites$Y1.bar), cov( sites$Y0.bar, sites$Y1.bar))
  scat( "\n\tcor( Y0.bar, beta ) = %.3f", cor(sites$Y0.bar, sites$beta))
  scat( "\n\tmean( n ) = %.3f\tsd( n ) = %.3f", mean(sites$n), sd( sites$n))
  scat( "\n\tmean( Z.bar ) = %.3f\t(%.3f)\n", mean(sites$p.Tx), sd(sites$p.Tx))
  invisible(sites)
}
