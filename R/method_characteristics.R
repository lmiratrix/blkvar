#' Get list of methods in package along with characteristics of those methods
#'
#' Each method name comes with whether it targets average impact for individuals
#' or average impact of sites, and also whether it appears to be a finite-sample
#' method (assuming fixed sites, even if individuals are sampled within site) or
#' superpopulation method (assuming sites/blocks are themselves sampled).
#'
#' @return A tibble of characteristics.
#'
#' @export

 method_characteristics <- function() {
   # Code to make the hard-coded list of characteristics
   if (FALSE) {
     dat <- make_obs_data(n_k = 4:10, p = 0.2)
     a <- compare_methods(data = dat[c("Yobs", "Z", "B")])
     a <- a[1]
     a$fullname <- a$method
     a$method <- c("hybrid_m", "hybrid_p", "plug_in_big", "DB-FP-Persons", "DB-FP-Sites", "DB-SP-Persons", "DB-SP-Sites", "FE",
       "FE-Het", "FE-CR", "FE-IPTW(n)", "FE-IPTW", "FE-IPTW-Sites", "FE-Int-Sites", "FE-Int-Persons", "RICC", "FIRC", "RIRC")
     a$finite <- c(1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0)
     a$site <- c(0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1)
     a$weight <- ifelse(a$site == 1, "site", "person")
     a$population <- ifelse(a$finite, "finite", "superpop")
     a$site <- a$finite <- NULL
     a$biased <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1)
     datapasta::tribble_paste(a)
   }
   tibble::tribble(
        ~fullname,            ~method,  ~weight, ~population, ~biased,
        "hybrid_m",       "hybrid_m", "person",    "finite",       0,
        "hybrid_p",       "hybrid_p", "person",    "finite",       0,
        "plug_in_big",    "plug_in_big", "person",    "finite",       0,
        "DB (individual-finite)",  "DB-FP-Persons", "person",    "finite",       0,
        "DB (site-finite)",    "DB-FP-Sites",   "site",    "finite",       0,
        "DB (individual-superpop)",  "DB-SP-Persons", "person",  "superpop",       0,
        "DB (site-superpop)",    "DB-SP-Sites",   "site",  "superpop",       0,
        "FE",             "FE", "person",    "finite",       1,
        "FE (sand)",         "FE-Het", "person",    "finite",       1,
        "FE (cluster)",          "FE-CR", "person",  "superpop",       1,
        "FE (club)",          "FE-Club", "person",  "superpop",       1,
        "IPTW weighted regression (naive)",     "FE-IPTW(n)", "person",    "finite",       0,
        "IPTW weighted regression",        "FE-IPTW", "person",    "finite",       0,
        "IPTW weighted regression (site)",  "FE-IPTW-Sites",   "site",    "finite",       0,
        "IPTW weighted regression (site, naive)",  "FE-IPTW-Sites(n)",   "site",    "finite",       0,
        "FE interact (site)",   "FE-Int-Sites",   "site",    "finite",       0,
        "FE interact (indiv)", "FE-Int-Persons", "person",    "finite",       0,
        "RICC",           "RICC", "person",    "finite",       1,
        "FIRC",           "FIRC",   "site",  "superpop",       1,
        "RIRC",           "RIRC",   "site",  "superpop",       1
    )
}