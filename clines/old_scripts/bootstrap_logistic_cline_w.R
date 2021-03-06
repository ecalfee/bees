#!/usr/bin/env Rscript
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(nlstools)
library(nls.multstart)

# This script estimates logistic clines for NA and SA,
# bootstrapping individuals within populations
# to run:
# Rscript bootstrap_logistic_cline.R 10 100 # test with just 10 bootstrap replicates

# arguments
args = commandArgs(trailingOnly=TRUE)
BOOT_SIZE = as.integer(args[1])
SEED = as.integer(args[2])

# cline functions:
source("cline_functions.R") # loads logistic and stepped clines etc.

# cline fit:
fit_logistic_cline_w <- function(data, supp_errors = 'N'){
  cline <- nls_multstart(alpha ~ logistic4(x = abs_lat,
                                           b = -4/(w + w_SA*S_America),
                                           mu = mu + mu_SA*S_America,
                                           K = 0.84),
                         start_lower = list(w = 0, mu = 25, b_SA = -1, mu_SA = -5),
                         start_upper = list(w = 15, mu = 40, b_SA = 1, mu_SA = 5),
                         supp_errors = supp_errors,
                         iter = 250,
                         convergence_count = 100,
                         data = data)
  return(cline)
}

# load data:
load("results/d_A.RData") # load data

# run bootstrap
set.seed(SEED)
boots_cline_w <- lapply(1:BOOT_SIZE, function(i)
  d_A %>%
    split(.$population) %>%
    purrr::map(~ dplyr::sample_frac(., size = 1, replace = T)) %>%
    bind_rows(.) %>%
    fit_logistic_cline_w(data = ., supp_errors = 'Y'))
save(file = paste0("results/bootstrap_logistic_cline_w_seed", SEED, ".RData"),
     list = c("boots_cline_w"))

# save summary of results as txt file
do.call(rbind, lapply(boots_cline_w, function(x)
  if (is.null(x)){ # no model fit
    return(NULL)
  } else{ # model fit
    bind_cols(broom::glance(x),
              broom::tidy(x) %>%
              dplyr::select(c("term", "estimate")) %>%
              tidyr::pivot_wider(data = ., names_from = "term", values_from = "estimate"))

  })) %>%
write.table(., file = paste0("results/bootstrap_logistic_cline_w_seed", SEED, ".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

            # # change parameterization to estimate w directly with bootstrap
            # 
            # 
            # CI_NA_w = 4/coef_lat_zone_mu_b_.84$b - quantile(4/boots_logistic$b - 4/coef_lat_zone_mu_b_.84$b, c(0.975, 0.025))
            # CI_NA_w
            # # SA slope
            # CI_SA_w = (4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA)) - quantile(4/(boots_logistic$b + boots_logistic$b_SA) - 4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA), c(0.975, 0.025))
            # CI_SA_w
            # CI_diff_w2 = (4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA) - 4/(coef_lat_zone_mu_b_.84$b)) -
            #                 quantile(4/(boots_logistic$b + boots_logistic$b_SA) - 4/boots_logistic$b -
            #                            (4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA) - 4/coef_lat_zone_mu_b_.84$b), c(0.975, 0.025))
            # CI_diff_w2
            # curve(-4/(-0.5 + x), from = -3, to = 0)
            # CI_SA_w-CI_NA_w
