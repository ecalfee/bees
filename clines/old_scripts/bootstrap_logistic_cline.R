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
fit_logistic_cline <- function(data, supp_errors = 'N'){
  cline <- nls_multstart(alpha ~ logistic4(x = abs_lat,
                                           b = b + b_SA*S_America,
                                           mu = mu + mu_SA*S_America,
                                           K = 0.84),
                         start_lower = list(b = -1, mu = 25, b_SA = -1, mu_SA = -5),
                         start_upper = list(b = 1, mu = 40, b_SA = 1, mu_SA = 5),
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
boots_cline <- lapply(1:BOOT_SIZE, function(i)
  d_A %>%
    split(.$population) %>%
    purrr::map(~ dplyr::sample_frac(., size = 1, replace = T)) %>%
    bind_rows(.) %>%
    fit_logistic_cline(data = ., supp_errors = 'Y'))
save(file = paste0("results/bootstrap_logistic_cline_seed", SEED, ".RData"),
     list = c("boots_cline"))

# save summary of results as txt file
do.call(rbind, lapply(boots_cline, function(x)
  if (is.null(x)){ # no model fit
    return(NULL)
  } else{ # model fit
    bind_cols(broom::glance(x),
              broom::tidy(x) %>%
              dplyr::select(c("term", "estimate")) %>%
              tidyr::pivot_wider(data = ., names_from = "term", values_from = "estimate"))

  })) %>%
write.table(., file = paste0("results/bootstrap_logistic_cline_seed", SEED, ".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)


# # to visualize:
# # bootstraps for uncertainty in cline parameter estimates:
# # (individuals re-sampled within populations)
# # bootstraps:
# boots_logistic <- read.table("results/bootstrap_logistic_cline_seed100.txt",
#                              sep = "\t", header = T, stringsAsFactors = F)
# table(boots_logistic$isConv) # all 1000 converged -- NOTE re-doing for 10000 bootstraps
# # pivot=basic bootstrap CI's= estimate - quantile(boots - estimate, c(0.975, 0.025))
# # NA center
# CI_NA_mu = coef_lat_zone_mu_b_.84$mu - quantile(boots_logistic$mu - coef_lat_zone_mu_b_.84$mu, c(0.975, 0.025))
# # SA center
# CI_SA_mu = (coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA) - quantile(boots_logistic$mu + boots_logistic$mu_SA - (coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA), c(0.975, 0.025))
# # diff center
# CI_diff_mu = coef_lat_zone_mu_b_.84$mu_SA - quantile(boots_logistic$mu_SA - coef_lat_zone_mu_b_.84$mu_SA, c(0.975, 0.025))
# # NA slope
# CI_NA_b = coef_lat_zone_mu_b_.84$b - quantile(boots_logistic$b - coef_lat_zone_mu_b_.84$b, c(0.975, 0.025))
# # SA slope
# CI_SA_b = (coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA) - quantile(boots_logistic$b + boots_logistic$b_SA - (coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA), c(0.975, 0.025))
# # diff slope
# CI_diff_b = coef_lat_zone_mu_b_.84$b_SA - quantile(boots_logistic$b_SA - coef_lat_zone_mu_b_.84$b_SA, c(0.975, 0.025))
# 
# coef_lat_zone_mu_b_.84
# CI_SA_mu
# CI_NA_mu
# CI_diff_mu
# CI_SA_b
# CI_NA_b
# CI_diff_b
# 
# 4/CI_NA_b
# 4/CI_SA_b
# 4/CI_diff_b # not right
