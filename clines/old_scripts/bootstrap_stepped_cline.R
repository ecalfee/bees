#!/usr/bin/env Rscript
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(nlstools)
library(nls.multstart)

# This script estimates a stepped cline, bootstrapping individuals within populations
# to run:
# Rscript bootstrap_stepped_cline.R 10 100 # test with just 10 bootstrap replicates

# arguments
args = commandArgs(trailingOnly=TRUE)
BOOT_SIZE = as.integer(args[1])
SEED = as.integer(args[2])

# cline functions:
source("cline_functions.R") # loads logistic and stepped clines etc.

# cline fit:
fit_stepped_cline_d <- function(data, rescale = 0.84, supp_errors = 'N'){
  cline_stepped <- nls_multstart(alpha_rescaled ~ ifelse(lat > y + dR,
                                                         1 - stepped_cline_edge(x = lat, 
                                                                                y = y, 
                                                                                theta = thetaR, 
                                                                                w = -w, 
                                                                                d = dR),
                                                         ifelse(lat < y - dL,
                                                                stepped_cline_edge(x = lat, 
                                                                                   y = y, 
                                                                                   theta = thetaL, 
                                                                                   w = w, 
                                                                                   d = -dL),
                                                                stepped_cline_center(x = lat, 
                                                                                     y = y, 
                                                                                     w = w))), 
                                 data =  data %>%
                                   mutate(alpha_rescaled = alpha/rescale), 
                                 lower = c(y = -40, dR = 0, thetaR = 0, w = 0, dL = 0, thetaL = 0),
                                 upper = c(y = -25, dR = 5, thetaR = 1, w = 15, dL = 5, thetaL = 1),
                                 start_lower = list(y = -40, dR = 0, thetaR = 0, w = 0, dL = 0, thetaL = 0),
                                 start_upper = list(y = -25, dR = 5, thetaR = 1, w = 15, dL = 5, thetaL = 1), # width here is in degrees latitude
                                 iter = 100,
                                 supp_errors = supp_errors)
  return(cline_stepped)
}

# load data:
load("results/d_A.RData") # load data
d_A_SA <- filter(d_A, continent == "S. America")

# run bootstrap
set.seed(SEED)
boots_cline_stepped_d <- lapply(1:BOOT_SIZE, function(i)
  d_A_SA %>%
    split(.$population) %>%
    purrr::map(~ dplyr::sample_frac(., size = 1, replace = T)) %>%
    bind_rows(.) %>%
    fit_stepped_cline_d(data = ., rescale = 0.84, supp_errors = 'Y'))
save(file = paste0("results/bootstrap_stepped_cline_SA_seed", SEED, ".RData"),
     list = c("boots_cline_stepped_d"))

# save summary of results as txt file
do.call(rbind, lapply(boots_cline_stepped_d, function(x)
  if (is.null(x)){ # no model fit
    return(NULL)
  } else{ # model fit
    bind_cols(broom::glance(x), 
              broom::tidy(x) %>% 
              dplyr::select(c("term", "estimate")) %>% 
              tidyr::pivot_wider(data = ., names_from = "term", values_from = "estimate"))

  })) %>%
write.table(., file = paste0("results/bootstrap_stepped_cline_SA_seed", SEED, ".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

