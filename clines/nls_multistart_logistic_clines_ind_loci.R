# script uses nonlinear least squares to calculate
# best fitting logistic cline for individual loci
# using population A ancestry frequencies in AR

#!/usr/bin/env Rscript
library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)

# to run first 100 snps from observed A ancestry data, creates output results/ind_snp_clines/A_1_100.txt:
# Rscript nls_multistart_logistic_clines_ind_loci.R A 1 100 A_1_100
# args = c("A", "1", "100", "A_1_10")
# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
DATA_NAME = args[1]
INDEX_START = as.integer(args[2])
INDEX_END = as.integer(args[3])
NAME_OUT = args[4]

# output to clines/results
# make sure output directory exists first (or create)!
dir_output = paste0("results/ind_snp_nls_multistart_clines")
if (!file.exists(dir_output)){ 
  dir.create(file.path(dir_output), recursive = T)
}

# load population metadata
load("../local_ancestry/results/pops_by_lat.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 

# load mean ancestry by population
# .RData file contains objects 'sites' which has columns chr and pos and 'DATA_NAME' which has ancestry where sites are rows and pops are columns
load(paste0("../local_ancestry/results/", DATA_NAME, ".RData")) 

# subset data to pops of interest, only SA
A_snp_all <- get(DATA_NAME)[ , meta.AR.order.by.lat$population]

# complain if index is out of range
if (INDEX_END > nrow(A_snp_all)) print(paste("! Warning end index", INDEX_END, 
                                             "exceeds rows in file:", nrow(A_snp_all), 
                                             ". This script will do calculations to the end of the rows"))
INDEX_END <- min(INDEX_END, nrow(A_snp_all))

# define logistic curve
logistic3 <- function(x, mu, b){
  1/(1 + exp(-b*(x - mu)))
}

fits <- data.frame(snp_index = INDEX_START:INDEX_END) %>%
  mutate(fit = lapply(snp_index, 
                      function(i) nls_multstart(A ~ logistic3(x = lat, 
                                                              b = b, mu = mu),
                                                start_lower = list(b = 0, mu = -40),
                                                start_upper = list(b = 1, mu = -25),
                                                supp_errors = 'Y',
                                                iter = 250,
                                                convergence_count = 100,
                                                data = data.frame(A = unname(t(A_snp_all[i, ])), 
                                                                  lat = meta.AR.order.by.lat$lat, 
                                                                  stringsAsFactors = F))))


# get model fit info
info <- fits %>%
  mutate(summary = purrr::map(fit, glance)) %>%
  unnest(summary) %>%
  dplyr::select(-fit) # don't leave in the model fit itself

# write to file
write.table(info, 
            paste0(dir_output, "/", NAME_OUT, ".info"),
            col.names = T, row.names = F, quote = F, sep = "\t")

# get model fit parameters
params <- fits %>%
  mutate(., p = purrr::map(fit, tidy)) %>%
  unnest(p) %>%
  dplyr::select(., c("snp_index", "term", "estimate", "std.error", "statistic", "p.value"))

# calculate confidence intervals
CI <- fits %>%
  mutate(., cis = purrr::map(fit, confint2),
         cis = purrr::map(cis, data.frame)) %>%
  unnest(cis) %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(., snp_index) %>%
  mutate(., term = c('mu', 'b')) %>%
  ungroup() %>%
  dplyr::select(., c("snp_index", "conf.low", "conf.high", "term"))

# merge parameters and CI estimates
params_CI <- left_join(params, CI, by = intersect(names(params), names(CI)))

# write to file parameters with confidence intervals
write.table(params_CI, 
           paste0(dir_output, "/", NAME_OUT, ".params"),
            col.names = T, row.names = F, quote = F, sep = "\t")
