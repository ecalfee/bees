# bootstrap steep clines enrichment

#!/usr/bin/env Rscript
# Author: Erin Calfee 2018
# This script calculates n bootstrap statistics
# of % low r snps are within the set of the steepest clines
# across bootstrap resamples of 0.2cM windows genomewide
library(dplyr)
library(tidyr)

# take in command line arguments
args = commandArgs(trailingOnly=TRUE)

# to run: Rscript bootstrap_steep_clines.R test_boot 10 
PREFIX = args[1] # prefix for output file
N = as.numeric(args[2]) # number of bootstraps
dir.create("results/BOOT_r_steep_clines")
file_out = paste0("results/BOOT_r_steep_clines/", PREFIX, ".txt")

# get original data:
load("results/sites_r_windows_steep_slopes.RData") # c("sites_r_windows", "top_1perc", "top_5perc"))

# block bootstrap with genomic data: Bickel, P. J., Boley, N., Brown, J. B., Huang, H., & Zhang, N. R. (2010). Subsampling Methods for Genomic Inference.The Annals ofApplied Statistics,4(4), 1660-1697.http://dx.doi.org/10.1214/10-AOAS363
# this follows ENCODE Strategy..generally though they deal with homogeneous blocks..and it's overly conservative if blocks are not homogeneous...better to segment the blocks in this case.

# set up bootstrap 
# the number of items in the indep. sample should be the number of windows, so I'll resample windows and keep all their points
windows_with_data = sites_r_windows %>% # split data into a list grouped by 0.2cM windows
  split(., sites_r_windows$window_0.2cM)

n_windows = length(windows_with_data)# some windows are excluded (have no data within 0.2cM)

# define bootstrap function
bootstrap_steep_clines <- function(windows_with_data, n_windows, n_boot, top_1perc, top_5perc, name = PREFIX){
  boot = do.call(rbind, # get bootstrap sample
                 sample(x = windows_with_data, size = n_windows, replace = T))
  top_1 = quantile(boot$b, .99)
  top_5 = quantile(boot$b, .95)
  boot %>%
    group_by(., r_bin5) %>%
    summarise(top1_orig = sum(b > top_1perc), # original cutoff
              top5_orig = sum(b > top_5perc), 
              top1_sample = sum(b > top_1), # bootstrap sample cutoff
              top5_sample = sum(b > top_5),
              mean_b = mean(b),
              mean_w = mean(4/b),
              n = n()) %>%
    pivot_wider(data = ., names_from = r_bin5, values_from = c(top1_orig, top5_orig, top1_sample, top5_sample, mean_b, mean_w, n)) %>%
    mutate(name = name, 
           n_boot = n_boot,
           top1_cutoff = top_1,
           top5_cutoff = top_5)
}
# run N bootstraps
boots = do.call(rbind,
                lapply(1:N, function(i) bootstrap_steep_clines(windows_with_data = windows_with_data,
                                           n_windows = n_windows,
                                           n_boot = i,
                                           top_1perc = top_1perc,
                                           top_5perc = top_5perc,
                                           name = PREFIX)))
# write output to file
write.table(boots, file_out,
              quote = F,
              sep = "\t", col.names = T, row.names = F)



