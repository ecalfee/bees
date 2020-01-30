#!/usr/bin/env Rscript
library(boot)
# This script estimates pi within ancestry, bootstrapped over SNPs
# to run:
# Rscript bootstrap_pi_within.R A 100 100

# arguments
args = commandArgs(trailingOnly=TRUE)
ANCESTRY = args[1]
BOOT_SIZE = as.integer(args[2])
SEED = as.integer(args[3])

if (ANCESTRY == "Combined"){
  load("results/hets_small_sample_combined.RData")
  hets = hets_small_sample
} else{
  ACM = c("A", "C", "M")
  load("results/hets_small_sample_by_ancestry.RData")
  hets = hets_small_sample_by_ancestry[[which(ACM == ANCESTRY)]]
  rm(hets_small_sample_by_ancestry) # get rid of other ancestries from memory for now
}

# function to get mean heterozygosity, ignoring NAs
mean_one_pop <- function(data, indices) {
  d <- data[indices] # allows boot to select sample
  return(mean(d))
}
bootstrap <- function(x, boot_size = BOOT_SIZE){
  x0 = x[!is.na(x)] # bootstrap sample will be equal to the non-NA sample size of SNP allele frequencies
  if (length(x0) == 0){
    return(NULL) # no data, no bootstrap
  }else{
    # bootstrapping replicates
    results <- boot(data = x0, statistic = mean_one_pop, R = boot_size)
    return(results)
  }
}
# test:
#bootstrap(y0)
#bootstrap(c())

# set random seed and run bootstrap:
set.seed(SEED)
boots <- apply(hets, 2, bootstrap)
save(list = "boots", file = paste0("results/bootstrap_pi_within_", ANCESTRY, "_", SEED, "_boots.RData"))
