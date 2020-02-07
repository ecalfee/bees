#!/usr/bin/env Rscript
library(dplyr)
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
  rm(hets_small_sample)
} else{
  ACM = c("A", "C", "M")
  load("results/hets_small_sample_by_ancestry.RData")
  hets = hets_small_sample_by_ancestry[[which(ACM == ANCESTRY)]]
  rm(hets_small_sample_by_ancestry) # get rid of other ancestries from memory for now
}

# get recombination position for all SNP sites:
snp_sites <- do.call(rbind, 
                      lapply(1:16, function(i) read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group", i, 
                                                                 ".rpos"),
                                                          header = F, stringsAsFactors = F))) %>%
  data.table::setnames(c("scaffold", "pos", "major", "minor", "pos_cM"))
# divide into 1 cm bins:
snp_sites$bin_1cM <- paste(snp_sites$scaffold, floor(snp_sites$pos_cM), sep = "-")

# function to bin heterozygosity into 1cM bins
hets_by_bin <- function(pop) {
  snp_sites %>%
  mutate(het = hets[ , pop]) %>%
  group_by(scaffold, bin_1cM) %>%
  summarise(n = sum(!is.na(het)),
            het_per_snp = mean(het, na.rm = T))}
# at end, scale heterozygosity per snp by snps per bp

# function to get mean heterozygosity, ignoring NAs
one_boot <- function(d){
  sum(d$n*d$het_per_snp)/sum(d$n) # mean het per snp over all snps
}
# bootstrap function
bootstrap <- function(x, boot_size = BOOT_SIZE){
  x0 = x[!is.na(x$het_per_snp), ] # single bootstrap sample size will be equal to the non-NA sample size of SNP allele frequencies
  if (length(x0) == 0){
    return(NULL) # no data, no bootstrap
  }else{
    # original estimate
    estimate <- one_boot(d = x0)
    
    # number of bins and chromosomes with data
    n_bins <- nrow(x0)
    n_chr <- length(unique(x$scaffold))
    
    # bootstrapping replicates - resample all rows to same size as original data
    boots <- sapply(1:boot_size, function(i) one_boot(d = dplyr::sample_frac(x0, 1, replace = TRUE)))
    ci = estimate - quantile(boots - estimate, c(0.975, 0.025))
    return(data.frame(estimate = estimate, lower = ci[1], upper = ci[2], 
                      chr = n_chr, n_bins = n_bins, stringsAsFactors = F))
  }
}

# set random seed and run bootstrap:
set.seed(SEED)
POPS = colnames(hets)
boots = do.call(rbind, lapply(POPS, function(p) bootstrap(x = hets_by_bin(p), boot_size = BOOT_SIZE) %>%
                        mutate(pop = p)))
# save results
save(list = "boots", file = paste0("results/bootstrap_pi_within_", ANCESTRY, "_", SEED, "_boots.RData"))
