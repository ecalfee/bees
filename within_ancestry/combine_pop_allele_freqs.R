#!/usr/bin/env Rscript
# this script loads population allele frequencies 
# for a pop, and combines them with the sites file,
# producing a new file with a frequency for every site (or NA if no data)
# and number of individuals with data for each site

library(dplyr)

# to run:
# Rscript combine_pop_allele_freqs.R AR01 combined_sept19 A Group1

# arguments
args = commandArgs(trailingOnly=TRUE)

POP = args[1]
PREFIX = args[2]
ANCESTRY = args[3]
SITES_PREFIX = args[4]


# paths to input and output directories
sites_file = paste0("../geno_lik_and_SNPs/results/", PREFIX, "/variant_sites/", SITES_PREFIX, ".var.sites")
path_allele_freqs = paste0("results/", PREFIX, "/", ANCESTRY, "/allele_freq/", SITES_PREFIX)


# read in sites
sites0 <- read.table(sites_file, 
                    stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("scaffold", "pos", "major", "minor"))

allele_freq <- read.table(paste0(path_allele_freqs, "/", POP, ".mafs.gz"),
                          stringsAsFactors = F, header = T,
                          # means will join on empty data.frame
                          colClasses = c("character", "numeric", 
                                         "character", "character",
                                         "character", "numeric",
                                         "integer")) %>%
  left_join(sites0, ., by = c("scaffold"="chromo", "pos"="position",
                              "major", "minor"))

# separate out allele freqs and write to file
f <- allele_freq %>%
  mutate(p = ifelse(phat > 1, 1, phat)) %>% # gets rid of weird angsd rounding error
  dplyr::select(p) %>%
  data.table::setnames(POP)
write.table(f, paste0(path_allele_freqs, "/", POP, ".freqs.txt"),
            col.names = T, row.names = F, quote = F)

# separate out number of individuals with data and write to file
n <- allele_freq %>%
  mutate(n = ifelse(is.na(nInd), 0, nInd)) %>% # NA means no individuals with data
  dplyr::select(n) %>%
  data.table::setnames(POP)
write.table(n, paste0(path_allele_freqs, "/", POP, ".nInd"),
            col.names = T, row.names = F, quote = F)