#!/usr/bin/env Rscript

# this script converts var sites used in hmm ancestry inference
# into a bed file of tracts around those sites

# to run:

# Rscript ./ancestry_sites_to_tracts.R results/SNPs/thin1kb_common3/included


# arguments
args = commandArgs(trailingOnly=TRUE)
# directory with chri.var.sites files
PATH = args[1]

# helper function makes start and end positions of tracts around a SNP site
# as the bp position halfway between that marker and the next closest marker.
# or, if marker has no nearest neighbor (at either end of chromosome),
# tract ends at the marker on one side
get_start_end <- function(path){
  sites0 <- read.table(paste0(path, ".var.sites"),
                      stringsAsFactors = F, header = F)
  snp_list <- read.table(paste0(path, ".snplist"),
  stringsAsFactors = F, header = F)
  colnames(sites0) <- c("chrom", "pos", "allele1", "allele2")
  colnames(snp_list) <- c("snp_id")
  sites <- cbind(sites0, snp_list)
  sites$start <- floor(sites$pos - diff(c(sites$pos[1]-1, sites$pos), lag = 1)/2)
  sites$end <- floor(sites$pos + diff(c(sites$pos, sites$pos[length(sites$pos)]), lag = 1)/2)
  sites$isStartChr <- c("hi", sites$chrom[-length(sites$chrom)]) != sites$chrom
  sites$isEndChr <- sites$chrom != c(sites$chrom[-1], "there")
  # special cases at start and end of the chromosome/scaffold 
  sites$start[sites$isStartChr] <- sites$pos[sites$isStartChr] - 1
  sites$end[sites$isEndChr] <- sites$pos[sites$isEndChr]
  return(sites[, c("chrom", "start", "end", "snp_id")])
}
# find start and end positions for all sites across 10 chromosomes
all_sites <- get_start_end(path = PATH)

# write sites to output file
write.table(all_sites, paste0(PATH, ".var.sites.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)
