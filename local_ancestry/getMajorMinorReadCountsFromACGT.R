#!/usr/bin/env Rscript

# input:
# Bee ID, and directory where the countsACGT subdirectory is located

# output (in subdirectory countsMajMin): 
# a counts.txt file with two columns: read counts for major allele,
# and read counts for minor allele for all positions
# (even for positions where the individual has zero coverage, ie. 0 0 counts)


library(dplyr)
# to run:
# Rscript count_reads_major_minor.R AR0302 thin1kb_common3

print(getwd()) # print current directory

# arguments
args = commandArgs(trailingOnly=TRUE)
# Bee ID
ID = args[1]
# DIR name for set of SNPs
DIR = args[2]
# full path to directory
path = paste0("results/SNPs/", DIR)

acgt_file = paste0(path, "/countsACGT/", ID, ".counts.gz")
pos_file = paste0(path, "/countsACGT/", ID, ".pos.gz")
sites_file = paste0(path, "/included.var.sites")
output_directory = paste0(path, "/countsMajMin")
output_file = paste0(output_directory, "/", ID, ".counts.txt")

# input data
SNPs = read.table(sites_file,
                  header = F, stringsAsFactors = F, sep = "\t")
colnames(SNPs) = c("chr", "pos", "major", "minor")
acgt = read.table(acgt_file, stringsAsFactors = F, header = T)
colnames(acgt) = substr(colnames(acgt), 4, 4) # make totA -> A column names
pos = read.table(pos_file, stringsAsFactors = F, header = T)
d = cbind(pos, acgt) %>%
  left_join(SNPs, ., by = c("chr", "pos")) %>%
  tidyr::gather(., "allele", "n", c("A", "C", "G", "T"))
maj_counts = filter(d, major==allele) %>%
  select(-allele) # get read counts for alleles matching major allele
min_counts = filter(d, minor==allele) %>%
  select(-allele) # get read counts for allele matching minor allele
counts = full_join(maj_counts, min_counts,
                   suffix = c("_major", "_minor"),
                   by = c("chr", "pos", "major", "minor", "totDepth")) %>%
left_join(SNPs, ., by = c("chr", "pos", "major", "minor")) # match order of original SNPs var.sites file

# make output directory if it doesn't exist
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = T) # will simply warn if directory already exists

# write output file
write.table(select(counts, c(n_major, n_minor)), 
            output_file,
            na = "0", # write NA's (no coverage) as zero counts
            sep = "\t", row.names = F, col.names = F, quote = F) # just counts, no headers
warnings() # print any warnings
