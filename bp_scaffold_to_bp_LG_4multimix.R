#!/usr/bin/env Rscript

# this script takes bp positions on a scaffold from a file
# and converts them to (1) physical bp positions 
# relative to the whole chromosome,
# assuming 50,000 bp gaps between scaffolds

# note: run one chromosome at a time (split bim file in plink by chromosome first)

library(dplyr)

# to run: Rscript bp_scaffold_to_bp_LG_4multimix.R data/bees_new_positions/ALL_splitbychr/ALL_chr1 1 data/sims_downsample/multimix
# the result is a set of map files for input to MULTIMIX
# OUTPUT is saved in $out/input/haplotypes/legend_files/chr*.legend files where * is the chromosome

# user input:
args <- commandArgs(trailingOnly=TRUE)
BIM_FILE = args[1] # looks for a prefix.bim file made by plink
CHR = args[2]
OUT_DIR = args[3]

# create output directory
if (!dir.exists(paste0(OUT_DIR, "/input/haplotypes/legend_files"))){
  dir.create(paste0(OUT_DIR, "/input/haplotypes/legend_files"), recursive=TRUE)
}

# get relative positions (bim file)
d0 = read.table(paste0(BIM_FILE, ".bim"), stringsAsFactors = F, header = F)
colnames(d0) <- c("CHR", "SNP_NAME", "R_POS_FALSE", "BP_POS_FALSE", "A1", "A2")
if (d0$CHR[1] != CHR){
  stop("error: chromosome in bim file doesn't match chromosome in input argument; use plink to subset bim file to chromosome of interest")
}
# reformat to access SNP position information
SNPs = data.frame(t(sapply(d0$SNP_NAME, 
                             function(i) 
                               unlist(strsplit(i, split = "[.]")))), stringsAsFactors = F)
colnames(SNPs) = c("snpN", "LG", "scaffold", "scaff_pos")
rownames(SNPs) = NULL 
# script only accepts one chromosome at a time
if (length(unique(SNPs$LG)) > 1){
  stop("error: please subset input data file in plink first to only include one LG/chromosome")
}

# convert characters to integers
class(SNPs$scaffold) = "integer"
class(SNPs$scaff_pos) = "integer"

d1 = cbind(d0[ , c("SNP_NAME", "A1", "A2")], SNPs)
# find position relative to chromosome, not scaffold, where we assume there is a 50,000 bp gap between any adjacent scaffolds
# first get absolute start positions for each scaffold
starts <- read.csv("data/honeybee_genome/Amel_4.5_scaffold_start_positions.csv", stringsAsFactors = F)
starts$pos_start = starts$pos - 1 # start of scaffold
starts$LG = paste0("Group", starts$lg) # reformat 4 -> Group4 to match d
# now I just add the scaffold-relative position to their start position to get an absolute chromosome position per marker
d = left_join(d1, starts[ , c("LG", "scaffold", "pos_start")], 
              by = c("LG", "scaffold")) %>%
  mutate(chr_pos = scaff_pos + pos_start)

# write final output SNP positions file for use by multimix
write.table(d[, c("SNP_NAME", "chr_pos", "A1", "A2")], file = paste0(OUT_DIR, "/input/haplotypes/legend_files/chr", CHR, ".legend"),
            quote = F, row.names = F, col.names = F, sep = " ")


