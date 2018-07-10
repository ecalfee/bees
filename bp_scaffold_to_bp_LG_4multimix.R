#!/usr/bin/env Rscript

# this script takes bp positions on a scaffold from a file
# and converts them to (1) physical bp positions 
# relative to the whole chromosome,
# assuming 50,000 bp gaps between scaffolds
# (2) local recombination rate at that position
# and (3) recombination positions in cM relative to the beginning of the chromosome

# note: will produce error if more than one chromosome is included in input file


library(dplyr)
#library(ggplot2)

# to run: Rscript bp_to_r_Wallberg2015_4multimix.R data/bees_new_positions/lowLDA4/A.frq.counts
# the result is a set of map files for input to MULTIMIX
# OUTPUT is split by chromosome and saved in OUT_DIR/chr*.map files where * is the chromosome

# user input:
args <- commandArgs(trailingOnly=TRUE)
POS_FILE = args[1]
OUT_DIR = args[2]

# start of new chromosome/Linkage Group default distance (in cM; will be divided by 100 to convert to morgans) 
# this number cannot be negative but is otherwise ignored by ancestry_hmm; 1 Morgan is recommended value
startingBP = 100

# Wallberg map with resolution 100kb scale -- note: blanks are NAs, no recomb. rate estimated for that window
rmap = read.csv("data/recomb_map/Wallberg_2015/recombination_rates/A.rates.1000.201.low_penalty.csv.cM_Mb.windows.100000.csv",
                stringsAsFactors = F, header = F, na.strings = "", sep = "\t")
colnames(rmap) = c("LG", "pos", "window", "cM_Mb")
# This map has blanks (now NAs) where recombination rates cannot be estimated due to lack of SNPs or gaps (e.g. between scaffolds)
# for these, I impute a chromosome-level average recombination rate
# calculate from the average rates across all windows for each Linkage Group
meanLG = tapply(rmap$cM_Mb, rmap$LG, function(r) mean(r, na.rm = T))
for (lg in paste0("Group", 1:16)){
  rmap[rmap$LG==lg & is.na(rmap$cM_Mb), "cM_Mb"] <- meanLG[lg]
} 

# get position file
d0 = read.table(POS_FILE, stringsAsFactors = F, header = T)
# reformat to access SNP position information
SNPs = data.frame(t(sapply(d0$SNP, 
                             function(i) 
                               unlist(strsplit(i, split = "[.]")))), stringsAsFactors = F)
colnames(SNPs) = c("snpN", "LG", "scaffold", "scaff_pos")
# script only accepts one chromosome at a time
if (length(unique(SNPs$LG)) > 1){
  stop("error: please subset input data file to only include one LG/chromosome")
}

# convert characters to integers
class(SNPs$scaffold) = "integer"
class(SNPs$scaff_pos) = "integer"
rownames(SNPs) = NULL 
d1 = cbind(d0, SNPs)
# find position relative to chromosome, not scaffold, where we assume there is a 50,000 bp gap between any adjacent scaffolds
# first get absolute start positions for each scaffold
starts <- read.csv("data/honeybee_genome/Amel_4.5_scaffold_start_positions.csv", stringsAsFactors = F)
starts$pos_start = starts$pos - 1 # start of scaffold
starts$LG = paste0("Group", starts$lg) # reformat 4 -> Group4 to match d
# now I just add the scaffold-relative position to their start position to get an absolute chromosome position per marker
d = left_join(d1, starts[ , c("LG", "scaffold", "pos_start")], 
              by = c("LG", "scaffold")) %>%
  mutate(chr_pos = scaff_pos + pos_start)

# create helper function to calculate recombination distance between 2 positions or return 1 if across chromosomes
rBetween = function(pos1, pos2, lg1, lg2, rmap){
  if (lg1 != lg2){
    dist = startingBP # start of new chromosome/Linkage Group 
  }else{ # within-chromosome
    if (pos1 > pos2){# error; position 1 is greater than position 2
      stop(paste("out of order --", lg1, pos1, "seen before", lg2, pos2))
    }else{
      # only consider recomb. map for current linkage group
      rmapLG = rmap[rmap$LG==lg1, ]
      # index of start of window for position 1 and 2
      i1 = max(which(rmapLG$pos <= pos1))
      i2 = max(which(rmapLG$pos <= pos2))
      if (i1==i2){
        # within same recombination window, calculate cM distance between
        dist = (pos2 - pos1) * rmapLG[i1, "cM_Mb"] / 1000000 # division to convert bp to Mbp
      } else{ # otherwise, solve by recursion
        # add distance to end of pos 1's window to remaining distance between positions
        dist = (rmapLG[i1 + 1, "pos"] - pos1) * rmapLG[i1, "cM_Mb"] / 1000000 + 
          rBetween(pos1 = rmapLG[i1 + 1, "pos"], pos2 = pos2, 
                             lg1 = lg1, lg2 = lg2, rmap = rmap)
      }
    }
  }
return(dist)
}
# test cases - all good
#rBetween(pos1 = 1, pos2 = 1, lg1 = "Group1", lg2 = "Group2", rmap = rmap) == -9
#rBetween(pos1 = 1, pos2 = 1, lg1 = "Group1", lg2 = "Group1", rmap = rmap) == 0
#rBetween(pos1 = 100, pos2 = 99, lg1 = "Group1", lg2 = "Group1", rmap = rmap) # error
#rBetween(pos1 = 1, pos2 = 10001, lg1 = "Group1", lg2 = "Group1", rmap = rmap) # .418
#rBetween(pos1 = 1, pos2 = 100001, lg1 = "Group1", lg2 = "Group1", rmap = rmap) # 4.18
#rBetween(pos1 = 1, pos2 = 350001, lg1 = "Group1", lg2 = "Group1", rmap = rmap) # 13.16 = 4.18+4.24+3.46+2.56/2

# apply function to calculate recombination distance
# /100 so distance is converted from cM to Morgans as specified by ancestry_hmm
d$rDist = c(startingBP, sapply(2:length(d$chr_pos), function(i) 
  rBetween(pos1 = d[i-1, "chr_pos"], pos2 = d[i, "chr_pos"],
           lg1 = d[i-1, "LG"], lg2 = d[i, "LG"], rmap = rmap)))/100

# write final output SNP positions file for use by ancestry_hmm
write.table(d[, c("LG", "chr_pos", "rDist")], file = paste0(POS_FILE, ".rmap"),
            quote = F, row.names = F, col.names = F)


