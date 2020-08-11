#!/usr/bin/env Rscript

# this script takes bp positions on a scaffold from a file
# and converts them to (1) physical bp positions 
# relative to the whole chromosome,
# assuming 50,000 bp gaps between scaffolds
# and (2) recombination positions in Morgans between each position
# where the first position of each chromosome has r = 1 Morgan
# (and will be ignored by ancestry_hmm because it has no preceding SNP)

library(dplyr)
library(tidyr)

# to run: Rscript ./bp_to_r_Wallberg2015.R thin1kb_common3 A.rates.1000.201.low_penalty.csv.cM_Mb.windows.100000.csv
# the result is a position file for input to ancestry_hmm, results/SNPs/thin1kb_common3/included.rmap
# with recombination rates based on a map resolution at 100kb scale

# user input:
args <- commandArgs(trailingOnly=TRUE)
PREFIX_FILE_OUT = paste0("results/SNPs/", args[1], "/included")
VAR_SITES_FILE = paste0(PREFIX_FILE_OUT, ".var.sites")
MAP_FILE = paste0("../data/recomb_map/Wallberg_2015/recombination_rates/", args[2])
SCAFFOLD_STARTS = "../data/honeybee_genome/Amel_4.5_scaffold_start_positions.csv"


# start of new chromosome/Linkage Group default distance (in cM; will be divided by 100 to convert to morgans) 
# this number cannot be negative but is otherwise ignored by ancestry_hmm; 1 Morgan is recommended value
startingBP = 100

# Wallberg map with different resolutions can be chosen -- note: blanks are NAs, no recomb. rate estimated for that window
rmap = read.table(MAP_FILE,
                stringsAsFactors = F, header = F, na.strings = "", sep = "\t")
colnames(rmap) = c("LG", "pos", "window", "cM_Mb")
# This map has blanks (now NAs) where recombination rates cannot be estimated due to lack of SNPs or gaps (e.g. between scaffolds)
# for these, I impute a chromosome-level average recombination rate
# calculate from the average rates across all windows for each Linkage Group
meanLG = tapply(rmap$cM_Mb, rmap$LG, function(r) mean(r, na.rm = T))
for (lg in paste0("Group", 1:16)){
  rmap[rmap$LG==lg & is.na(rmap$cM_Mb), "cM_Mb"] <- meanLG[lg]
} 

# get SNP position file
d0 = read.table(VAR_SITES_FILE, stringsAsFactors = F, header = F, sep = "\t")
colnames(d0) <- c("id", "scaff_pos", "ref_allele", "alt_allele") # id is the scaffold id, e.g. Group1.1
d1 = d0 %>%
separate(., id, c("LG", "scaffold"), sep = "[.]", remove = FALSE) # separate id Group1.7 into LG Group1 and scaffold 7
mode(d1$scaffold) <- "numeric" # change type from character to number

# find position relative to chromosome, not scaffold, where we assume there is a 50,000 bp gap between any adjacent scaffolds
# first get absolute start positions for each scaffold
starts <- read.csv(SCAFFOLD_STARTS, stringsAsFactors = F)
starts$pos_start = starts$pos - 1 # where 0 is the first bp of the scaffold
starts$LG = paste0("Group", starts$lg) # reformat 4 -> Group4 to match d
# now I just add the scaffold-relative position to their start position to get an absolute chromosome position per marker
d = left_join(d1, starts[ , c("LG", "scaffold", "pos_start")], 
              by = c("LG", "scaffold")) %>%
  mutate(chr_pos = scaff_pos + pos_start)
  
# visualize scaffold start positions (jumps = longer scaffold)
#d %>% ggplot(., aes(y = pos_start, x = scaffold, color = LG)) + geom_point()
# visualize chromosome-relative positions together
#d %>% ggplot(., aes(y = chr_pos, x = seq_along(chr_pos), color = LG)) + geom_point()

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
write.table(d[, c("LG", "chr_pos", "rDist")], file = paste0(PREFIX_FILE_OUT, ".rmap"),
            quote = F, row.names = F, col.names = F, sep = "\t")

write.table(d[, c("LG", "chr_pos")], file = paste0(PREFIX_FILE_OUT, ".pos"),
            quote = F, row.names = F, col.names = F, sep = "\t")

write.table(d[, c("rDist")], file = paste0(PREFIX_FILE_OUT, ".rdist"),
            quote = F, row.names = F, col.names = F, sep = "\t")

