#!/usr/bin/env Rscript

# this script takes in a Wallberg recombination map by windows,
# a chromosome number,
# and a multimix main directory within which there's a legends file with physical bp positions 
# relative to the whole chromosome,
# (assuming 50,000 bp gaps between scaffolds)
# and returns a multimix map file with 
# (1) recombination rate in cM at each SNP position
# and (2) in cM position of each SNP relative to the chromosome/LG


# to run: Rscript bp_LG_to_r_Wallberg2015_4multimix.R data/sims_downsample/multimix 16 data/recomb_map/Wallberg_2015/recombination_rates/A.rates.1000.201.low_penalty.csv.cM_Mb.windows.100000.csv
# the result is a map file for multimix, e.g. data/multimix/input/haplotypes/genetic_maps/chr*.map
# sample output file:
# position COMBINED_rate(cM/Mb) Genetic_Map(cM)
# 45413 -1 0
# 711153 2.6858076690 0
# 713682 2.8222713027 0.0067924076
# 713754 2.9813105581 0.0069956111
# 718105 2.9806151254 0.0199672934

# user input:
#args <- c("data/sims_downsample/multimix", "1", "data/recomb_map/Wallberg_2015/recombination_rates/A.rates.1000.201.low_penalty.csv.cM_Mb.windows.100000.csv")
args <- commandArgs(trailingOnly=TRUE)
MULTIMIX_DIR = args[1]
CHR = args[2]
RMAP = args[3]
# find 'legends' file in multimix directory
POS_FILE = paste0(MULTIMIX_DIR, "/input/haplotypes/legend_files/chr", CHR, ".legend")

# create output directory
if (!dir.exists(paste0(MULTIMIX_DIR, "/input/haplotypes/genetic_maps"))){
  dir.create(paste0(MULTIMIX_DIR, "/input/haplotypes/genetic_maps"), recursive=TRUE)
}

# Wallberg map with resolution 100kb scale -- note: blanks are NAs, no recomb. rate estimated for that window
rmap = read.csv(RMAP,
                stringsAsFactors = F, header = F, na.strings = "", sep = "\t")
colnames(rmap) = c("LG", "pos", "window", "cM_Mb")
class(rmap$pos) = "integer"
class(rmap$cM_Mb) = "numeric"
# This map has blanks (now NAs) where recombination rates cannot be estimated due to lack of SNPs or gaps (e.g. between scaffolds)
# for these, I impute a chromosome-level average recombination rate
# calculate from the average rates across all windows for each Linkage Group
meanLG = tapply(rmap$cM_Mb, rmap$LG, function(r) mean(r, na.rm = T))
for (lg in paste0("Group", 1:16)){
  rmap[rmap$LG==lg & is.na(rmap$cM_Mb), "cM_Mb"] <- meanLG[lg]
} 
# subset to only linkage group for relevant chromosome
rmapLG = rmap[rmap$LG == paste0("Group", CHR), ]



# get position file
d0 = read.table(POS_FILE, stringsAsFactors = F, header = F)
colnames(d0) <- c("SNP_NAME", "position", "A1", "A2")
class(d0$position) = "integer"

# create helper function to calculate recombination distance between 2 positions
# and recombination rate of 2nd position
rBetween = function(pos1, pos2, rmapLG){
    if (pos1 > pos2){# error; position 1 is greater than position 2
      stop(paste("out of order --", pos1, "seen before", pos2))
    }else{
      # index of start of window for position 1 and 2
      i1 = max(which(rmapLG$pos <= pos1))
      i2 = max(which(rmapLG$pos <= pos2))
      rate2 = rmapLG[i2, "cM_Mb"]
      if (i1==i2){
        # within same recombination window, calculate cM distance between
        dist = (pos2 - pos1) * rmapLG[i1, "cM_Mb"] / 1000000 # division to convert bp to Mbp
      } else{ # otherwise, solve by recursion
        # add distance to end of pos 1's window to remaining distance between positions
        dist = (rmapLG[i1 + 1, "pos"] - pos1) * rmapLG[i1, "cM_Mb"] / 1000000 + 
          rBetween(pos1 = rmapLG[i1 + 1, "pos"], pos2 = pos2, 
                             rmap = rmap)[1]
      }
    }
return(c(dist, rate2))
}
# test cases - all good
#rBetween(pos1 = 1, pos2 = 1, rmapLG = rmap) #0
#rBetween(pos1 = 100, pos2 = 99, rmap) # error
#rBetween(pos1 = 1, pos2 = 10001, rmap) # .418
#rBetween(pos1 = 1, pos2 = 100001, rmap) # 4.18
#rBetween(pos1 = 1, pos2 = 350001, rmap = rmap) # 13.16 = 4.18+4.24+3.46+2.56/2

# apply function to calculate recombination distance
d1 <- rbind(data.frame(SNP_NAME="FIRST_SNP", position=1, A1="Z", A2="Z", stringsAsFactors = F), d0) # add a first SNP at position zero
rBtwn = do.call(rbind, lapply(2:nrow(d1), function(i) 
  rBetween(pos1 = d1[i-1, "position"], pos2 = d1[i, "position"],
           rmapLG = rmapLG)))
d = data.frame("position" = d0[ , "position"], "Genetic_Map" = cumsum(rBtwn[ , 1]), 
               "COMBINED_rate" = c(rBtwn[ , 2])) # ignores rate for 1st SNP (because I added FIRST_SNP at position 1)

# write final output SNP positions file for use by ancestry_hmm
write.table(d[, c("position", "COMBINED_rate", "Genetic_Map")], 
            file = paste0(MULTIMIX_DIR, "/input/haplotypes/genetic_maps/chr", CHR, ".map"),
            # must input names by hand because they have special characters ( that r will otherwise convert "(" -> "."))
            quote = F, row.names = F, col.names = c("position", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)"))


