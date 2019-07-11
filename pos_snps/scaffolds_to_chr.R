library(dplyr)
library(tidyr)
library(ggplot2)
# this script translates local scaffold positions
# to chromosome positions from Wallberg 2015
# and it calculates recombination rates for each position based on a recomb. map


# I don't trust the orientation but I need it to 
# match to the genetic map and get map positions for all SNPs

# + means oriented forwards
# 0 means not oriented but assumed forwards in their genome/map
# e.g. scaffold 2.1/Group2 start same way
# - means oriented backwards
scaffold_orientation <- read.csv("../data/recomb_map/Wallberg_2015/genome/Amel_4.5.SCAFtoCHROM.AGP",
                                   sep = "\t", header = F, stringsAsFactors = F, na.strings = c("")) %>%
  filter(!V5=="N") # filter out 'gap' scaffolds of 50k N's

colnames(scaffold_orientation) <- c("chr", "start", "end", "n_scaffold", "Ws", 
                                    "scaffold", "ones", "length", "orientation")

# helper function to walk forwards of backwards from the scaffold
# start position, based on whether the orientation is + or -
get_new_pos <- function(scaffold, scaffold_pos){
  if (scaffold_orientation$orientation[scaffold_orientation$scaffold == scaffold] == "-"){
    pos = scaffold_orientation$end[scaffold_orientation$scaffold == scaffold] - (scaffold_pos - 1)
  } else{
    pos = scaffold_orientation$start[scaffold_orientation$scaffold == scaffold] + (scaffold_pos - 1)
  }
  return(pos)
}

# get SNP id, scaffold and scaffold position 
# for each included SNP in the Ramirez SNP set
SNPs <- read.table("../data/bees_new_positions/ALL.map",
                   stringsAsFactors = F, header = F, sep = " ") %>%
  rename(snp_id = V2) %>%
  separate(., snp_id, c("snp1", "chr", "scaffold_name", "scaffold_pos"), remove = F) %>%
  mutate(scaffold = paste(chr, scaffold_name, sep = ".")) %>%
  mutate(scaffold_pos = as.numeric(scaffold_pos)) %>%
  dplyr::select(snp_id, chr, scaffold, scaffold_pos) 

SNPs$pos <- sapply(1:nrow(SNPs), function(i) 
  get_new_pos(scaffold = SNPs[i, "scaffold"], scaffold_pos = SNPs[i, "scaffold_pos"]))

# test cases
#get_new_pos(scaffold = "Group1.1", scaffold_pos = 800)
#get_new_pos(scaffold = SNPs[1,"scaffold"], scaffold_pos = SNPs[1,"scaffold_pos"])
#scaffold_orientation[scaffold_orientation$scaffold=="Group6.1", ]
#get_new_pos(scaffold = "Group6.1", scaffold_pos = 792224) # backwards scaffold

# write out file with new positions
options(scipen=999) # turns off scientific notation
SNPs %>%
  dplyr::select(snp_id, pos) %>%
  write.table(., "results/Ramirez_SNPs_pos_on_Wallberg_Amel4.5_chr.txt",
              quote = F, col.names = F, row.names = F)
options(scipen=0)
map <- read.csv("../data/recomb_map/Wallberg_2015/recombination_rates/A.rates.1000.201.low_penalty.csv.cM_Mb.windows.100000.csv",
                  stringsAsFactors = F, header = F, sep = "\t")
colnames(map) <- c("chr", "pos_start", "map_region_n", "r")

# there's no map estimate for one 100kb region on chr2
# which has two 50kb gaps around a short scaffold: Group2.16 length = 1185
# but there's also no SNPs in my set in this region
# SNPs[SNPs$chr == "Group2" & SNPs$pos > 8600000 & SNPs$pos < 8700000, ]
# so I'll just give the gap a genomewide mean recombination rate
# which won't be used anyways since there are no SNPs and it's a really big 100kb gap
map$r[is.na(map$r)] <- mean(map$r, na.rm = T)
map$rlength <- r*100000
# r has units cM/Mb and a mean ~ 25cM/Mb. plink also uses cM
# I want cM positions relative to the chromosome
# I'll can use approxfun() to do this

# write bed file for recombination rates
options(scipen=999) # turn off scientific notation
map %>%
  mutate(bed_start = pos_start - 1) %>%
  mutate(bed_end = bed_start + 100000) %>%
  dplyr::select("chr", "bed_start", "bed_end", "r") %>%
  write.table("results/rmap_Wallberg2015_100kb.bed",
              quote = F, col.names = F, row.names = F, sep = "\t")
options(scipen=0)

# get smaller 10kb map that doesn't go over gaps between scaffolds
map10 <- read.csv("../data/recomb_map/Wallberg_2015/recombination_rates/A.rates.1000.201.low_penalty.csv.cM_Mb.windows.10000.csv",
                stringsAsFactors = F, header = F, sep = "\t")
colnames(map10) <- c("chr", "pos_start", "map_region_n", "r")

# write befile for map
# first classify category recomb bin 1-5 or 1-10
map10_bed <- map10 %>%
  mutate(bed_start = pos_start - 1) %>%
  mutate(bed_end = bed_start + 10000) %>%
  filter(!is.na(r)) # filter out any 10kb regions (about 5%) without map estimates
  
map10_bed$r_bin5 <- cut(map10_bed$r,  # note need to load map from scaffolds_to_chr.R
      breaks = unique(quantile(map10_bed$r, # I extend bounds so that every value gets in a bin
                               p = seq(0, 1, by = .2))),
      right = T,
      include.lowest = T)
map10_bed$r_bin10 <- cut(map10_bed$r,  # note need to load map from scaffolds_to_chr.R
                        breaks = unique(quantile(map10_bed$r, # I extend bounds so that every value gets in a bin
                                                 p = seq(0, 1, by = .1))),
                        right = T,
                        include.lowest = T)
map10_bed$r_bin5_n <- as.numeric(map10_bed$r_bin5)
map10_bed$r_bin10_n <- as.numeric(map10_bed$r_bin10)
table(map10_bed$r_bin10)
table(map10_bed$r_bin5_n)
options(scipen=999) # turn off scientific notation
map10_bed %>% 
  dplyr::select("chr", "bed_start", "bed_end", "r", "r_bin5_n", "r_bin10_n") %>%
  write.table(., "results/rmap_Wallberg2015_10kb.bed",
              quote = F, col.names = F, row.names = F, sep = "\t")
options(scipen=0)

# write files with the recombination rate values for 5 and 10 recomb bins
write.table(levels(map10_bed$r_bin10), "results/rmap_Wallberg2015_10kb.rbin10.levels",
            col.names = F, row.names = F)
write.table(levels(map10_bed$r_bin5), "results/rmap_Wallberg2015_10kb.rbin5.levels",
            col.names = F, row.names = F)


SNPs_NGSAdmix <- read.table("../global_ancestry/results/input/ordered_scaffolds_CA_AR_MX_harpur_sheppard_kohn_wallberg_NoGroupUn_prunedBy10.snplist",
                   stringsAsFactors = F, header = F, sep = " ") %>%
  rename(snp_id = V1) %>%
  separate(., snp_id, c("scaffold", "scaffold_pos"), sep = "_", remove = F) %>%
  mutate(scaffold_pos = as.numeric(scaffold_pos)) %>%
  dplyr::select(scaffold, scaffold_pos, snp_id) 

SNPs_NGSAdmix$pos <- sapply(1:nrow(SNPs_NGSAdmix), function(i) 
  get_new_pos(scaffold = SNPs[i, "scaffold"], scaffold_pos = SNPs[i, "scaffold_pos"]))
SNPs_NGSAdmix$orig_order <- 1:nrow(SNPs_NGSAdmix)
# write .bed file for NGSAdmix SNP positions
# ! sorts - does NOT have original order matching .beagle.gz GL file
# so it needs to be matched on snp
options(scipen=999) # turns off scientific notation
SNPs_NGSAdmix %>%
  mutate(start = pos - 1, end = pos) %>%
  separate(scaffold, c("chr", "scaffold_n")) %>%
  mutate(chr = factor(chr, levels = unique(.$chr))) %>% # levels keeps original chr order (not alphabetical)
  dplyr::select(chr, start, end, snp_id, orig_order) %>%
  arrange(chr, start) %>%
  write.table(., "../global_ancestry/results/input/ordered_scaffolds_CA_AR_MX_harpur_sheppard_kohn_wallberg_NoGroupUn_prunedBy10.bed",
              sep = "\t", quote = F, col.names = F, row.names = F)
options(scipen=0)

