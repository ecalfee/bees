#!/usr/bin/env Rscript

# this script takes bp positions and converts to cM position
# relative to the first bp=1 of the scaffold
# for the map, it uses 10kb resolution recombination rates
# and then linearly interpolates the position of each SNP

library(dplyr)
library(tidyr)

MAP_FILE = "../data/recomb_map/Wallberg_HAv3.1/map_rates_extended_10kb.bed"

# Wallberg map based on LD
rmap = read.table(MAP_FILE,
                stringsAsFactors = F, header = T, sep = "\t")

# split into diff chromosomes
chromosomes = paste0("Group", 1:16)
rmap_chr0 = lapply(chromosomes, function(x) filter(rmap, chr == x))
# make positions for map coordinates based on rates
rmap_chr = lapply(1:16, function(i) rbind(data.frame(chr = paste0("Group", i),
                                                     start = 0,
                                                     end = 1,
                                                     cM_Mb = 0,
                                                     off_end_map = TRUE), # add a first bp position to each chr at map pos = 0
                                          rmap_chr0[[i]]) %>%
  mutate(., rpos = cumsum((end - start)*cM_Mb/1000000)) %>% # get cumulative map position for the end of each window
    rename(pos = end) %>%
    dplyr::select(chr, pos, rpos)) 
# dividing by 10^6 is for the Mb -> bp rate conversion

interp_rmap_chr <- lapply(rmap_chr, function(w)
                          approxfun(x = w$pos, y = w$rpos, method = "linear"))
# simple check:
#interp_rmap_chr[[1]](c(1, 5000, 10000))

# read in variant sites and calculate map position for each: (takes ~ 1min)
for (i in 1:16){
  read.table(paste0("results/combined_sept19/variant_sites/Group", i, ".var.sites"),
                  header = F, stringsAsFactors = F,
                  sep = "\t") %>%
    data.table::setnames(c("scaffold", "pos", "major", "minor")) %>%
    mutate(rpos = round(interp_rmap_chr[[i]](pos), 16)) %>% # round to # decimal places
    dplyr::select(scaffold, pos, major, minor, rpos) %>%
    write.table(paste0("results/combined_sept19/variant_sites/Group", i, ".rpos"),
                quote = F, col.names = F, row.names = F, sep = "\t")
}

