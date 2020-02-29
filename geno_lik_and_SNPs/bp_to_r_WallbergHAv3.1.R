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

# read in variant sites and calculate map position for each in cM: (takes ~ 1min)
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

# now do reverse (go from cM position to bp position)
interp_rbp_chr <- lapply(rmap_chr, function(w)
  approxfun(x = w$rpos, y = w$pos, method = "linear"))
# test converting back and forth bp -> cM -> bp:
# interp_rbp_chr[[6]](interp_rmap_chr[[6]](5000)) == 5000

# get chromosome ends (max bp)
genome <- read.table("../data/honeybee_genome/chr.list",
                     header = F, sep = "\t") %>%
  data.table::setnames(c("scaffold", "length", "chr", "LG")) %>%
  mutate(chr_n = 1:16)

# extend map to ends of chromosomes (not relevant to our SNP set, but to make full bins):
tail_rmap <- do.call(rbind, 
                     lapply(rmap_chr, function(r) tail(r, n = 1))) %>%
  left_join(., rmap[ , c("chr", "end", "cM_Mb")], by = c("chr", "pos"="end")) %>%
  left_join(., genome, by = c("chr")) %>%
  filter(length > pos) %>%
  mutate(length_cM = rpos + (length-pos)/10^6*cM_Mb) %>% # extend map from last measured cM/Mb bin to end of sequenced scaffold
  dplyr::select(chr, length, length_cM) %>%
  rename(pos = length, rpos = length_cM) # ok now I have rpos for the end of the few chromosomes longer than the map

# add map extensions
rmap_ext <- lapply(1:16, function(i)
  rbind(rmap_chr[[i]], filter(tail_rmap, chr == paste0("Group", i))))
  
# make extended interpolation functions
# bp to cM position
interp_rmap_ext <- lapply(rmap_ext, function(w)
  approxfun(x = w$pos, y = w$rpos, method = "linear"))
# now do reverse (go from cM position to bp position)
interp_rbp_ext <- lapply(rmap_ext, function(w)
  approxfun(x = w$rpos, y = w$pos, method = "linear"))

# divide genome into 1cM bins
n_bins <- sapply(1:16, function(i) ceiling(max(rmap_ext[[i]]$rpos)))
bins0 <- do.call(rbind, 
                lapply(1:16, function(i) 
                  data.frame(chr = paste0("Group", i),
                             end_cM = 1:n_bins[i],
                             end = c(sapply(1:(n_bins[i]-1), function(x) round(interp_rbp_ext[[i]](x))),
                                     genome$length[genome$chr_n == i])) %>% # last bin ends at end of chromosome
                             mutate(., start = c(0, end[1:(n_bins[i]-1)]))))
# name 1cM bins
bins <- bins0 %>%
  left_join(., genome[ , c("chr", "scaffold")], by = "chr") %>%
  mutate(bin_1cM = paste(scaffold, end_cM - 1, sep = "-")) %>%
  dplyr::select(scaffold, start, end, chr, bin_1cM) %>%
  mutate(region = paste0(scaffold, ":", start, "-", end))

# write results
write.table(bins, "results/1cM_bins.bed", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(bins$region, "results/1cM_bins.regions",
            col.names = F, row.names = F, quote = F)
write.table(bins$bin_1cM, "results/1cM_bins.names",
            col.names = F, row.names = F, quote = F)
