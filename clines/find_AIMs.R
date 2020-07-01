library(dplyr)
library(tidyr)
library(ggplot2)
# identify ancestry informative sites (AIMs),
# e.g. A-ancestry aim has > 0.95 freq diff with both other ancestries (C & M),
# where all ancestry panels have at least 5 individuals with data.
# just printed output files for two chromosomes with outliers - chr 1 and 11
# output: .sites file with AIM positions (for ANGSD input) and .ACM.freqs file that stores sites & freqs


# load allele freqs A/C/M at all sites:
ACM <- c("A", "C", "M")

find_AIMs_by_chr <- function(chr){
  sites <- read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group",
                             chr, ".rpos"),
                      header = F, stringsAsFactors = F, sep = "\t") %>%
    data.table::setnames(c("scaffold", "pos", "major", "minor", "rpos"))
  ACM_GL <- lapply(ACM, function(a)
    read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19/genotypes/", a, "/Group", chr, ".mafs.gz"),
               header = T, stringsAsFactors = F, sep = "\t") %>%
      dplyr::select(chromo, position, knownEM, nInd) %>%
      data.table::setnames(c("scaffold", "pos", paste0("freq_", a), paste0("nInd_", a))) %>%
      left_join(sites, ., by = c("scaffold", "pos")))
  GL <- cbind(ACM_GL[[1]], ACM_GL[[2]][ , 6:7], ACM_GL[[3]][ , 6:7]) %>%
    filter(nInd_A >= 5 & nInd_C >= 5 & nInd_M >= 5) # remove SNPs with low amounts of data


  # filter for ancestry informative markers
  # A ancestry:
  AIMs_A = GL %>%
    filter(abs(freq_A - freq_M) > .95 & abs(freq_A - freq_C) > .95)

  write.table(AIMs_A, paste0("results/AIMs/A/Group", chr, ".ACM.freqs"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  # write sites file for ANGSD
  AIMs_A %>%
    dplyr::select(c("scaffold", "pos", "major", "minor")) %>%
    write.table(., paste0("results/AIMs/A/Group", chr, ".var.sites"),
                col.names = F, row.names = F, quote = F, sep = "\t")
  # C and M:
  # filter for ancestry informative markers
  AIMs_C = GL %>%
    filter(abs(freq_C - freq_M) > .95 & abs(freq_C - freq_A) > .95)
  write.table(AIMs_C, paste0("results/AIMs/C/Group", chr, ".ACM.freqs"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  # write sites file for ANGSD
  AIMs_C %>%
    dplyr::select(c("scaffold", "pos", "major", "minor")) %>%
    write.table(., paste0("results/AIMs/C/Group", chr, ".var.sites"),
                col.names = F, row.names = F, quote = F, sep = "\t")
  # filter for ancestry informative markers
  AIMs_M = GL %>%
    filter(abs(freq_M - freq_A) > .95 & abs(freq_M - freq_C) > .95)
  write.table(AIMs_M, paste0("results/AIMs/M/Group", chr, ".ACM.freqs"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  # write sites file for ANGSD
  AIMs_M %>%
    dplyr::select(c("scaffold", "pos", "major", "minor")) %>%
    write.table(., paste0("results/AIMs/M/Group", chr, ".var.sites"),
                col.names = F, row.names = F, quote = F, sep = "\t")

}
# create output directories if they don't already exist:
for (a in ACM){
  if(!dir.exists(file.path("results/AIMs/", a))){
    dir.create(file.path("results/AIMs/", a))
  }
}

# find AIMs for regions of interest:
find_AIMs_by_chr(chr = 1)
find_AIMs_by_chr(chr = 11)
# note: AIMs for A are the least common (becasue it's the most diverse lineage)
lapply(c(1:16), function(x) find_AIMs_by_chr(chr = x)) # make files for all the others too
