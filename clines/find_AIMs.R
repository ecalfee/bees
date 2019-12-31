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

chr = 1
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
GL <- cbind(ACM_GL[[1]], ACM_GL[[2]][ , 6:7], ACM_GL[[3]][ , 6:7])
GL %>%
  filter(abs(freq_A - freq_M) > .95 & abs(freq_A - freq_C) > .95) %>%
  filter(nInd_A >= 5 & nInd_C >= 5 & nInd_M >= 5) %>%
  #dim()
  #head()
  pivot_longer(data = ., cols = c("freq_A", "freq_C", "freq_M"), 
               names_to = "pop", values_to = "freq") %>%
  #filter(pos > 7.75*10^6 & pos < 1.25*10^7) %>%
  ggplot(., aes(x = pos, y = freq, color = pop)) +
  geom_point()
# filter for ancestry informative markers
AIMs_A = GL %>%
  filter(abs(freq_A - freq_M) > .95 & abs(freq_A - freq_C) > .95) %>%
  filter(nInd_A >= 5 & nInd_C >= 5 & nInd_M >= 5)
write.table(AIMs_A, paste0("results/AIMs/A/Group", chr, ".ACM.freqs"),
            col.names = T, row.names = F, quote = F, sep = "\t")
# write sites file for ANGSD
AIMs_A %>%
  dplyr::select(c("scaffold", "pos", "major", "minor")) %>%
  write.table(., paste0("results/AIMs/A/Group", chr, ".var.sites"),
              col.names = F, row.names = F, quote = F, sep = "\t")


# now chr 11, M informative sites:
chr = 11
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
GL <- cbind(ACM_GL[[1]], ACM_GL[[2]][ , 6:7], ACM_GL[[3]][ , 6:7])
GL %>%
  filter(abs(freq_A - freq_M) > .95 & abs(freq_M - freq_C) > .95) %>%
  filter(nInd_A >= 5 & nInd_C >= 5 & nInd_M >= 5) %>%
  #dim()
  #head()
  pivot_longer(data = ., cols = c("freq_A", "freq_C", "freq_M"), 
               names_to = "pop", values_to = "freq") %>%
  ggplot(., aes(x = pos, y = freq, color = pop)) +
  geom_point()
# filter for ancestry informative markers
AIMs_M = GL %>%
  filter(abs(freq_A - freq_M) > .95 & abs(freq_M - freq_C) > .95) %>%
  filter(nInd_A >= 5 & nInd_C >= 5 & nInd_M >= 5)
dim(AIMs_M)
write.table(AIMs_M, paste0("results/AIMs/M/Group", chr, ".ACM.freqs"),
            col.names = T, row.names = F, quote = F, sep = "\t")
# write sites file for ANGSD
AIMs_M %>%
  dplyr::select(c("scaffold", "pos", "major", "minor")) %>%
  write.table(., paste0("results/AIMs/M/Group", chr, ".var.sites"),
              col.names = F, row.names = F, quote = F, sep = "\t")

