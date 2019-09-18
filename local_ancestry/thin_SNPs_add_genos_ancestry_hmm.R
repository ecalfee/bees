library(dplyr)
library(tidyr)

# script to thin SNPs and create sites and A/C/M counts part of input files for ancestry_hmm

ACM <- c("A", "C", "M")
PREFIX = "combined_sept19"
dir.create(paste0("results/SNPs/", PREFIX), recursive = T)
# repeat for each chromosome:
for (chr in paste0("Group", 1:16)){
  sites <- read.table(paste0("../geno_lik_and_SNPs/results/", PREFIX, "/variant_sites/", chr, ".rpos"),
                      header = F, stringsAsFactors = F, sep = "\t") %>%
    data.table::setnames(c("scaffold", "pos", "major", "minor", "rpos"))
  ACM_GL <- lapply(ACM, function(a) 
    read.table(paste0("../geno_lik_and_SNPs/results/", PREFIX, "/genotypes/", a, "/", chr, ".mafs.gz"),
               header = T, stringsAsFactors = F, sep = "\t") %>%
      dplyr::select(chromo, position, knownEM, nInd) %>%
      data.table::setnames(c("scaffold", "pos", paste0("freq_", a), paste0("nInd_", a))) %>%
      left_join(sites, ., by = c("scaffold", "pos")))
  GL <- cbind(ACM_GL[[1]], ACM_GL[[2]][ , 6:7], ACM_GL[[3]][ , 6:7])

  ACM_genos0 <- lapply(ACM, function(a) 
    read.table(paste0("../geno_lik_and_SNPs/results/", PREFIX, "/genotypes/", a, "/", chr, ".geno.gz"),
               header = F, stringsAsFactors = F, sep = "\t", na.strings = "-1") %>%
      rename(scaffold = V1, pos = V2) %>%
      left_join(sites[ , c("scaffold", "pos")], ., by = c("scaffold", "pos"))) # so doesn't skip SNPs with no genotypes called
  ACM_genos1 <- lapply(ACM_genos0, function(acm_genos){ # count allele from genotype calls
    acm_genos[ , c("scaffold", "pos")] %>%
      mutate(n_tot = apply(acm_genos[ , -c(1,2)], 
                           1, 
                           function(row) 2*(sum(!is.na(row))))) %>%
      mutate(n_minor = apply(acm_genos[ , -c(1,2)], 1, function(row) 
        sum(row, na.rm = T))) %>%
      mutate(n_major = n_tot - n_minor)
  })
  CMA_genos <- data.frame(sites, 
                          C_major = ACM_genos1[[2]][ , "n_major"],
                          C_minor = ACM_genos1[[2]][ , "n_minor"],
                          M_major = ACM_genos1[[3]][ , "n_major"],
                          M_minor = ACM_genos1[[3]][ , "n_minor"],
                          A_major = ACM_genos1[[1]][ , "n_major"],
                          A_minor = ACM_genos1[[1]][ , "n_minor"])

  # enrich for more ancestry-informative sites:
  # I filter SNPs for # individuals with data and 1 pop with > 0.3 MAF:
  GL1 <- filter(GL, complete.cases(GL)) %>%
    mutate(freq_tot = (freq_A + freq_C + freq_M)/3) %>%
    mutate(het_tot = 2*freq_tot*(1-freq_tot)) %>%
    filter(het_tot != 0) %>%
    mutate(fst_A = 1 - (2 * freq_A * (1 - freq_A))/het_tot,
           fst_C = 1 - (2 * freq_C * (1 - freq_C))/het_tot,
           fst_M = 1 - (2 * freq_M * (1 - freq_M))/het_tot)
  GL2 <- GL1 %>%
    filter(., nInd_A >= 6, nInd_C >= 6, nInd_M >= 6)
  GL3 <- GL2 %>%
    filter(., (freq_A >= .3 | freq_C >= .3 | freq_M >= .3) &
             (freq_A <= .7 | freq_C <= .7 | freq_M <= .7))


# thin SNPs to reduce within-ancestry LD:
  spacing = .005 # min cM spacing
  keep_snp <- rep(F, nrow(GL3))
  last_rpos <- -100
  for (i in 1:nrow(GL3)) {
    this_rpos <- GL3[i, "rpos"]
    if (this_rpos - last_rpos >= spacing){
      keep_snp[i] <- T
      last_rpos <- this_rpos 
    } # else do nothing  
  }
  GL3_CMA_genos_low <- left_join(GL3[keep_snp, c("scaffold", "pos")],
                                 CMA_genos, by = c("scaffold", "pos")) %>%
    mutate(rdiff = c(1, format(round(diff(rpos)/100, 12), scientific = F))) # 12 decimals, units Morgans, no scientific number
  write.table(GL3_CMA_genos_low[ , c("scaffold", "pos", 
                                     "C_major", "C_minor", 
                                     "M_major", "M_minor", 
                                     "A_major", "A_minor", "rdiff")], 
              paste0("results/SNPs/", PREFIX, "/", chr, ".CMA.genos"),
              col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(GL3_CMA_genos_low[ , 1:4], 
              paste0("results/SNPs/", PREFIX, "/", chr, ".var.sites"),
              col.names = F, row.names = F, quote = F, sep = "\t")
  
}
