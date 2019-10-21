library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
library(gridExtra)
source("../colors.R") # get color palette
library(zoo)
# plot pi and Fst across largest ancestry outlier regions
# this script takes in within-ancestry allele freqs from allele_freq_within_ancestry_outliers.sh
# and reference A/C/M allele freqs generated per chromosome and mapped to all sites with combine_pop_allele_freqs.R
ACM = c("A", "C", "M")
hybrid_pops <- c("NA.A", "NA.C", "NA.M", "SA.A", "SA.C", "SA.M")

outlier_regions = data.frame(outlier_names = c("outlier1_chr1", "outlier2_chr11"), 
                             ranges = c("NC_037638.1:5000000-15000000", "NC_037648.1:12500000-17500000"),
                             start = c(5000000, 12500000), end = c(15000000, 17500000),
                             scaffold = c("NC_037638.1", "NC_037648.1"),
                             chr = c("Group1", "Group11"),
                             chr_n = c(1, 11), stringsAsFactors = F)
# for outlier 1, get A C and M allele freqs from ref. panels for all SNPs:
outlier_ACM_freqs <- lapply(1:nrow(outlier_regions), function(r)
  read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/",
                    outlier_regions$chr[r], ".var.sites"),
             header = F, stringsAsFactors = F) %>%
    data.table::setnames(c("scaffold", "pos", "ref", "alt")) %>%
    mutate(chr = outlier_regions$chr[r]) %>%
    dplyr::select(chr, scaffold, pos) %>%
    cbind(., do.call(cbind, lapply(ACM, function(i) 
      read.table(paste0("results/combined_sept19/combined/allele_freq/", outlier_regions$chr[r], "/", i, ".freqs.txt"),
                 header = T)))) %>%
    filter(pos >= outlier_regions$start[r] & pos <= outlier_regions$end[r]))
outlier_ACM_ns <- lapply(1:nrow(outlier_regions), function(r)
  read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/",
                    outlier_regions$chr[r], ".var.sites"),
             header = F, stringsAsFactors = F) %>%
    data.table::setnames(c("scaffold", "pos", "ref", "alt")) %>%
    mutate(chr = outlier_regions$chr[r]) %>%
    dplyr::select(chr, scaffold, pos) %>%
    cbind(., do.call(cbind, lapply(ACM, function(i) 
      read.table(paste0("results/combined_sept19/combined/allele_freq/", outlier_regions$chr[r], "/", i, ".nInd"),
                 header = T)))) %>%
    filter(pos >= outlier_regions$start[r] & pos <= outlier_regions$end[r]))
head(outlier_ACM_freqs[[2]])
head(outlier_ACM_ns[[1]])

k = 1000 # window size 1kb
smoother = .1
k = 5000 
smoother = .5 # larger = more smooth
# now read in allele freqs from hybid zones:
hybrid_freqs <- lapply(1:nrow(outlier_regions), function(r)
  cbind(outlier_ACM_freqs[[r]], 
        do.call(cbind,
                lapply(hybrid_pops, function(p)
                  read.table(paste0("results/combined_sept19/outliers/", outlier_regions$outlier_names[r], "/", p, ".mafs.gz"),
                             header = T, stringsAsFactors = F) %>%
                    left_join(outlier_ACM_freqs[[r]], ., by = c("scaffold" = "chromo", "pos" = "position")) %>%
                    dplyr::select(phat) %>%
                    data.table::setnames(p)))))
hybrid_ns <- lapply(1:nrow(outlier_regions), function(r)
  cbind(outlier_ACM_ns[[r]], 
        do.call(cbind,
                lapply(hybrid_pops, function(p)
                  read.table(paste0("results/combined_sept19/outliers/", outlier_regions$outlier_names[r], "/", p, ".mafs.gz"),
                             header = T, stringsAsFactors = F) %>%
                    left_join(outlier_ACM_ns[[r]], ., by = c("scaffold" = "chromo", "pos" = "position")) %>%
                    dplyr::select(nInd) %>%
                    data.table::setnames(p)))))
head(hybrid_freqs[[1]])
head(hybrid_ns[[2]])
dim(hybrid_ns[[1]]) == dim(hybrid_freqs[[1]])
dim(hybrid_ns[[2]]) == dim(hybrid_freqs[[2]])
# calculate heterozygosity at each site, with small sample size correction
# het_small_sample_correction() and frac_snps are loaded from plot_pi_fst_from_allele_freqs
hets_small_sample <- lapply(1:nrow(outlier_regions), function(r)
  do.call(cbind, lapply(4:ncol(hybrid_freqs[[r]]), 
                        function(i) het_small_sample_correction(p = hybrid_freqs[[r]][ , i],
                                                                n = hybrid_ns[[r]][ , i])*frac_snps)))

#meanA <- zoo::rollmean(na.pad = T, hets_small_sample[ , "A"], k = 1000)
for (r in 1:nrow(outlier_regions)){
  colnames(hets_small_sample[[r]]) <- colnames(hybrid_freqs[[r]])[4:ncol(hybrid_freqs[[r]])]
  hybrid_freqs[[r]]$round_pos <- trunc(hybrid_freqs[[r]]$pos/k)*k+k/2
  hybrid_freqs[[r]]$snp_n = 1:nrow(hybrid_freqs[[r]])
  hybrid_freqs[[r]]$group_snp <- trunc(hybrid_freqs[[r]]$snp_n/100) # 100 snps per group
}
length(unique(hybrid_freqs[[1]]$round_pos))
# plot heterozygosity with a smooth-over after takign mean in 100 snp intervals    


for (a in 1:2){
  hets_small_sample[[a]] %>%
    cbind(hybrid_freqs[[a]][ , c("group_snp", "round_pos")], .) %>%
    group_by(round_pos) %>%
    summarise(A = mean(SA.A, na.rm = T), M = mean(SA.M, na.rm = T),
              C = mean(SA.C, na.rm = T)) %>%
    tidyr::gather(., "ancestry", "het", ACM) %>%
    ggplot(., aes(x = round_pos, y = het, color = ancestry)) +
    geom_point(size = .5) +
    xlab("Position (bp)") +
    ylab("pi") +
    scale_color_manual(values = col_ACM) +
    geom_smooth(span = smoother, method = "loess", fill = "black") +
    theme_classic() +
    ggtitle(paste0(outlier_regions$outlier_names[a], " - S. America pi within ancestries"))
  ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[a], "_SA_meanOver", k, "bp.png"),
         height = 4, width = 8, device = "png")
  hets_small_sample[[a]] %>%
    cbind(hybrid_freqs[[a]][ , c("group_snp", "round_pos")], .) %>%
    group_by(round_pos) %>%
    summarise(A = mean(NA.A, na.rm = T), M = mean(NA.M, na.rm = T),
              C = mean(NA.C, na.rm = T)) %>%
    tidyr::gather(., "ancestry", "het", ACM) %>%
    ggplot(., aes(x = round_pos, y = het, color = ancestry)) +
    geom_point(size = .5) +
    xlab("Position (bp)") +
    ylab("pi") +
    scale_color_manual(values = col_ACM) +
    geom_smooth(span = smoother, method = "loess", fill = "black") +
    theme_classic() +
    ggtitle(paste0(outlier_regions$outlier_names[a], " - N. America pi within ancestries"))
  ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[a], "_NA_meanOver", k, "bp.png"),
         height = 4, width = 8, device = "png")
  hets_small_sample[[a]] %>%
    cbind(hybrid_freqs[[a]][ , c("group_snp", "round_pos")], .) %>%
    group_by(round_pos) %>%
    summarise(A = mean(A, na.rm = T), M = mean(M, na.rm = T),
              C = mean(C, na.rm = T)) %>%
    tidyr::gather(., "ancestry", "het", ACM) %>%
    ggplot(., aes(x = round_pos, y = het, color = ancestry)) +
    geom_point(size = .5) +
    xlab("Position (bp)") +
    ylab("pi") +
    scale_color_manual(values = col_ACM) +
    geom_smooth(span = smoother, method = "loess", fill = "black") +
    theme_classic() +
    ggtitle(paste0(outlier_regions$outlier_names[a], " - pi within ancestries ref. panel"))
  ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[a], "_ACM_meanOver", k, "bp.png"),
         height = 4, width = 8, device = "png")
}
hets_small_sample[[1]] %>%
  cbind(hybrid_freqs[[1]][ , c("group_snp", "round_pos")], .) %>%
  group_by(round_pos) %>%
  summarise(A = mean(A, na.rm = T), SA.A = mean(SA.A, na.rm = T),
            NA.A = mean(NA.A, na.rm = T)) %>%
  tidyr::gather(., "group", "het", c("A", "SA.A", "NA.A")) %>%
  mutate(group = ifelse(group=="SA.A", "S. America", ifelse(group=="NA.A", "N. America", group))) %>%
  ggplot(., aes(x = round_pos, y = het, color = group)) +
  geom_point(size = .5) +
  xlab("Position (bp)") +
  ylab("pi") +
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "Population") +
  geom_smooth(span = smoother, method = "loess", fill = "black") +
  theme_classic()
ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[1], "_Api_meanOver", k, "bp.png"),
       height = 4, width = 8, device = "png")
hets_small_sample[[2]] %>%
  cbind(hybrid_freqs[[2]][ , c("group_snp", "round_pos")], .) %>%
  group_by(round_pos) %>%
  summarise(M = mean(M, na.rm = T), SA.M = mean(SA.M, na.rm = T),
            NA.M = mean(NA.M, na.rm = T)) %>%
  tidyr::gather(., "group", "het", c("M", "SA.M", "NA.M")) %>%
  mutate(group = ifelse(group=="SA.M", "S. America", ifelse(group=="NA.M", "N. America", group))) %>%
  ggplot(., aes(x = round_pos, y = het, color = group)) +
  geom_point(size = .5) +
  xlab("Position (bp)") +
  ylab("pi") +
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "Population") +
  geom_smooth(span = smoother, method = "loess", fill = "black") +
  theme_classic()
ggsave(paste0("../../bee_manuscript/figures/outlier_", outlier_regions$outlier_names[2], "_Mpi_meanOver", k, "bp.png"),
       height = 4, width = 8, device = "png")


