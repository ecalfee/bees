library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(gridExtra)
source("../colors.R") # get color palette
source("het_fst_functions.R") # get pi and fst functions
library(zoo)
library(slider) # for rolling window averages by index (bp)
# plot pi and Fst across largest ancestry outlier regions
# this script takes in within-ancestry allele freqs from allele_freq_within_ancestry_outliers.sh
# and reference A/C/M allele freqs generated per chromosome and mapped to all sites with combine_pop_allele_freqs.R
ACM = c("A", "C", "M")
hybrid_pops <- c("NA.A", "NA.C", "NA.M", "SA.A", "SA.C", "SA.M")
outliers_all = read.table("../functional_analysis/results/outlier_regions/all.bed",
                           header = T, sep = "\t") # list of all outlier peaks, more precisely

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
                    mutate(phat = ifelse(phat > 1, 1, phat)) %>% # fixes weird ANGSD glitch that allele freq estimates can be slightly > 1
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
hets_small_sample <- lapply(1:nrow(outlier_regions), function(r)
  do.call(cbind, lapply(4:ncol(hybrid_freqs[[r]]), 
                        function(i) het_small_sample_correction_nind(p = hybrid_freqs[[r]][ , i],
                                                                n_ind = hybrid_ns[[r]][ , i])*frac_snps)))

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

##### ------------------------- plot SNP-wise Fst across outlier regions -------------------------- ######
outlier_windows <- bind_rows(data.frame(chr = "Group1", scaffold = "NC_037638.1", start = 1.025*10^7, end = 1.225*10^7, ancestry = "A", stringsAsFactors = F),
                  data.frame(chr = "Group11", scaffold = "NC_037648.1", start = 1.3*10^7, end = 1.6*10^7, ancestry = "M", stringsAsFactors = F),
                  data.frame(chr = "Group11", scaffold = "NC_037648.1", start = 1.5*10^6, end = 2.5*10^6, ancestry = "A", stringsAsFactors = F))
# buffer to make outlier shading show up:
buff = 0
# windowed Fst:
# summarise mean SNP-wise Fst using sliding means
step_size = 10^3
window_size = 5*10^4
min_sites = 10 # minimum sites with data to calculate Fst
# sliding windows:
# outliers chr 1
slideFst1_windows <- data.frame(starts = seq(outlier_windows$start[1], 
                                             outlier_windows$end[1], 
                                             by = step_size)) %>%
  mutate(stops = starts + window_size) %>%
  mutate(pos = (starts + stops)/2) # middle of window
# high M outlier chr 11
slideFst11_windows <- data.frame(starts = seq(outlier_windows$start[2], 
                                              outlier_windows$end[2], 
                                              by = step_size)) %>%
  mutate(stops = starts + window_size) %>%
  mutate(pos = (starts + stops)/2) # middle of window
# high shared A window chr 11
slideFst11_A_windows <- data.frame(starts = seq(outlier_windows$start[3], 
                                              outlier_windows$end[3], 
                                              by = step_size)) %>%
  mutate(stops = starts + window_size) %>%
  mutate(pos = (starts + stops)/2) # middle of window
# nonoutlier region chr3
slideFst3_windows <- data.frame(starts = seq(min(nonoutliers3$pos), 
                                             max(nonoutliers3$pos), 
                                             by = step_size)) %>%
  mutate(stops = starts + window_size) %>%
  mutate(pos = (starts + stops)/2) # middle of window



# get allele frequencies from just outlier regions:
out1 <- hybrid_freqs[[1]] %>%
  left_join(., hybrid_ns[[1]],
            by = c("chr", "scaffold", "pos"),
            suffix = c("", ".n")) %>%
  filter(., scaffold == outlier_windows$scaffold[1] & 
           pos > outlier_windows$start[1] &
           pos < outlier_windows$end[1]) %>%
  left_join(., c1_sites, by = c("scaffold", "pos")) %>%
  filter(!is.na(NA.A.n) & !is.na(SA.A.n) & !is.na(A.n)) %>% # not sure what causes small # of NAs
  mutate(fst.A.NA_SA = sapply(1:nrow(.), function(i) hudson_Fst(p1 = NA.A[i], p2 = SA.A[i], 
                                                                n_ind1 = NA.A.n[i], n_ind2 = SA.A.n[i], 
                                                                min_ind = 2)),
         fst.A.NA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = NA.A[i], p2 = A[i], 
                                                               n_ind1 = NA.A.n[i], n_ind2 = A.n[i], 
                                                               min_ind = 2)),
         fst.A.SA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = SA.A[i], p2 = A[i], 
                                                               n_ind1 = SA.A.n[i], n_ind2 = A.n[i], 
                                                               min_ind = 2)),
         fst.M_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M[i], p2 = A[i], 
                                                               n_ind1 = M.n[i], n_ind2 = A.n[i], 
                                                               min_ind = 2)),
         fst.C_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = A[i], 
                                                               n_ind1 = C.n[i], n_ind2 = A.n[i], 
                                                               min_ind = 2)),
         fst.C_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = M[i], 
                                                            n_ind1 = C.n[i], n_ind2 = M.n[i], 
                                                            min_ind = 2))) %>%
  mutate(pi.A = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = A[i], n_ind = A.n[i], min_ind = 2)),
         pi.C = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = C[i], n_ind = C.n[i], min_ind = 2)),
         pi.M = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = M[i], n_ind = M.n[i], min_ind = 2)),
         pi.A.SA = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = SA.A[i], n_ind = SA.A.n[i], min_ind = 2)),
         pi.A.NA = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = NA.A[i], n_ind = NA.A.n[i], min_ind = 2)))

out1 %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
  ggplot(.) +
  geom_smooth(aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison)) +
  ggtitle("SNP wise Fst") + # add shading for outliers
  geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
              filter(outlier_type %in% c("high_AR", "high_CA")) %>%
              filter(., scaffold == outlier_windows$scaffold[1] & 
                       start >= outlier_windows$start[1] &
                       end <= outlier_windows$end[1]),
            #filter(outlier_type == "high_shared2" & scaffold == outlier_windows$scaffold[1]),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2)
ggplot() +
  geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
              filter(outlier_type %in% c("high_AR", "high_CA")) %>%
              filter(., scaffold == outlier_windows$scaffold[1] & 
                       start >= outlier_windows$start[1] &
                       end <= outlier_windows$end[1]),
            #filter(outlier_type == "high_shared2" & scaffold == outlier_windows$scaffold[1]),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_point(data = pivot_longer(out1, cols = starts_with("fst.A"), names_to = "pop_comparison", values_to = "fst"), 
             aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison), size = 0.4, alpha = 0.75) +
  geom_smooth(data = pivot_longer(out1, cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst"), 
              aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison))

slideFst1 <- bind_cols(slideFst1_windows,
                        apply(select(out1, starts_with("fst")), 2,
                              function(fst) unlist(slider::hop_index(.x = fst, 
                                                                     .i = out1$pos, 
                                                                     .starts = slideFst1_windows$starts,
                                                                     .stops = slideFst1_windows$stops,
                                                                     .f = ~ifelse(sum(!is.na(.x)) < min_sites, NA,
                                                                                  mean(.x, na.rm = T))))) %>%
                          data.frame()) %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst")

ggplot(data = slideFst1) +
  geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
              filter(outlier_type %in% c("high_AR", "high_CA")) %>%
              #mutate(outlier = ifelse(outlier_type == "high_AR", "S. America", "N. America")) %>%
              filter(., scaffold == outlier_windows$scaffold[1] & 
                       start >= outlier_windows$start[1] &
                       end <= outlier_windows$end[1]),
            #filter(outlier_type == "high_shared2" & scaffold == outlier_windows$scaffold[1]),
            aes(xmin = (start - buff)/10^6, xmax = (end + buff)/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_line(aes(x = pos/10^6, y = fst, color = pop_comparison)) +
  ylab(expression(F[ST])) +
  scale_color_manual(name = element_blank(), values = col_fst, labels = labels_col_fst) +
  theme_classic() +
  ylim(c(-.05, 1)) +
  xlab("Chromosome 1 position (Mbp)") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))
ggsave("plots/Fst_across_chr1_outliers.png",
       device = "png", dpi = 600,
       height = 4, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures/Fst_across_chr1_outliers.png",
       device = "png", dpi = 600,
       height = 4, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures_supp/Fst_across_chr1_outliers.tiff",
       device = "tiff", dpi = 600,
       height = 4, width = 5.2, units = "in")

# chr 11 outlier:
# continent level allele frequencies and number of individuals, n
c11_sites <- read.table("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group11.rpos",
                       header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))

out11 <- hybrid_freqs[[2]] %>%
  left_join(., hybrid_ns[[2]],
            by = c("chr", "scaffold", "pos"),
            suffix = c("", ".n")) %>%
  filter(., scaffold == outlier_windows$scaffold[2] & 
           pos > outlier_windows$start[2] &
           pos < outlier_windows$end[2]) %>%
  left_join(., c11_sites, by = c("scaffold", "pos")) %>%
  filter(!is.na(NA.M.n) & !is.na(SA.M.n) & !is.na(M.n)) %>% # later find out what causes these NA's!!
  mutate(fst.M.NA_SA = sapply(1:nrow(.), function(i) hudson_Fst(p1 = NA.M[i], p2 = SA.M[i], 
                                                                n_ind1 = NA.M.n[i], n_ind2 = SA.M.n[i], 
                                                                min_ind = 2)),
         fst.M.NA_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = NA.M[i], p2 = M[i], 
                                                               n_ind1 = NA.M.n[i], n_ind2 = M.n[i], 
                                                               min_ind = 2)),
         fst.M.SA_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = SA.M[i], p2 = M[i], 
                                                               n_ind1 = SA.M.n[i], n_ind2 = M.n[i], 
                                                               min_ind = 2)),
         fst.M_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M[i], p2 = A[i], 
                                                            n_ind1 = M.n[i], n_ind2 = A.n[i], 
                                                            min_ind = 2)),
         fst.C_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = A[i], 
                                                            n_ind1 = C.n[i], n_ind2 = A.n[i], 
                                                            min_ind = 2)),
         fst.C_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = M[i], 
                                                            n_ind1 = C.n[i], n_ind2 = M.n[i], 
                                                            min_ind = 2))) %>%
  mutate(pi.A = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = A[i], n_ind = A.n[i], min_ind = 2)),
         pi.C = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = C[i], n_ind = C.n[i], min_ind = 2)),
         pi.M = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = M[i], n_ind = M.n[i], min_ind = 2)),
         pi.A.SA = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = SA.A[i], n_ind = SA.A.n[i], min_ind = 2)),
         pi.A.NA = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = NA.A[i], n_ind = NA.A.n[i], min_ind = 2)),
         pi.M.SA = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = SA.M[i], n_ind = SA.M.n[i], min_ind = 2)),
         pi.M.NA = sapply(1:nrow(.), function(i) het_small_sample_correction_nind(p = NA.M[i], n_ind = NA.M.n[i], min_ind = 2)))
# raw Fst each SNP
out11 %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
  ggplot(.) +
  geom_smooth(aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison)) +
  ggtitle("SNP wise Fst") + # add shading for outliers
  geom_rect(data = rename(outliers_all, scaffold = chr) %>%
              filter(outlier_type == "low_AR") %>%
              filter(scaffold == "NC_037648.1"),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2)
# sliding window
slideFst11 <- bind_cols(slideFst11_windows,
                        apply(select(out11, starts_with("fst")), 2,
                              function(fst) unlist(slider::hop_index(.x = fst, 
                                                                     .i = out11$pos, 
                                                                     .starts = slideFst11_windows$starts,
                                                                     .stops = slideFst11_windows$stops,
                                                                     .f = ~ifelse(sum(!is.na(.x)) < min_sites, NA,
                                                                                  mean(.x, na.rm = T))))) %>%
                          data.frame()) %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst")
slideFst11 %>%
  ggplot(.) +
  ggtitle("SNP wise Fst") + # add shading for outliers
  geom_rect(data = rename(outliers_all, scaffold = chr) %>%
              filter(outlier_type == "low_AR") %>%
              filter(scaffold == "NC_037648.1"),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_line(aes(x = pos, y = fst, color = pop_comparison)) +
  scale_color_manual(values = col_fst, labels = labels_col_fst)
ggplot(data = slideFst11) +
  geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
              filter(outlier_type %in% c("low_AR")) %>%
              filter(., scaffold == outlier_windows$scaffold[2] & 
                       start >= outlier_windows$start[2] &
                       end <= outlier_windows$end[2]),
            aes(xmin = (start - buff)/10^6, xmax = (end + buff)/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_line(aes(x = pos/10^6, y = fst, color = pop_comparison)) +
  ylab(expression(F[ST])) +
  scale_color_manual(name = element_blank(), values = col_fst, labels = labels_col_fst) +
  theme_classic() +
  ylim(c(-.05, 1)) +
  xlab("Chromosome 11 position (Mbp)") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))
ggsave("plots/Fst_across_chr11_highM_outlier.png",
       device = "png", dpi = 600,
       height = 4, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures/Fst_across_chr11_highM_outlier.png",
       device = "png", dpi = 600,
       height = 4, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures_supp/Fst_across_chr11_highM_outlier.tiff",
       device = "tiff", dpi = 600,
       height = 4, width = 5.2, units = "in")

# ratio of averages alternative:
pairs_A <- data.frame(pop1 = c("A", "A", "C", "A", "A", "NA.A"), 
                      pop2 = c("M", "C", "M", "NA.A", "SA.A", "SA.A"), stringsAsFactors = F)
rownames(pairs_A) <- paste("FST", pairs_A$pop1, pairs_A$pop2, sep = "_")
pairs_M <- data.frame(pop1 = c("A", "A", "C", "M", "M", "NA.M"), 
                      pop2 = c("M", "C", "M", "NA.M", "SA.M", "SA.M"), stringsAsFactors = F)
rownames(pairs_M) <- paste("FST", pairs_M$pop1, pairs_M$pop2, sep = "_")

slideFst1_ratio <- bind_cols(slideFst1_windows,
                             apply(pairs_A, 1, function(r)
                               unlist(slider::hop_index(.x = 1:nrow(out1), 
                                                        .i = out1$pos, 
                                                        .starts = slideFst1_windows$starts,
                                                        .stops = slideFst1_windows$stops,
                                                        .f = ~hudson_Fst_window(p1 = out1[.x, r["pop1"]], 
                                                                                p2 = out1[.x, r["pop2"]], 
                                                                                n_ind1 = out1[.x, paste0(r["pop1"], ".n")], 
                                                                                n_ind2 = out1[.x, paste0(r["pop2"], ".n")], 
                                                                                min_ind = 2, min_snps = 10)))) %>%
                               data.frame(., stringsAsFactors = F) %>%
                               data.table::setnames(rownames(pairs_A)))

ggplot(data = pivot_longer(slideFst1_ratio, cols = starts_with("FST"), names_to = "pop_comparison", values_to = "fst")) +
  geom_line(aes(x = pos, y = fst, color = pop_comparison)) +
  ggtitle("SNP wise Fst - ratio") + # add shading for outliers
  geom_rect(data = rename(outliers_all, scaffold = chr) %>%
              filter(outlier_type == "high_shared2") %>%
              filter(scaffold == "NC_037638.1"),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2)

slideFst11_ratio <- bind_cols(slideFst11_windows,
                              apply(pairs_M, 1, function(r)
                                unlist(slider::hop_index(.x = 1:nrow(out11), 
                                                         .i = out11$pos, 
                                                         .starts = slideFst11_windows$starts,
                                                         .stops = slideFst11_windows$stops,
                                                         .f = ~hudson_Fst_window(p1 = out11[.x, r["pop1"]], 
                                                                                 p2 = out11[.x, r["pop2"]], 
                                                                                 n_ind1 = out11[.x, paste0(r["pop1"], ".n")], 
                                                                                 n_ind2 = out11[.x, paste0(r["pop2"], ".n")], 
                                                                                 min_ind = 2, min_snps = 10)))) %>%
                                data.frame(., stringsAsFactors = F) %>%
                                data.table::setnames(rownames(pairs_M)))
# different ascertainment -- must be polymorphic in both pops in the contrast
slideFst11_ratio2 <- bind_cols(slideFst11_windows,
                              apply(pairs_M, 1, function(r)
                                unlist(slider::hop_index(.x = 1:nrow(out11), 
                                                         .i = out11$pos, 
                                                         .starts = slideFst11_windows$starts,
                                                         .stops = slideFst11_windows$stops,
                                                         .f = ~hudson_Fst_window2(p1 = out11[.x, r["pop1"]], 
                                                                                 p2 = out11[.x, r["pop2"]], 
                                                                                 n_ind1 = out11[.x, paste0(r["pop1"], ".n")], 
                                                                                 n_ind2 = out11[.x, paste0(r["pop2"], ".n")], 
                                                                                 min_ind = 2, min_snps = 10)))) %>%
                                data.frame(., stringsAsFactors = F) %>%
                                data.table::setnames(rownames(pairs_M)))

ggplot(data = pivot_longer(slideFst11_ratio, cols = starts_with("FST"), names_to = "pop_comparison", values_to = "fst")) +
  geom_line(aes(x = pos, y = fst, color = pop_comparison)) +
  ggtitle("SNP wise Fst - ratio") + # add shading for outliers
  geom_rect(data = rename(outliers_all, scaffold = chr) %>%
              filter(outlier_type == "low_AR") %>%
              filter(scaffold == "NC_037648.1"),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2)

ggplot(data = pivot_longer(slideFst11_ratio2, cols = starts_with("FST"), names_to = "pop_comparison", values_to = "fst")) +
  geom_line(aes(x = pos, y = fst, color = pop_comparison)) +
  ggtitle("SNP wise Fst - ratio, ascertainment both") + # add shading for outliers
  geom_rect(data = rename(outliers_all, scaffold = chr) %>%
              filter(outlier_type == "low_AR") %>%
              filter(scaffold == "NC_037648.1"),
            aes(xmin = start - buff, xmax = end + buff,
                ymin = -Inf, ymax = Inf),
            alpha = .2)


# quick background Fst using chr 3 (no outliers):
# continent level allele frequencies and number of individuals, n
c3_A_freqs <- cbind(read.table("results/combined_sept19/A/allele_freq/Group3/SA.freqs.txt", 
                             header = T, stringsAsFactors = F),
                  read.table("results/combined_sept19/A/allele_freq/Group3/NA.freqs.txt", 
                             header = T, stringsAsFactors = F))
c3_A_ns <- cbind(read.table("results/combined_sept19/A/allele_freq/Group3/SA.nInd", 
                          header = T, stringsAsFactors = F),
               read.table("results/combined_sept19/A/allele_freq/Group3/NA.nInd", 
                          header = T, stringsAsFactors = F))
c3_M_freqs <- cbind(read.table("results/combined_sept19/M/allele_freq/Group3/SA.freqs.txt", 
                               header = T, stringsAsFactors = F),
                    read.table("results/combined_sept19/M/allele_freq/Group3/NA.freqs.txt", 
                               header = T, stringsAsFactors = F))
c3_M_ns <- cbind(read.table("results/combined_sept19/M/allele_freq/Group3/SA.nInd", 
                            header = T, stringsAsFactors = F),
                 read.table("results/combined_sept19/M/allele_freq/Group3/NA.nInd", 
                            header = T, stringsAsFactors = F))
c3_sites <- read.table("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group3.rpos",
                       header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))
c3_ACM_freqs <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group3/", a, ".freqs.txt"),
             header = T, stringsAsFactors = F)))
c3_ACM_ns <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group3/", a, ".nInd"),
             header = T, stringsAsFactors = F))) %>%
  data.table::setnames(paste0("n.", colnames(.)))
nonoutliers3 <- cbind(c3_sites, c3_ACM_freqs, c3_ACM_ns,
                  data.frame(A.NA = c3_A_freqs$NA., A.SA = c3_A_freqs$SA, n.A.NA = c3_A_ns$NA., n.A.SA = c3_A_ns$SA,
                  M.NA = c3_M_freqs$NA., M.SA = c3_M_freqs$SA, n.M.NA = c3_M_ns$NA., n.M.SA = c3_M_ns$SA)) %>%
  mutate(fst.A.NA_SA = sapply(1:nrow(.), function(i) hudson_Fst(p1 = A.NA[i], p2 = A.SA[i], 
                                                                n_ind1 = n.A.NA[i], n_ind2 = n.A.SA[i], 
                                                                min_ind = 2)),
         fst.A.NA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = A.NA[i], p2 = A[i], 
                                                               n_ind1 = n.A.NA[i], n_ind2 = n.A[i], 
                                                               min_ind = 2)),
         fst.A.SA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = A.SA[i], p2 = A[i], 
                                                               n_ind1 = n.A.SA[i], n_ind2 = n.A[i], 
                                                               min_ind = 2)),
         fst.M.NA_SA = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M.NA[i], p2 = M.SA[i], 
                                                                n_ind1 = n.M.NA[i], n_ind2 = n.M.SA[i], 
                                                                min_ind = 2)),
         fst.M.NA_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M.NA[i], p2 = M[i], 
                                                               n_ind1 = n.M.NA[i], n_ind2 = n.M[i], 
                                                               min_ind = 2)),
         fst.M.SA_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M.SA[i], p2 = M[i], 
                                                               n_ind1 = n.M.SA[i], n_ind2 = n.M[i], 
                                                               min_ind = 2)))
head(nonoutliers3)
nonoutliers3 %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
  ggplot(., aes(x = fst, color = pop_comparison)) +
  geom_density()

nonoutliers3 %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
  group_by(pop_comparison) %>%
  filter(!is.na(fst)) %>%
  summarise(median = median(fst),
            perc_over0.25 = sum(fst > 0.25)/n(),
            perc_over0.5 = sum(fst > 0.5)/n())

# that's a bit unfair comparison. plot rolling avg. instead for chr3 to compare visually:

slideFst3 <- bind_cols(slideFst3_windows,
                        apply(select(nonoutliers3, starts_with("fst")), 2,
                              function(fst) unlist(slider::hop_index(.x = fst, 
                                                                     .i = nonoutliers3$pos, 
                                                                     .starts = slideFst3_windows$starts,
                                                                     .stops = slideFst3_windows$stops,
                                                                     .f = ~mean(.x, na.rm = T)))) %>%
                          data.frame()) %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst")
pairs_3 <- data.frame(pop1 = c("A", "A", "C", "M", "M", "M.NA", "A", "A", "A.NA"), 
                      pop2 = c("M", "C", "M", "M.NA", "M.SA", "M.SA", "A.NA", "A.SA", "A.SA"), stringsAsFactors = F)
rownames(pairs_3) <- paste("FST", pairs_3$pop1, pairs_3$pop2, sep = "_")
slideFst3_ratio <- bind_cols(slideFst3_windows,
                              apply(pairs_3, 1, function(r)
                                unlist(slider::hop_index(.x = 1:nrow(nonoutliers3), 
                                                         .i = nonoutliers3$pos, 
                                                         .starts = slideFst3_windows$starts,
                                                         .stops = slideFst3_windows$stops,
                                                         .f = ~hudson_Fst_window(p1 = nonoutliers3[.x, r["pop1"]], 
                                                                                 p2 = nonoutliers3[.x, r["pop2"]], 
                                                                                 n_ind1 = nonoutliers3[.x, paste0("n.", r["pop1"])], 
                                                                                 n_ind2 = nonoutliers3[.x, paste0("n.", r["pop2"])], 
                                                                                 min_ind = 2, min_snps = 10)))) %>%
                                data.frame(., stringsAsFactors = F) %>%
                                data.table::setnames(rownames(pairs_3)))

ggplot(data = pivot_longer(slideFst3_ratio, cols = starts_with("FST"), names_to = "pop_comparison", values_to = "fst")) +
  geom_line(aes(x = pos, y = fst, color = pop_comparison)) +
  ggtitle("SNP wise Fst - 3 ratio, ascertainment both")

ggplot(slideFst3, aes(x = pos, y = fst, color = pop_comparison)) +
  geom_line()
slideFst3 %>%
  group_by(pop_comparison) %>%
  filter(!is.na(fst)) %>%
  summarise(median = median(fst),
            perc_over0.1 = sum(fst > 0.1)/n(),
            perc_over0.15 = sum(fst > 0.15)/n(),
            perc_over0.25 = sum(fst > 0.25)/n(),
            perc_over0.5 = sum(fst > 0.5)/n())

### --------------- old, not using ----------- #
# continent level allele frequencies and number of individuals, n
c1_freqs <- cbind(read.table("results/combined_sept19/A/allele_freq/Group1/SA.freqs.txt", 
                             header = T, stringsAsFactors = F),
                  read.table("results/combined_sept19/A/allele_freq/Group1/NA.freqs.txt", 
                             header = T, stringsAsFactors = F))
c1_ns <- cbind(read.table("results/combined_sept19/A/allele_freq/Group1/SA.nInd", 
                          header = T, stringsAsFactors = F),
               read.table("results/combined_sept19/A/allele_freq/Group1/NA.nInd", 
                          header = T, stringsAsFactors = F))
c1_sites <- read.table("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group1.rpos",
                       header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))
c1_ACM_freqs <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group1/", a, ".freqs.txt"),
             header = T, stringsAsFactors = F)))
c1_ACM_ns <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group1/", a, ".nInd"),
             header = T, stringsAsFactors = F))) %>%
  data.table::setnames(paste0("n.", colnames(.)))
outlier1 <- cbind(c1_sites, c1_ACM_freqs, c1_ACM_ns,
                  data.frame(A.NA = c1_freqs$NA., A.SA = c1_freqs$SA, n.NA = c1_ns$NA., n.SA = c1_ns$SA)) %>%
  mutate(fst.A.NA_SA = sapply(1:nrow(.), function(i) hudson_Fst(p1 = A.NA[i], p2 = A.SA[i], 
                                                                n_ind1 = n.NA[i], n_ind2 = n.SA[i], 
                                                                min_ind = 2)),
         fst.A.NA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = A.NA[i], p2 = A[i], 
                                                               n_ind1 = n.NA[i], n_ind2 = n.A[i], 
                                                               min_ind = 2)),
         fst.A.SA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = A.SA[i], p2 = A[i], 
                                                               n_ind1 = n.SA[i], n_ind2 = n.A[i], 
                                                               min_ind = 2)))

outlier1 %>%
  filter(., scaffold == outlier_windows$scaffold[1] & 
           pos > outlier_windows$start[1] &
           pos < outlier_windows$end[1]) %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
  ggplot(., aes(x = pos, y = fst, color = pop_comparison)) +
  geom_point()
# oops! the angsd file must have failed -- no frequency calls for A within S. America after the first part of the chromosome
outlier1 %>%
  ggplot(., aes(x = pos, y = n.SA)) + 
  geom_line()


# chr 11 all allele freqs:
# continent level allele frequencies and number of individuals, n
c11_freqs <- cbind(read.table("results/combined_sept19/A/allele_freq/Group11/SA.freqs.txt", 
                             header = T, stringsAsFactors = F),
                  read.table("results/combined_sept19/A/allele_freq/Group11/NA.freqs.txt", 
                             header = T, stringsAsFactors = F))
c11_ns <- cbind(read.table("results/combined_sept19/A/allele_freq/Group11/SA.nInd", 
                          header = T, stringsAsFactors = F),
               read.table("results/combined_sept19/A/allele_freq/Group11/NA.nInd", 
                          header = T, stringsAsFactors = F))
c11_sites <- read.table("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group11.rpos",
                       header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))
c11_ACM_freqs <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group11/", a, ".freqs.txt"),
             header = T, stringsAsFactors = F)))
c11_ACM_ns <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group11/", a, ".nInd"),
             header = T, stringsAsFactors = F))) %>%
  data.table::setnames(paste0(colnames(.), ".n"))
cbind(c11_sites, c11_ns) %>%
  ggplot(., aes(x = pos, y = NA.)) + 
  geom_line()


outlier11_A <- cbind(c11_sites, c11_ACM_freqs, c11_ACM_ns,
                  data.frame(NA.A = c11_freqs$NA., SA.A = c11_freqs$SA, NA.A.n = c11_ns$NA., SA.A.n = c11_ns$SA))
out11_A <- outlier11_A %>%
  filter(., scaffold == outlier_windows$scaffold[3] & 
           pos > outlier_windows$start[3] &
           pos < outlier_windows$end[3]) %>%
  filter(!is.na(NA.A.n) & !is.na(SA.A.n) & !is.na(A.n)) %>% # not sure what causes small # of NAs
  mutate(fst.A.NA_SA = sapply(1:nrow(.), function(i) hudson_Fst(p1 = NA.A[i], p2 = SA.A[i], 
                                                                n_ind1 = NA.A.n[i], n_ind2 = SA.A.n[i], 
                                                                min_ind = 2)),
         fst.A.NA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = NA.A[i], p2 = A[i], 
                                                               n_ind1 = NA.A.n[i], n_ind2 = A.n[i], 
                                                               min_ind = 2)),
         fst.A.SA_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = SA.A[i], p2 = A[i], 
                                                               n_ind1 = SA.A.n[i], n_ind2 = A.n[i], 
                                                               min_ind = 2)),
         fst.M_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M[i], p2 = A[i], 
                                                            n_ind1 = M.n[i], n_ind2 = A.n[i], 
                                                            min_ind = 2)),
         fst.C_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = A[i], 
                                                            n_ind1 = C.n[i], n_ind2 = A.n[i], 
                                                            min_ind = 2)),
         fst.C_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = M[i], 
                                                            n_ind1 = C.n[i], n_ind2 = M.n[i], 
                                                            min_ind = 2)))
slideFst11_A <- bind_cols(slideFst11_A_windows,
                       apply(select(out11_A, starts_with("fst")), 2,
                             function(fst) unlist(slider::hop_index(.x = fst, 
                                                                    .i = out11_A$pos, 
                                                                    .starts = slideFst11_A_windows$starts,
                                                                    .stops = slideFst11_A_windows$stops,
                                                                    .f = ~ifelse(sum(!is.na(.x)) < min_sites, NA,
                                                                                 mean(.x, na.rm = T))))) %>%
                         data.frame()) %>%
  pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst")

ggplot(data = slideFst11_A) +
  geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
              filter(outlier_type %in% c("high_AR", "high_CA")) %>%
              #mutate(outlier = ifelse(outlier_type == "high_AR", "S. America", "N. America")) %>%
              filter(., scaffold == outlier_windows$scaffold[3] & 
                       start >= outlier_windows$start[3] &
                       end <= outlier_windows$end[3]),
            #filter(outlier_type == "high_shared2" & scaffold == outlier_windows$scaffold[1]),
            aes(xmin = (start - buff)/10^6, xmax = (end + buff)/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_line(aes(x = pos/10^6, y = fst, color = pop_comparison)) +
  ylab(expression(F[ST])) +
  scale_color_manual(name = element_blank(), values = col_fst, labels = labels_col_fst) +
  theme_classic() +
  ylim(c(-.05, 1)) +
  xlab("Chromosome 11 position (Mbp)") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))
ggsave("plots/Fst_across_chr11_highA_outlier.png",
       device = "png", dpi = 600,
       height = 4, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures/Fst_across_chr11_highA_outlier.png",
       device = "png", dpi = 600,
       height = 4, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures_supp/Fst_across_chr11_highA_outlier.tiff",
       device = "tiff", dpi = 600,
       height = 4, width = 5.2, units = "in")  

# what are the fst quantiles for A-C, M-C, A-M genome-wide? are any of these outliers?
all_chrs_ACM <- do.call(rbind,
                        lapply(1:16, function(chr){
  s <- read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group", chr, ".rpos"),
                          header = F, stringsAsFactors = F) %>%
    data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))
  s_ACM_freqs <- do.call(cbind, lapply(ACM, function(a) 
    read.table(paste0("results/combined_sept19/combined/allele_freq/Group", chr, "/", a, ".freqs.txt"),
               header = T, stringsAsFactors = F)))
  s_ACM_ns <- do.call(cbind, lapply(ACM, function(a) 
    read.table(paste0("results/combined_sept19/combined/allele_freq/Group", chr, "/", a, ".nInd"),
               header = T, stringsAsFactors = F))) %>%
    data.table::setnames(paste0(colnames(.), ".n"))
  slide_windows <- data.frame(starts = seq(min(s$pos), 
                                           max(s$pos), 
                                           by = step_size)) %>%
    mutate(stops = starts + window_size) %>%
    mutate(pos = (starts + stops)/2)
  s_ACM <- cbind(s, s_ACM_freqs, s_ACM_ns) %>%
    mutate(fst.M_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = M[i], p2 = A[i], 
                                                       n_ind1 = M.n[i], n_ind2 = A.n[i], 
                                                       min_ind = 2)),
           fst.C_A = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = A[i], 
                                                     n_ind1 = C.n[i], n_ind2 = A.n[i], 
                                                     min_ind = 2)),
           fst.C_M = sapply(1:nrow(.), function(i) hudson_Fst(p1 = C[i], p2 = M[i], 
                                                     n_ind1 = C.n[i], n_ind2 = M.n[i], 
                                                     min_ind = 2)))
  s_slide <- bind_cols(slide_windows,
                          apply(select(s_ACM, starts_with("fst")), 2,
                                function(fst) unlist(slider::hop_index(.x = fst, 
                                                                       .i = s_ACM$pos, 
                                                                       .starts = slide_windows$starts,
                                                                       .stops = slide_windows$stops,
                                                                       .f = ~ifelse(sum(!is.na(.x)) < min_sites, NA,
                                                                                    mean(.x, na.rm = T))))) %>%
                            data.frame())
  return(s_slide)
}))
save(all_chrs_ACM, file = "results/all_chrs_ACM_fst_windows.RData")
# how big of outliers are the outliers on chr 11 for A-C-M Fst?
# pretty big outliers
max_fst_vals <- slideFst11 %>%
  group_by(pop_comparison) %>%
  summarise(max = max(fst, na.rm = T)) %>%
  .[1:3,]
sapply(1:nrow(max_fst_vals), function(i) 
  sum(all_chrs_ACM[ , max_fst_vals$pop_comparison[i]] < max_fst_vals$max[i], na.rm = T)/sum(!is.na(all_chrs_ACM[ , max_fst_vals$pop_comparison[i]])))


apply(all_chrs_ACM[ , 4:6], 2, function(x) quantile(x, c(0.9, 0.99, 0.999), na.rm = T))



