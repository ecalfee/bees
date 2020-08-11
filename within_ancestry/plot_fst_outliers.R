library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(gridExtra)
source("../colors.R") # get color palette
source("het_fst_functions.R") # get pi and fst functions
source("calc_frac_snps.R") # calculates snp density (to scale pi)
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


#k = 1000 # window size 1kb
#smoother = .1
# k = 5000 
# smoother = .5 # larger = more smooth
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


# # calculate heterozygosity at each site, with small sample size correction
# hets_small_sample <- lapply(1:nrow(outlier_regions), function(r)
#   do.call(cbind, lapply(4:ncol(hybrid_freqs[[r]]), 
#                         function(i) het_small_sample_correction_nind(p = hybrid_freqs[[r]][ , i],
#                                                                 n_ind = hybrid_ns[[r]][ , i])*frac_snps)))
# 
# #meanA <- zoo::rollmean(na.pad = T, hets_small_sample[ , "A"], k = 1000)
# for (r in 1:nrow(outlier_regions)){
#   colnames(hets_small_sample[[r]]) <- colnames(hybrid_freqs[[r]])[4:ncol(hybrid_freqs[[r]])]
#   hybrid_freqs[[r]]$round_pos <- trunc(hybrid_freqs[[r]]$pos/k)*k+k/2
#   hybrid_freqs[[r]]$snp_n = 1:nrow(hybrid_freqs[[r]])
#   hybrid_freqs[[r]]$group_snp <- trunc(hybrid_freqs[[r]]$snp_n/100) # 100 snps per group
# }
# length(unique(hybrid_freqs[[1]]$round_pos))
# plot heterozygosity with a smooth-over after takign mean in 100 snp intervals    




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
# # nonoutlier region chr3
# slideFst3_windows <- data.frame(starts = seq(min(nonoutliers3$pos), 
#                                              max(nonoutliers3$pos), 
#                                              by = step_size)) %>%
#   mutate(stops = starts + window_size) %>%
#   mutate(pos = (starts + stops)/2) # middle of window



# get allele frequencies from just outlier regions:
c1_sites <- read.table("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group1.rpos",
                       header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))

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

# out1 %>%
#   pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
#   ggplot(.) +
#   geom_smooth(aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison)) +
#   ggtitle("SNP wise Fst") + # add shading for outliers
#   geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
#               filter(outlier_type %in% c("high_AR", "high_CA")) %>%
#               filter(., scaffold == outlier_windows$scaffold[1] & 
#                        start >= outlier_windows$start[1] &
#                        end <= outlier_windows$end[1]),
#             #filter(outlier_type == "high_shared2" & scaffold == outlier_windows$scaffold[1]),
#             aes(xmin = start - buff, xmax = end + buff,
#                 ymin = -Inf, ymax = Inf),
#             alpha = .2)
# ggplot() +
#   geom_rect(data = rename(outliers_all, scaffold = chr) %>% # add shading for outliers
#               filter(outlier_type %in% c("high_AR", "high_CA")) %>%
#               filter(., scaffold == outlier_windows$scaffold[1] & 
#                        start >= outlier_windows$start[1] &
#                        end <= outlier_windows$end[1]),
#             #filter(outlier_type == "high_shared2" & scaffold == outlier_windows$scaffold[1]),
#             aes(xmin = start - buff, xmax = end + buff,
#                 ymin = -Inf, ymax = Inf),
#             alpha = .2) +
#   geom_point(data = pivot_longer(out1, cols = starts_with("fst.A"), names_to = "pop_comparison", values_to = "fst"), 
#              aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison), size = 0.4, alpha = 0.75) +
#   geom_smooth(data = pivot_longer(out1, cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst"), 
#               aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison))
# 
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
              filter(., scaffold == outlier_windows$scaffold[1] & 
                       start >= outlier_windows$start[1] &
                       end <= outlier_windows$end[1]),
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
# # raw Fst each SNP
# out11 %>%
#   pivot_longer(., cols = starts_with("fst"), names_to = "pop_comparison", values_to = "fst") %>%
#   ggplot(.) +
#   geom_smooth(aes(x = pos, y = fst, color = pop_comparison, fill = pop_comparison)) +
#   ggtitle("SNP wise Fst") + # add shading for outliers
#   geom_rect(data = rename(outliers_all, scaffold = chr) %>%
#               filter(outlier_type == "low_AR") %>%
#               filter(scaffold == "NC_037648.1"),
#             aes(xmin = start - buff, xmax = end + buff,
#                 ymin = -Inf, ymax = Inf),
#             alpha = .2)
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
# slideFst11 %>%
#   ggplot(.) +
#   ggtitle("SNP wise Fst") + # add shading for outliers
#   geom_rect(data = rename(outliers_all, scaffold = chr) %>%
#               filter(outlier_type == "low_AR") %>%
#               filter(scaffold == "NC_037648.1"),
#             aes(xmin = start - buff, xmax = end + buff,
#                 ymin = -Inf, ymax = Inf),
#             alpha = .2) +
#   geom_line(aes(x = pos, y = fst, color = pop_comparison)) +
#   scale_color_manual(values = col_fst, labels = labels_col_fst)
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

# Get 2nd outlier on chr 11 -- chr 11 all allele freqs:
# continent level allele frequencies and number of individuals, n
c11_sites <- read.table("../geno_lik_and_SNPs/results/combined_sept19/variant_sites/Group11.rpos",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "pos", "allele1", "allele2", "rpos"))

c11_freqs <- cbind(read.table("results/combined_sept19/A/allele_freq/Group11/SA.freqs.txt", 
                              header = T, stringsAsFactors = F),
                   read.table("results/combined_sept19/A/allele_freq/Group11/NA.freqs.txt", 
                              header = T, stringsAsFactors = F))
c11_ns <- cbind(read.table("results/combined_sept19/A/allele_freq/Group11/SA.nInd", 
                           header = T, stringsAsFactors = F),
                read.table("results/combined_sept19/A/allele_freq/Group11/NA.nInd", 
                           header = T, stringsAsFactors = F))
c11_ACM_freqs <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group11/", a, ".freqs.txt"),
             header = T, stringsAsFactors = F)))
c11_ACM_ns <- do.call(cbind, lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/combined/allele_freq/Group11/", a, ".nInd"),
             header = T, stringsAsFactors = F))) %>%
  data.table::setnames(paste0(colnames(.), ".n"))


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
              filter(., scaffold == outlier_windows$scaffold[3] & 
                       start >= outlier_windows$start[3] &
                       end <= outlier_windows$end[3]),
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
