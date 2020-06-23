library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
#library(rethinking)
library(Hmisc)
source("../colors.R") # for color palette
source("k_matrix.R") # import useful functions

# this script looks at the effect of recombination rate on the K matrix
# in an effort to disentangle demography and selection
# in creating ancestry covariances between populations, especially the negative covariances observed

# I need to load the local ancestry data from plotLocalAncestryTracts.R:
# objects: A, meta.ind, meta.pop, sites
load("results/A.RData") # c("A", "sites", "meanA")
load("results/meta.RData") # c("meta.ind", "meta.pop", "pops_by_lat", "meta.AR.order.by.lat", "NS_segments")

# load recombination rates:
rmap <- read.table("../data/recomb_map/Wallberg_HAv3.1/map_rates_extended_10kb.bed",
                   header = T,
                   stringsAsFactors = F)

# get genomewide r quintiles:
rmap$r_bin5 <- cut(rmap$cM_Mb,
                   breaks = unique(quantile(c(0, rmap$cM_Mb, max(rmap$cM_Mb)),
                                            p = seq(0, 1, by = .2))),
                   right = T,
                   include.lowest = T)
table(rmap$r_bin5) # good, each bin includes 20% of windows
levels(rmap$r_bin5)

# map recombination rates (and bins/quintiles) onto sites
rerun_sites_r = F
if (rerun_sites_r){
  # get SNP sites where ancestry was called
  sites0 <- read.table("results/SNPs/combined_sept19/chr.var.sites", stringsAsFactors = F,
                       sep = "\t", header = F)[ , 1:2]
  colnames(sites0) <- c("scaffold", "pos")
  chr_lengths <- cbind(read.table("../data/honeybee_genome/chr.names", stringsAsFactors = F),
                       read.table("../data/honeybee_genome/chr.lengths", stringsAsFactors = F)) %>%
    data.table::setnames(c("chr", "scaffold", "chr_length")) %>%
    mutate(chr_n = 1:16) %>%
    mutate(chr_end = cumsum(chr_length)) %>%
    mutate(chr_start = chr_end - chr_length) %>%
    mutate(chr_mid = (chr_start + chr_end)/2)
  
  sites <- left_join(sites0, chr_lengths[ , c("chr", "scaffold", "chr_n", "chr_start")], by = "scaffold") %>%
    mutate(cum_pos = pos + chr_start)
  sites_r <- bedr(
    engine = "bedtools", 
    input = list(a = dplyr::mutate(sites, start = pos - 1, end = pos) %>%  # make 0 index for bedtools
                   dplyr::select(chr, start, end, pos, scaffold, chr_n, chr_start, cum_pos),
                 b = rmap), 
    method = "map", 
    # cM_Mb is col 4 and r_bin5 is col 6
    params = "-g ../data/honeybee_genome/chr_names.lengths -c 4,6 -o collapse", 
    check.chr = F
  ) 
  colnames(sites_r) <- c("chr", "start", "end", "pos", "scaffold", 
                         "chr_n", "chr_start", "cum_pos", "cM_Mb", "r_bin5")
  sites_r$cM_Mb <- as.numeric(sites_r$cM_Mb)
  sites_r$pos <- as.numeric(sites_r$pos)
  sites_r$chr_n <- as.numeric(sites_r$chr_n)
  sites_r$cum_pos <- as.numeric(sites_r$cum_pos)
  sites_r$chr_start <- as.numeric(sites_r$chr_start)
  sites_r$r_bin5_factor <- factor(sites_r$r_bin5, levels = levels(rmap$r_bin5),
                                  ordered = T)
  sites_r$r_bin5 <- as.numeric(sites_r$r_bin5_factor) # translate to 1-5 numbers
  
  #head(sites_r)
  #str(sites_r)
  #table(sites_r$r_bin5_factor)
  #table(is.na(sites_r$r_bin5_factor))
  #levels(sites_r$r_bin5_factor)
  save(sites_r, file = "results/sites_r.RData")
  # also add cM position for each site & save:
  sites_rpos <- sites_r %>%
    left_join(., do.call(rbind, # add recombination position
                         lapply(paste0("Group", 1:16), function(chr) 
                           read.table(paste0("../geno_lik_and_SNPs/results/combined_sept19", 
                                             "/variant_sites/", chr, ".rpos"),
                                      header = F, stringsAsFactors = F, sep = "\t") %>%
                             data.table::setnames(c("scaffold", "pos", "major", "minor", "rpos")) %>%
                             mutate(chr = chr))),
              by = c("scaffold", "pos", "chr")) %>%
    rename(pos_cM = rpos)
  save(sites_rpos, file = "results/sites_rpos.RData")
  
} else{
  load("results/sites_r.RData")
}



# look at mean A ancestry across recombination bins:
cbind(sites_r, meanA) %>%
  group_by(r_bin5_factor) %>%
  summarise(mean = mean(meanA))

# get mean A ancestry for each population for the 5 recombination bins
A_r <- cbind(sites_r, A) %>%
  tidyr::gather(., "population", "A", colnames(A)) %>%
  group_by(r_bin5, r_bin5_factor, population) %>%
  summarise(A = mean(A)) %>%
  data.frame() %>%
  left_join(., meta.pop, by = "population") %>%
  mutate(abs_lat_c = abs(lat) - mean(abs(lat)))
A_r$r_bin5_number = as.numeric(A_r$r_bin5_factor)
A_r$r_bin5_number_c = A_r$r_bin5_number - mean(A_r$r_bin5_number)
#unique(A_r[ , c("r_bin5", "r_bin5_number")]) # in order low to high, good
#str(A_r)
levels(A_r$r_bin5_factor)
A_r %>%
  #filter(r_bin5_factor %in% c("[0,2.92]", "(38.6,66.9]")) %>%
  ggplot(., aes(x = abs(lat), y = A, color = r_bin5_factor)) +
  geom_point() +
  facet_wrap(~zone)

# I test and find that mean cline is the tiniest bit steeper for low recombination rate regions (but high uncertainty(:
# but center of the clines show no consistent pattern
# includes N and S America jointly..maybe better to split them up
m_lat_r_mean_cline <- map2stan( # quadratic approximation of the posterior MAP
  alist(
    A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu[r_bin5_number] + b_lat[r_bin5_number]*abs_lat_c,
    mu[r_bin5_number] ~ dnorm(0, 5),
    theta ~ dcauchy(0, 2), # allows it to sample large #s as needed
    b_lat[r_bin5_number] ~ dnorm(0, 5)
  ),
  data = A_r,
  iter = 2000, warmup = 1000, chains = 2, cores = 2)
write.table(precis(m_lat_r_mean_cline, depth = 2)@output,
            "results/m_lat_r_mean_cline_model_summary.txt",
            row.names = T, col.names = T)
precis(m_lat_r_mean_cline, depth = 2)
precis(m_lat_r_mean_cline, depth = 2)@output
pairs(m_lat_r_mean_cline)
plot(m_lat_r_mean_cline) # there are warnings from STAN about 1 divergent transition, but chains look ok
# HOWEVER, probably not how I want to analyze these with a beta dist...perhaps better just using nls()

# 
# 
# # do we see a significant effect of recombination rate on the shape of the cline?
# # short answer: no
# d_A_r <- cbind(sites_r, A) %>%
#   tidyr::gather(., "population", "A", colnames(A)) %>%
#   left_join(., meta.pop, by = "population") %>%
#   mutate(abs_lat_c = abs(lat) - mean(abs(lat))) %>%
#   mutate(cM_Mb_c = cM_Mb - mean(cM_Mb))
# dim(d_A_r) # very large dataset, I'll do the fit with a smaller portion
# d_A_r_test <- d_A_r[c(T, rep(F, 100)),]
# dim(d_A_r_test)
# table(d_A_r_test$population)
# 
# m_lat_r_bypop <- map( # quadratic approximation of the posterior MAP
#   alist(
#     A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
#     # beta distribution provided by the 'rethinking' package
#     # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
#     # are used instead of shape1 = a and shape2 = b
#     logit(p) <- mu + b_lat*abs_lat_c + b_r*cM_Mb_c + b_r_lat*cM_Mb_c*abs_lat_c,
#     mu ~ dnorm(0, 5),
#     theta ~ dunif(0, 10),
#     b_lat ~ dnorm(0, 5),
#     b_r ~ dnorm(0,5),
#     b_r_lat ~ dnorm(0, 5)
#   ),
#   data = d_A_r_test,
#   start = list(mu = -1, theta = 2, b_lat = 0, b_r = 0, b_r_lat = 0))
# precis(m_lat_r_bypop)
# pairs(m_lat_r_bypop)
# 
# 

# make K matrix for each recombination rate bin
# using the genomewide alpha rather than alpha
# for that recombination rate bin
pop_alpha_genomewide <- apply(A[ , pops_by_lat], 2, mean)
K_byr2 <- lapply(1:5, function(r)
  calcK(ancFreqMatrix = t(A[sites_r$r_bin5 == r, pops_by_lat]),
        alpha = pop_alpha_genomewide[pops_by_lat]))
# # make new plots:
# for (a in 1:5){
#   melt(K_byr2[[a]]) %>%
#     filter(Var1 != Var2) %>% # don't plot diagonal
#     ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#     geom_tile() +
#     xlab("") +
#     ylab("") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     scale_fill_viridis(begin = 0, end = 1, direction = 1,
#                        limits = c(-.005, .01)) +
#     ggtitle(paste0("K matrix covariances rbin ", a, ": ", levels(rmap$r_bin5)[a], " genomewide alpha"))
#   ggsave(paste0("plots/k_matrix_covariance_r", a, "_genomewide_alpha.png"), 
#          height = 6, width = 8, 
#          units = "in", device = "png")
#   melt(cov2cor(K_byr2[[a]])) %>%
#     filter(Var1 != Var2) %>% # don't plot diagonal
#     ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#     geom_tile() +
#     xlab("") +
#     ylab("") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     scale_fill_viridis(begin = 0, end = 1, direction = 1,
#                        limits = c(-.35, .55)) +
#     ggtitle(paste0("K matrix correlations rbin ", a, ": ", levels(rmap$r_bin5)[a], " genomewide alpha"))
#   ggsave(paste0("plots/k_matrix_correlation_r", a, "_genomewide_alpha.png"), 
#          height = 6, width = 8, 
#          units = "in", device = "png")
# }

# calculate mean correlations for different recombination rate bins:
# function get_mean_from_K needs to be loaded from plotLocalAncestryTracts.R:
zAnc_bees <- make_K_calcs(t(A[ , pops_by_lat]))
mean_corr_k <- get_mean_from_K(cov2cor(zAnc_bees$K), m = meta.pop)
mean_corrs0 <- do.call(rbind,
                      lapply(1:5, function(k) get_mean_from_K(cov2cor(K_byr2[[k]]), m = meta.pop) %>%
                               mutate(r_bin5 = k))) %>%
  left_join(., distinct(dplyr::select(sites_r, c("r_bin5", "r_bin5_factor"))), by = "r_bin5")

pair_types <- data.frame(label = c("High-A SA", "Low-A SA vs. High-A SA", "Low-A SA", "Low-A NA vs. High-A SA", "Low-A NA vs. Low-A SA", "Low-A NA"),
                         type = c("ARN_ARN", "ARN_ARS", "ARS_ARS", "CA_ARN",  "CA_ARS",  "CA_CA"),
                         stringsAsFactors = F) %>%
  mutate(label = factor(label, levels = c("High-A SA", "Low-A SA", "Low-A NA", "Low-A NA vs. Low-A SA", "Low-A SA vs. High-A SA", "Low-A NA vs. High-A SA"),
                        ordered = T))
mean_corrs = left_join(mean_corrs0, pair_types, by = "type")

mean_corrs %>%
  write.table(., "results/mean_anc_corr_grouped_by_r.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
mean_corrs_plot <- mean_corrs %>%
  ggplot(., aes(x = label, y = mean_anc_corr, color = r_bin5_factor)) +
  geom_point(alpha = .75) +
  geom_point(data = left_join(mean_corr_k, pair_types, by = "type"), 
             aes(x = label, y = mean_anc_corr), color = "black", pch = 4) +
  ylab("Mean Ancestry Correlation") +
  theme_classic() +
  xlab("") +
  scale_color_viridis_d(name = "Recombination bin (cM/Mb)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mean_corrs_plot
ggsave(paste0("plots/mean_k_corr_by_groups_and_r.png"), 
       height = 4, width = 5.2,
       plot = mean_corrs_plot,
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures/mean_k_corr_by_groups_and_r.png"), 
       height = 4, width = 5.2, dpi = 600,
       plot = mean_corrs_plot,
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures_supp/mean_k_corr_by_groups_and_r.tiff"), 
       height = 4, width = 5.2, dpi = 600,
       plot = mean_corrs_plot,
       units = "in", device = "tiff",
       compression = "lzw", type = "cairo")

# paired t-test for manuscript:
with(filter(mean_corrs, type %in%  c("CA_ARS", "ARN_ARS")) %>% arrange(r_bin5), 
     t.test(mean_anc_corr ~ type, paired = T), alternative = "two.sided")

# is this finding robust across all chromosomes?
K_by_chr <- lapply(1:16, function(chr)
  calcK(t(A[sites_r$chr_n == chr, pops_by_lat]),
               alpha = pop_alpha_genomewide[pops_by_lat]))

mean_corrs_chr <- do.call(rbind,
                       lapply(1:16, function(k) get_mean_from_K(cov2cor(K_by_chr[[k]])) %>%
                                mutate(chr_n = k))) %>%
  left_join(., pair_types, by = "type")
mean_covs_chr <- do.call(rbind,
                          lapply(1:16, function(k) get_mean_from_K(K_by_chr[[k]]) %>%
                                   mutate(chr_n = k))) %>%
  left_join(., pair_types, by = "type")


# plot with bars representing all chromosomes:
p_corrs_chr <- mean_corrs_chr %>%
  ggplot(., aes(x = label, y = mean_anc_corr)) +
  geom_jitter(aes(color = factor(chr_n), shape = factor(chr_n)), width = .2) +
  ylab("Mean Ancestry Correlation") +
  theme_classic() +
  xlab("") +
  stat_summary(data = mean_corrs_chr, #%>%
               #filter(!(chr_n %in% c(1, 11))),
               mapping = aes(x = label, y = mean_anc_corr),
               fun.data = "mean_cl_normal", 
               geom = "errorbar",
               color = "black",
               width = .3) +
  scale_color_viridis_d(name = "Chromosome", labels = 1:16, option = "D") +
  scale_shape_manual(name = "Chromosome",
                     labels = 1:16,
                     values = rep(c(1:3,19), 10)[1:16]) +
  guides(color=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
p_corrs_chr
ggsave(paste0("plots/mean_k_corr_by_groups_and_chr.png"), 
       height = 4, width = 5.2, dpi = 600,
       plot = p_corrs_chr + ggtitle("Ancestry correlations by chromosome"),
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures/mean_k_corr_by_groups_and_chr.png"), 
       height = 4, width = 5.2, dpi = 600,
       plot = p_corrs_chr,
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures_supp/mean_k_corr_by_groups_and_chr.tiff"), 
       height = 4, width = 5.2, dpi = 600,
       plot = p_corrs_chr,
       units = "in", device = "tiff",
       compression = "lzw", type = "cairo")

# plot with bars representing 95% CI from t-test excluding outlier chromosome 1 and 11:
mean_corrs_chr %>%
  ggplot(., aes(x = label, y = mean_anc_corr)) +
  #stat_summary(data = mean_corrs_chr %>%
  #             filter(!(chr_n %in% c(1, 11))),
  #             mapping = aes(x = label, y = mean_anc_corr),
  #             fun.y = "mean", 
  #             geom = "point",
  #             color = "black",
  #             fill = "black",
  #             shape = 23) +
  geom_jitter(aes(color = factor(chr_n), shape = factor(chr_n)), width = .2) +
  ylab("Mean Ancestry Correlation") +
  theme_classic() +
  xlab("") +
  stat_summary(data = mean_corrs_chr %>%
                 filter(!(chr_n %in% c(1, 11))),
               mapping = aes(x = label, y = mean_anc_corr),
               fun.data = "mean_cl_normal", 
               geom = "errorbar",
               color = "black",
               width = .3) +
  scale_color_viridis_d(name = "Chromosome", labels = 1:16, option = "D") +
  scale_shape_manual(name = "Chromosome",
                     labels = 1:16,
                     values = rep(c(1:3,19), 10)[1:16]) +
  guides(color=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Ancestry correlations by chromosome")
ggsave(paste0("plots/mean_k_corr_by_groups_and_chr_outliers_chr1and11.png"), 
       height = 4, width = 6, 
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures/mean_k_corr_by_groups_and_chr_outliers_chr1and11.pdf"), 
       height = 4, width = 6, 
       units = "in", device = "pdf")

# report paired t-test in main text (paired because chromosomes are 'resampled' for all comparisons):
with(filter(mean_corrs_chr, type %in%  c("CA_ARS", "ARN_ARS")) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T), alternative = "two.sided")
# sensitive to including chr 1 and chr 11 (in diff. directions), but result still valid without either large outlier chromosomes
with(filter(mean_corrs_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              !(chr_n %in% c(1, 11))) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))
with(filter(mean_corrs_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              chr_n != 1) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))
with(filter(mean_corrs_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              chr_n != 11) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))

# plot covariances:
# plot with bars representing all chromosomes:
mean_covs_chr %>%
  ggplot(., aes(x = label, y = mean_anc_corr)) +
  geom_jitter(aes(color = factor(chr_n), shape = factor(chr_n)), width = .2) +
  ylab("Mean Ancestry Covariance") +
  theme_classic() +
  xlab("") +
  stat_summary(data = mean_covs_chr %>%
               filter(!(chr_n %in% c(1, 11))),
               mapping = aes(x = label, y = mean_anc_corr),
               fun.data = "mean_cl_normal", 
               geom = "errorbar",
               color = "black",
               width = .3) +
  scale_color_viridis_d(name = "Chromosome", labels = 1:16, option = "D") +
  scale_shape_manual(name = "Chromosome",
                     labels = 1:16,
                     values = rep(c(1:3,19), 10)[1:16]) +
  guides(color=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Ancestry covariances by chromosome")
ggsave(paste0("plots/mean_k_cov_by_groups_and_chr.png"), 
       height = 4, width = 6, 
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures/mean_k_cov_by_groups_and_chr.pdf"), 
       height = 4, width = 6, 
       units = "in", device = "pdf")
# clearly these covariances aren't diff.
with(filter(mean_covs_chr, type %in%  c("CA_ARS", "ARN_ARS")), t.test(mean_anc_corr ~ type), alternative = "two.sided")
with(filter(mean_covs_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              !(chr_n %in% c(1, 11))), t.test(mean_anc_corr ~ type), alternative = "two.sided")
with(filter(mean_covs_chr, type %in%  c("CA_ARS", "ARN_ARS")) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))
with(filter(mean_covs_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              !(chr_n %in% c(1, 11))) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))

# plot standardized covariances instead of correlations: 
# i.e. dividing by sqrt(2*alpha1*(1-alpha1)*2*alpha2*(1-alpha2)) for mean ancestry in pops 1 and 2 
# like we would for extracting drift out of co-varying allele freqs:
standardize_k <- function(K, alpha){
  v = sqrt(2*alpha*(1-alpha))
  K/outer(v, v, "*") # divide to standardize
}
standardize_k(K = matrix(1, 2, 2), alpha = c(.25, .5))
mean_sd_cov_chr <- do.call(rbind,
                          lapply(1:16, function(k) get_mean_from_K(standardize_k(K_by_chr[[k]], 
                                                                                 alpha = pop_alpha_genomewide[colnames(K_by_chr[[k]])])) %>%
                                   mutate(chr_n = k))) %>%
  left_join(., pair_types, by = "type")
mean_sd_cov_chr %>%
  ggplot(., aes(x = label, y = mean_anc_corr)) +
  geom_jitter(aes(color = factor(chr_n), shape = factor(chr_n)), width = .2) +
  ylab("Mean Ancestry Covariance (Standardized)") +
  theme_classic() +
  xlab("") +
  stat_summary(data = mean_sd_cov_chr,# %>%
                 #filter(!(chr_n %in% c(1, 11))),
               mapping = aes(x = label, y = mean_anc_corr),
               fun.data = "mean_cl_normal", 
               geom = "errorbar",
               color = "black",
               width = .3) +
  scale_color_viridis_d(name = "Chromosome", labels = 1:16, option = "D") +
  scale_shape_manual(name = "Chromosome",
                     labels = 1:16,
                     values = rep(c(1:3,19), 10)[1:16]) +
  guides(color=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Ancestry covariances by chromosome")
ggsave(paste0("plots/mean_k_standardized_cov_by_groups_and_chr.png"), 
       height = 4, width = 6, 
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures/mean_k_standardized_cov_by_groups_and_chr.pdf"), 
       height = 4, width = 6, 
       units = "in", device = "pdf")
# t-tests
with(filter(mean_sd_cov_chr, type %in%  c("CA_ARS", "ARN_ARS")), t.test(mean_anc_corr ~ type), alternative = "two.sided")
with(filter(mean_sd_cov_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              !(chr_n %in% c(1, 11))), t.test(mean_anc_corr ~ type), alternative = "two.sided")
with(filter(mean_sd_cov_chr, type %in%  c("CA_ARS", "ARN_ARS")) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))
with(filter(mean_sd_cov_chr, type %in%  c("CA_ARS", "ARN_ARS") &
              !(chr_n %in% c(1, 11))) %>% arrange(chr_n), 
     t.test(mean_anc_corr ~ type, paired = T))
     
