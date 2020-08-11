library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
library(gridExtra)
library(boot)
source("../colors.R") # get color palette
source("het_fst_functions.R") # get pi and fst functions
source("calc_frac_snps.R") # calculates snp density (to scale pi)
# this script takes the population allele frequencies
# from combine_pop_allele_freqs.R
# and calculates pi and Fst between all populations

#------------------ get data ----------------------#
pops <- read.table("../bee_samples_listed/byPop/combined_sept19_pops.list",
                   stringsAsFactors = F)$V1
ACM <- c("A", "C", "M")
ACM_pops <- c(ACM, pops)
ancestries <- ACM
ancestries_combined <- c(ancestries, "combined")
ancestries_Combined <- c(ancestries, "Combined")

# get bee population meta data
load("../local_ancestry/results/meta.RData")


# load new data?
load_new_data = T

# read in allele freq and heterozygosity data:
if (load_new_data) {
  freqs_pops <- do.call(rbind, 
                   lapply(1:16, function(i) read.table(paste0("results/combined_sept19/combined/allele_freq/Group", i, 
                                                              "/pops_included_every_SNP.freqs.txt"),
                                                              header = T, stringsAsFactors = F)))
  freqs_ACM <- do.call(rbind, 
                   lapply(1:16, function(i) read.table(paste0("results/combined_sept19/combined/allele_freq/Group", i, 
                                                              "/ACM_every_SNP.freqs.txt"),
                                                       header = T, stringsAsFactors = F)))
  freqs <- cbind(freqs_ACM, freqs_pops)
  rm(freqs_ACM, freqs_pops)

  hets <- 2*freqs*(1-freqs)
  het_mean <- apply(hets, 2, function(x) mean(x, na.rm = T))*frac_snps
  rm(hets)
  
  # how many alleles sampled at a site?
  ns_pops <- do.call(rbind, 
                        lapply(1:16, function(i) read.table(paste0("results/combined_sept19/combined/allele_freq/Group", i, 
                                                                   "/pops_included_every_SNP.nInd"),
                                                                   #"/pops_included_every", n_snp, "th_SNP.nInd"),
                                                            header = T, stringsAsFactors = F)))*2 # times 2 for diploid
  ns_ACM <- do.call(rbind, 
                       lapply(1:16, function(i) read.table(paste0("results/combined_sept19/combined/allele_freq/Group", i, 
                                                                  "/ACM_every_SNP.nInd"),
                                                                  #"/ACM_every", n_snp, "th_SNP.nInd"),
                                                           header = T, stringsAsFactors = F)))*2
  ns <- cbind(ns_ACM, ns_pops) # x2 because diploid
  rm(ns_ACM, ns_pops)
  
  # mean across A/C/M
  ACM_het <- data.frame(ancestry = ancestries,
                        combined = het_mean[ACM],
                        ref_pop = T)
  
  
  # heterozygosity corrected for small sample size:
  # by default also drop any snps with fewer than 2 individuals with data (<= 2 alleles observed, i.e. 1 ind.)
  hets_small_sample <- do.call(cbind, lapply(1:ncol(freqs), 
                                             function(i) het_small_sample_correction(p = freqs[, i],
                                                                                     n = ns[, i])))
  colnames(hets_small_sample) <- ACM_pops
  save(list = "hets_small_sample", file = "results/hets_small_sample_combined.RData")
  het_small_sample_mean <- apply(hets_small_sample, 2, function(x) mean(x, na.rm = T))*frac_snps
  rm(hets_small_sample)
  
  ACM_het_small_sample <- data.frame(ancestry = ACM,
                                     combined = het_small_sample_mean[ACM],
                                     ref_pop = T)
  
  save(file = "results/freqs_combined.RData", 
       list = c("freqs", "ns", "het_mean", "ACM_het", 
                "het_small_sample_mean", "ACM_het_small_sample"))
} else {
  load("results/freqs_combined.RData")
}

# freqs and heterozygosity within AA CC MM ancestry
if (load_new_data){
  freqs_by_ancestry <- lapply(ancestries, function(a) 
    cbind(freqs[ , ACM], 
          do.call(rbind, 
                  lapply(1:16, function(i) read.table(paste0("results/combined_sept19/", a, "/allele_freq/Group", i, 
                                                             "/pops_included_every_SNP.freqs.txt"),
                                                      header = T, stringsAsFactors = F)))))
  hets_by_ancestry <- lapply(freqs_by_ancestry, function(f) 2*f*(1-f))
  het_mean_by_ancestry <- lapply(hets_by_ancestry, function(h)
    apply(h, 2, function(x) mean(x, na.rm = T)*frac_snps))
  rm(hets_by_ancestry)
  
  # how many alleles were sampled at each site?
  ns_by_ancestry <- lapply(ancestries, function(a) 
    cbind(ns[ , ACM], 
          do.call(rbind, 
                  lapply(1:16, function(i) read.table(paste0("results/combined_sept19/", a, "/allele_freq/Group", i,
                                                             "/pops_included_every_SNP.nInd"),
                                                      header = T, stringsAsFactors = F)*2)))) # x2 for diploid
  
d_het <- data.frame(population = ACM_pops,
                      A = het_mean_by_ancestry[[1]],
                      C = het_mean_by_ancestry[[2]],
                      M = het_mean_by_ancestry[[3]],
                      combined = het_mean,
                      ref_pop = c(rep(T, 3), rep(F, length(pops)))) %>%
    left_join(., meta.pop, by = "population")

    # use small sample size correction and filter sites with < 2 individuals with data
  hets_small_sample_by_ancestry <- lapply(1:3, function(a) do.call(cbind, 
                                                                   lapply(1:ncol(freqs_by_ancestry[[a]]), 
                                                                          function(i) het_small_sample_correction(p = freqs_by_ancestry[[a]][ , i], 
                                                                                                                  n = ns_by_ancestry[[a]][ , i])))) 
  # add columns/population names, then save
  for (i in 1:3){
    colnames(hets_small_sample_by_ancestry[[i]]) <- colnames(freqs_by_ancestry[[i]])
  }
  save(list = "hets_small_sample_by_ancestry",
       file = "results/hets_small_sample_by_ancestry.RData")
  het_small_sample_mean_by_ancestry <- lapply(hets_small_sample_by_ancestry, function(h) 
    apply(h, 2, function(x) mean(x, na.rm = T)*frac_snps))
  rm(hets_small_sample_by_ancestry)
  d_het_small_sample <- data.frame(population = ACM_pops,
                                   A = het_small_sample_mean_by_ancestry[[1]],
                                   C = het_small_sample_mean_by_ancestry[[2]],
                                   M = het_small_sample_mean_by_ancestry[[3]],
                                   Combined = het_small_sample_mean,
                                   ref_pop = c(rep(T, 3), rep(F, length(pops)))) %>%
    left_join(., meta.pop, by = "population") %>%
    mutate(year = factor(year, levels = c("2018", "2014", "2002", "1999")))

  save(file = "results/freqs_within_ancestry.RData", 
       list = c("freqs_by_ancestry", "ns_by_ancestry", "d_het", "het_mean_by_ancestry",
                "d_het_small_sample", "het_small_sample_mean_by_ancestry"))
} else {
  load("results/freqs_within_ancestry.RData")
}

# load bootstrap results for pi within ancestry:
bootstrap_pi <- do.call(rbind, lapply(c("Combined", ACM), function(a) 
  read.table(paste0("results/block_bootstrap_pi_within_", a, "_boots.txt"),
             header = T, sep = "\t", stringsAsFactors = F) %>%
    mutate(ancestry = a))) %>%
  # scale heterozygosity estimates per snp to het per bp
  mutate(estimate = frac_snps*estimate,
         lower = frac_snps*lower,
         upper = frac_snps*upper) %>%
  mutate(ancestry = factor(ancestry, levels = c("Combined", ACM), ordered = T)) %>% # for better plot order
  mutate(exclude = (chr < 15 | n_bins < 75))
save(file = "results/bootstrap_pi.RData", list = "bootstrap_pi")
#load("results/bootstrap_pi.RData")



# -----------------------Plot Fig 4 - pi bootstrap -------------------#

# Plot of diversity (pi) across each hybrid zone (w/ bootstrap):
bootstrap_pi %>%
  left_join(., meta.pop, by = c("pop"="population")) %>%
  filter(chr >= 15, n_bins >= 75) %>% # minimum thresholds to include a population (at least 15 chromosomes and 75 bins with data)
  filter(!(pop %in% ACM)) %>%
  ggplot(., aes(x = abs(lat), y = estimate, color = ancestry, 
                shape = factor(year, levels = c("2018", "2014"), ordered = T))) +
  geom_hline(data = filter(bootstrap_pi,
                            pop %in% ACM & ancestry == pop) %>%
               arrange(ancestry) %>%
               head(n = 3),
              aes(yintercept = estimate), color = c(col_ACM, col_ACM), # manually set. double for upper & lower panel 
             alpha = 0.75, show.legend = F) +
  geom_point(size = 1, alpha = .75, aes(color = ancestry)) + #, position = position_dodge2(width = 0.05)) +
  geom_linerange(aes(ymin = lower, ymax = upper), lwd = 0.25) +#, position = position_dodge2(width = 0.05)) +
  
  xlab("Degrees latitude from the equator") +
  ylab(expression(pi)) + # pi symbol
  facet_grid(zone ~ ., scales = "free_x") +
  theme_classic() +
  scale_color_manual(values = c(col_ACM_all)[c("Combined", ACM)], name = "Ancestry") +
  scale_shape_manual(values = c(1, 2)) +
  labs(shape = "Year") + 
  guides(shape = guide_legend(order = 2), 
         color = guide_legend(order = 1,
                              override.aes = list(shape = 15, linetype = "blank"))) # set order of legend so Ancestry type comes first

ggsave("plots/pi_by_latitude.png", device = "png",
       width = 5.2, height = 3)
ggsave("../../bee_manuscript/figures/pi_by_latitude.png", device = "png",
       width = 5.2, height = 3, dpi = 600)
ggsave("../../bee_manuscript/figures_main/pi_by_latitude.tiff", device = "tiff",
       width = 5.2, height = 3, dpi = 600)



# ----------------- Heterozygosity predictions --------------------#
# Compare heterozygosity observed to predicted heterozygosity based on admixture fractions alone
# what do we expect pi to be for these admixed populations?
# first I need admixture data for each population:
# get admixture data
load("../global_ancestry/results/NGSAdmix/ACM_K3_combined_sept19_chr_prunedBy250.rData")

admix.pops <- group_by(d_admix_ACM_labelled, population) %>%
  summarise(A = mean(A), M = mean(M), C = mean(C), n = n())


# predictions based on reference panel allele frequencies
freqs_predicted_ACM <- do.call(cbind, 
                               lapply(pops, function(p) 
                                 as.matrix(freqs[ , ACM]) %*% 
                                   t(as.matrix(admix.pops[admix.pops$population == p, 
                                                          ACM])))) 
colnames(freqs_predicted_ACM) <- pops


hets_small_sample_predicted_ACM <- do.call(cbind, lapply(1:ncol(freqs_predicted_ACM), 
                                                         function(i) het_small_sample_correction(p = freqs_predicted_ACM[, i],
                                                                                                 n = ns[, i])))
colnames(hets_small_sample_predicted_ACM) <- pops

het_small_sample_mean_predicted_ACM <- apply(hets_small_sample_predicted_ACM, 2, function(x) mean(x, na.rm = T))*frac_snps
rm(hets_small_sample_predicted_ACM)


p_predicted_pi <- d_het_small_sample %>%
  filter(!ref_pop) %>%
  left_join(., data.frame(population = pops, 
                          predicted_ref = het_small_sample_mean_predicted_ACM), 
            by = "population") %>%
  rename(observed = Combined) %>%
  tidyr::gather(., "diversity", "pi", c("observed", "predicted_ref")) %>%
  mutate(year = factor(year)) %>%
  ggplot(., aes(x = abs(lat), y = pi, color = diversity, shape = year)) +
  geom_point(alpha= 0.75) +
  xlab("Degrees latitude from the equator") +
  facet_grid(zone ~ ., scales = "free_x") + 
  theme_classic() +
  ylab(expression(pi)) +
  scale_color_manual(values = col_pi_predictions, labels = c("Observed", "Predicted"), name = "Diversity estimate") +
  labs(shape = "Year") +
  guides(shape = guide_legend(order = 2), color = guide_legend(order = 1))


p_predicted_pi_comparison <- d_het_small_sample %>%
  filter(!ref_pop) %>%
  left_join(., data.frame(population = pops, 
                          #predicted_admix = het_small_sample_mean_predicted,
                          predicted_ref = het_small_sample_mean_predicted_ACM), 
            by = "population") %>%
  rename(observed = Combined) %>%
  ggplot(., aes(x = predicted_ref, y = observed, color = zone, shape = population=="Avalon_2014")) +
  geom_point() +
  xlab(expression(paste("Predicted ", pi))) +
  ylab(expression(paste("Observed ", pi))) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  scale_color_manual(values = col_NA_SA_both, name = "Continent") +
  guides(shape = F) +
  coord_fixed()

g_predicted_pi <- grid.arrange(nrow = 2, grobs = list(p_predicted_pi + ggtitle("A"), p_predicted_pi_comparison + ggtitle("B")))
ggsave("plots/pi_observed_and_predicted_from_admixture.png",
       plot = g_predicted_pi,
       width = 5.4, height = 7, units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures/pi_observed_and_predicted_from_admixture.png",
       plot = g_predicted_pi,
       width = 5.4, height = 7, units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/pi_observed_and_predicted_from_admixture.tiff",
       plot = g_predicted_pi,
       width = 5.4, height = 7, units = "in", dpi = 600, device = "tiff")

# linear model for this prediction
d_het_small_sample %>%
  filter(!ref_pop) %>%
  left_join(., data.frame(population = pops, 
                          #predicted_admix = het_small_sample_mean_predicted,
                          predicted_ref = het_small_sample_mean_predicted_ACM), 
            by = "population") %>%
  rename(observed = Combined) %>%
  with(., lm(observed ~ predicted_ref)) %>%
  summary()

# write out mean small sample corrected within ancestry, combined, and predicted heterozygosities:
d_het_small_sample %>%
  left_join(., data.frame(population = pops, 
                          Predicted = het_small_sample_mean_predicted), 
            by = "population") %>%
  write.table(.,
              "results/mean_het_ACM_combined_predicted_small_sample_corrected.txt",
              sep = "\t", col.names = T, row.names = F)

# -------------------- other --------------------------#
# some other summary stats - basic genomewide Fst A-C-M:
# w/out small sample size correction
complete_freqs_ACM = na.omit(freqs[, ACM]) %>%
  mutate(AC = (A+C)/2,
         AM = (A+M)/2,
         CM = (C+M)/2)
frac_snps_complete_ACM <- nrow(complete_freqs_ACM)/genome_size
complete_het_ACM = 2*complete_freqs_ACM*(1-complete_freqs_ACM)
complete_het_ACM_mean = apply(complete_het_ACM, 2, mean)
# calc basic Fst:
1 - mean(complete_het_ACM_mean[c("A", "C")])/complete_het_ACM_mean[c("AC")]
1 - mean(complete_het_ACM_mean[c("A", "M")])/complete_het_ACM_mean[c("AM")]
1 - mean(complete_het_ACM_mean[c("C", "M")])/complete_het_ACM_mean[c("CM")]

# dxy: divergence isnt' that high, it's mostly bottlenecks that increase Fst here
complete_dxy_ACM = complete_freqs_ACM %>%
  mutate(dxyAC = (A*(1-C) + C*(1-A)),
         dxyAM = (A*(1-M) + M*(1-A)),
         dxyCM = (M*(1-C) + C*(1-M))) %>%
  dplyr::select(starts_with("dxy")) %>%
  apply(., 2, mean)
#pi
complete_dxy_ACM*frac_snps_complete_ACM
complete_het_ACM_mean*frac_snps_complete_ACM
# fixed differences are very rare
table(complete_freqs_ACM$A==1 & complete_freqs_ACM$C==0 & complete_freqs_ACM$M ==0)
table(complete_freqs_ACM$A==0 & complete_freqs_ACM$C==1 & complete_freqs_ACM$M ==1)

