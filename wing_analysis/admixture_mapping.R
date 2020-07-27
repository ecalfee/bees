# admixture mapping for wing traits, starting with length:
library(dplyr)
library(tidyr)
library(ggplot2)
library(devtools)
#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("SNPRelate")
#BiocManager::install("SeqArray")
#devtools::install_github("kegrinde/STEAM")
library(STEAM)

# load bee wing length data and model residuals for fit with latitude from plot_wing_lengths.R
load("results/wing_fits.RData") #d.wings, m_wing, m_wing_SA

# get individual local ancestry data for bees with wing data:
# these are A ancestry counts (0,1,2) based on MAP values of the posterior
A_wings_wide_counts <- do.call(cbind, 
                               lapply(d.wings$Bee_ID, function(id) 
                                 read.table(paste0("../local_ancestry/results/ancestry_hmm/combined_sept19/posterior/anc_count/", 
                                                   id, ".A.count"), stringsAsFactors = F)))
colnames(A_wings_wide_counts) <- d.wings$Bee_ID
A_wings_wide_counts <- as.matrix(A_wings_wide_counts)
save(A_wings_wide_counts, file = "results/A_wings_wide_counts.RData")
#load("results/A_wings_wide_counts.RData")


# admixture mapping of residuals after fitting prediction of wing length from global ancestry
fits_counts <- t(sapply(1:nrow(A_wings_wide_counts), function(i){
  m <- lm(d.wings$model_residual_cm ~ A_wings_wide_counts[i, ])
  cf <- summary(m)$coefficients
  results <- c(cf[2,], cf[1,])
  names(results) <- c("b", "se", "t.value", "p.value",
                      "intercept", "intercept_se", 
                      "intercept_t.value", "intercept_p.value")
  return(results)
}))

save(fits_counts, file = "results/ancestry_mapping_wing_length_counts.RData")
#summary(fits_counts[ , "p.value"])
#hist(fits_counts[ , "p.value"])
#load("results/ancestry_mapping_wing_length_counts.RData")


# get chromosome and sites data
chr_lengths <- cbind(read.table("../data/honeybee_genome/chr.names", stringsAsFactors = F),
                     read.table("../data/honeybee_genome/chr.lengths", stringsAsFactors = F)) %>%
  data.table::setnames(c("chr", "scaffold", "chr_length")) %>%
  mutate(chr_n = 1:16) %>%
  mutate(chr_end = cumsum(chr_length)) %>%
  mutate(chr_start = chr_end - chr_length) %>%
  mutate(chr_mid = (chr_start + chr_end)/2)

# SNP sites where ancestry was called
# with physical and cM position based on recombination map
# sites_rpos created in K_by_recomb_rate.R
load("../local_ancestry/results/sites_rpos.RData") # sites_rpos


sites_map <- sites_rpos[ , c("pos_cM", "chr_n")] %>%
   rename(cM = pos_cM, chr = chr_n)

# use R package from Grinde 2019 to get threshold sig. p-value for multiple testing
pval_thresh <- STEAM::get_thresh_analytic(g = 30,
                                          type = "pval",
                                          map = sites_map,
                                          alpha = 0.05)

# range of thresholds for a range of admixture generations g
gs <- c(22, 47.65, 60.4, 150) # min, median, mean, max
pval_threshs <- sapply(gs, function(g) STEAM::get_thresh_analytic(g = g,
                                                                  type = "pval",
                                                                  map = sites_map,
                                                                  alpha = 0.05))
#pval_threshs
#-log10(pval_threshs)

# plot admixture mapping results with
# threshold based on 47.6 gen. admixture (median across pops)
cbind(sites_rpos, data.frame(fits_counts, stringsAsFactors = F)) %>%
  mutate(even_chr = (chr_n %% 2 == 1)) %>%
  ggplot(., aes(x = cum_pos, y = -log10(p.value), color = even_chr)) +
  scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  geom_point(cex = .2) +
  theme_classic() +
  xlab("Chromosome") +
  ylab(expression('-log'[10]*'('*italic('P')*')')) +
  geom_hline(yintercept = -log10(pval_threshs[[2]]), 
             color = "red", lty = "dashed") +
  scale_color_manual(values = brewer.pal(4,"Paired")) +
  theme(legend.position = "none")
ggsave("plots/admixture_mapping_wing_length.png", 
       height = 3, width = 5.2)
ggsave("../../bee_manuscript/figures/admixture_mapping_wing_length.png", 
       height = 3, width = 5.2, dpi = 600)
ggsave("../../bee_manuscript/figures_supp/admixture_mapping_wing_length.tiff", 
       height = 3, width = 5.2, dpi = 600)
