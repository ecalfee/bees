# plots PCA results
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette
library(stringr)
source("../colors.R") # for color palette

# get metadata for individuals included in analysis
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")


prefix <- "combined_sept19"

# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table(paste0("../bee_samples_listed/", prefix, ".list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")

# get coverage estimates
coverage <- read.table(paste0("../geno_lik_and_SNPs/results/", prefix, 
                              "/coverage/mean_ind_coverage.chr.random_pos_1000.txt"),
                       header = T, stringsAsFactors = F, sep = "\t")

# join all data together
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") %>%
  dplyr::left_join(., coverage, by = "Bee_ID")

# pruned every nth snp
n <- 250

# get PCA data
cov_dir <- "results/PCA"
cov_file <- paste0(prefix, "_chr_prunedBy", n, ".cov")
cov_data <- read.table(file.path(cov_dir, cov_file),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
m = 10 # make small dataframe w/ only first m eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:m])
colnames(pca_small) = paste0("PC", 1:m)

#pca_princomp <- princomp(covmat = cov_data)

# join bams and firstPCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(bees, pca_small)  %>%
  arrange(., population) %>%
  arrange(., strain) %>%
  mutate(., label = ifelse(group %in% c("A", "C", "M"), 
                           group,
                           ifelse(group == "AR_2018", "S. America", "N. America")))


# rounded eigen values
PC_var_explained = 100*pca$values/sum(pca$values)


# plot first PC's
# PC1 and 2, all no filter
p12 = d %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", round(PC_var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(PC_var_explained[2], 2), "%)")) +
  theme_classic()
plot(p12 + geom_point(aes(color = strain, size = est_coverage), alpha = .5) +
       ggtitle(paste0("PCA ", prefix, " every ", n, "th snp")))
ggsave(paste0("plots/PCA_12_", prefix, ".png"), 
       plot = p12 + geom_point(aes(color = strain), alpha = .5), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
# by group, not strain
pca_for_manuscript <- p12 + 
  geom_point(aes(color = label), alpha = .5, size = 2) +
  scale_color_manual(values = c(col_ACM, col_NA_SA_both),
                     name = NULL) # + ggtitle("PCA")

plot(pca_for_manuscript)
ggsave(paste0("plots/PCA_12_", prefix, "_byGroup.png"), 
       plot = pca_for_manuscript, 
       device = "png", 
       width = 5.2, height = 4, 
       units = "in",
       dpi = 600)
ggsave(paste0("../../bee_manuscript/figures/PCA_12_", prefix, "_byGroup.png"), 
       plot = pca_for_manuscript, 
       device = "png", 
       width = 5.2, height = 4, 
       units = "in",
       dpi = 600)
ggsave(paste0("../../bee_manuscript/figures_supp/PCA_12_", prefix, "_byGroup.tiff"), 
       plot = pca_for_manuscript, 
       device = "tiff", 
       width = 5.2, height = 4, 
       units = "in",
       dpi = 600)



# look at reference bees more carefully:
p12_ref = d %>%
  filter(d$group %in% c("A", "C", "M", "O")) %>%
  ggplot(., aes(PC1, PC2)) + 
  geom_point(aes(shape = source, color = geographic_location_short)) +
  xlab(paste0("PC1 (", round(PC_var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(PC_var_explained[2], 2), "%)")) +
  ggtitle(paste0("PCA ref ", prefix, " every ", n, "th snp"))
plot(p12_ref)
ggsave(paste0("plots/PCA_REF_12_", prefix, "_byLocationSource.png"), 
       plot = p12_ref, 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
# 3rd and 4th PCs ref bees:
p34_ref = d %>%
  filter(d$group %in% c("A", "C", "M", "O")) %>%
  ggplot(., aes(PC3, PC4)) + 
  geom_point(aes(shape = source, color = geographic_location_short)) +
  xlab(paste0("PC3 (", round(PC_var_explained[3], 2), "%)")) +
  ylab(paste0("PC4 (", round(PC_var_explained[4], 2), "%)")) +
  ggtitle(paste0("PCA ref ", prefix, " every ", n, "th snp"))
plot(p34_ref)
ggsave(paste0("plots/PCA_REF_34_", prefix, "_byLocationSource.png"), 
       plot = p34_ref, 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

p34 = d %>%
  ggplot(., aes(PC3, PC4)) + 
  geom_point(aes(shape = source, color = geographic_location_short)) +
  xlab(paste0("PC3 (", round(PC_var_explained[3], 2), "%)")) +
  ylab(paste0("PC4 (", round(PC_var_explained[4], 2), "%)")) +
  ggtitle(paste0("PCA ref ", prefix, " every ", n, "th snp"))
plot(p34)
ggsave(paste0("plots/PCA_34_", prefix, "_byLocationSource.png"), 
       plot = p34, 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

