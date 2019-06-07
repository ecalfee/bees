# plots PCA results
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette
library(stringr)

# get metadata for individuals included in analysis
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")


#prefix <- "pass1" # NOTE: bee ap50 is a duplicate entry in pass1, but not any others below:
#prefix <- "pass1_and_kohn"
#prefix <- "pass1_plus_kohn_and_wallberg"
#prefix <- "pass1_plus_kohn_and_wallberg_noJordanO"
prefix <- "CA_AR_MX_harpur_sheppard_kohn_wallberg"

# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table(paste0("../bee_samples_listed/", prefix, ".list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")

# get coverage estimates
coverage <- read.table(paste0("../geno_lik_and_SNPs/results/", prefix, "/coverage/est_coverage_by_reads_Q30.txt"),
                       header = T, stringsAsFactors = F, sep = "\t")

# join all data together
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") %>%
  dplyr::left_join(., coverage, by = "Bee_ID")

# pruned every nth snp
#n <- 1000
n <- 250

# get PCA data
cov_dir <- "results/PCA"
# cov_file <- "ordered_scaffolds_prunedBy", n, ".cov" # just pass1
cov_file <- paste0("ordered_scaffolds_", prefix, "_prunedBy", n, ".cov")
cov_data <- read.table(file.path(cov_dir, cov_file),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
m = 10 # make small dataframe w/ only first m eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:m])
colnames(pca_small) = paste0("PC", 1:m)
  

# join bams and firstPCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(bees, pca_small)  %>%
  arrange(., population) %>%
  arrange(., strain)


# rounded eigen values
PC_var_explained = round(pca$values, 2)


# plot first PC's
# PC1 and 2, all no filter
p12 = d %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle(paste0("PCA ", prefix, " every ", n, "th snp"))
plot(p12 + geom_point(aes(color = strain, size = est_coverage), alpha = .5))
ggsave(paste0("plots/PCA_12_", prefix, ".png"), 
       plot = p12 + geom_point(aes(color = strain), alpha = .5), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
# by group, not strain
plot(p12 + geom_point(aes(color = group), alpha = .5))
ggsave(paste0("plots/PCA_12_", prefix, "_byGroup.png"), 
       plot = p12 + geom_point(aes(color = group), alpha = .5), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# look at reference bees more carefully:
p12_ref = d %>%
  filter(d$group %in% c("A", "C", "M", "O")) %>%
  ggplot(., aes(PC1, PC2)) + 
  geom_point(aes(shape = source, color = geographic_location_short)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
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
  xlab(paste0("PC3 (", PC_var_explained[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained[4], "%)")) +
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
  xlab(paste0("PC3 (", PC_var_explained[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained[4], "%)")) +
  ggtitle(paste0("PCA ref ", prefix, " every ", n, "th snp"))
plot(p34)
ggsave(paste0("plots/PCA_34_", prefix, "_byLocationSource.png"), 
       plot = p34, 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# Use vars() to supply faceting variables:
plot(p12 + geom_point(aes(color = population), alpha = .5) + facet_wrap(~group))
ggsave(paste0("plots/PCA_", prefix, "facet_by_population.png"), 
       plot = p12 + geom_point(aes(color = population), alpha = .5) + facet_wrap(~source), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# just look at old bees
p12_old = d %>%
  filter(., source != "Calfee") %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle(paste0("PCA ", prefix, " every ", n, "th snp"))
plot(p12_old + geom_point(aes(color = population), alpha = .5))
# PC3 is avalon and PC4 is driven by Riverside 1999 b/c of sample duplicate
p34_old = d %>%
  filter(., source != "Calfee") %>%
  ggplot(., aes(PC3, PC4)) + 
  xlab(paste0("PC3 (", PC_var_explained[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained[4], "%)")) +
  ggtitle(paste0("PC3 and PC4 (all ind's); every ", n, "th snp"))
plot(p34_old + geom_point(aes(color = population), alpha = .5))
ggsave("plots/PCA_3_driven_by_Avalon.png", 
       plot = p34_old+ geom_point(aes(color = population), alpha = .5),
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# no clear separation of Argentina vs. California in any of the first 10 PCs
p34 = d %>%
  ggplot(., aes(PC3, PC4)) + 
  xlab(paste0("PC3 (", PC_var_explained[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained[4], "%)")) +
  ggtitle(paste0("PCA 34 ", prefix, " every ", n, "th snp"))
plot(p34 + geom_point(aes(color = group), alpha = .5))
table(d$population)
plot(p34 + facet_wrap(~group) + geom_point(aes(color = population), alpha = .5))
ggsave(paste0("plots/PCA_34_", prefix, "by_population.png"), 
       plot = p34 + facet_wrap(~group) + geom_point(aes(color = population), alpha = .5),
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
plot(p12 + facet_wrap(~group) + geom_point(aes(color = population), alpha = .5))
ggsave(paste0("plots/PCA_12_", prefix, "by_population.png"), 
       plot = p12 + facet_wrap(~group) + geom_point(aes(color = population), alpha = .5),
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# PC3 and PC6 to some extent separate CA from AR hybrid zone samples
# PC 5,6
p56 = d %>%
  ggplot(., aes(PC5, PC6)) + 
  xlab(paste0("PC5 (", PC_var_explained[5], "%)")) +
  ylab(paste0("PC6 (", PC_var_explained[6], "%)")) +
  ggtitle(paste0("PC5 and PC6 (all ind's); every ", n, "th snp"))
plot(p56 + geom_point(aes(color = group), alpha = .5))
ggsave(paste0("plots/PCA_56_", prefix, ".png"), 
       plot = p56 + geom_point(aes(color = group), alpha = .5),
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC 7,8
p78_new = d %>%
  ggplot(., aes(PC7, PC8)) + 
  xlab(paste0("PC7 (", PC_var_explained[7], "%)")) +
  ylab(paste0("PC8 (", PC_var_explained[8], "%)")) +
  ggtitle(paste0("PC7 and PC8 (all ind's); every ", n, "th snp"))
plot(p78_new + geom_point(aes(color = group), alpha = .5))
# PC 9,10
p910_new = d %>%
  ggplot(., aes(PC9, PC10)) + 
  xlab(paste0("PC9 (", PC_var_explained[9], "%)")) +
  ylab(paste0("PC10 (", PC_var_explained[10], "%)")) +
  ggtitle(paste0("PC9 and PC10 (all ind's); every ", n, "th snp"))
plot(p910 + geom_point(aes(color = group), alpha = .5))

# do I see population structure? Not obviously
p12_AR = d %>%
  filter(., group == "AR_2018") %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle(paste0("PC1 and PC2 (AR_2018); every ", n , "th snp"))
plot(p12_AR + geom_point(aes(color = population), alpha = .5))
# even less in CA
p12_CA = d %>%
  filter(., group == "CA_2018") %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle(paste0("PC1 and PC2 (CA_2018); every ", n , "th snp"))
plot(p12_CA + geom_point(aes(color = population), alpha = .5))
