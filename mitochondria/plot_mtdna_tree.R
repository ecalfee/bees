# script for creating a haplotype tree and plotting it for my mtDNA haplotypes
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
source("../colors.R") # for color palette

# read in sample IDs:
ids <- read.table("../bee_samples_listed/combined_sept19.list", 
                  header = F, stringsAsFactors = F)$V1

# metadata
load("../local_ancestry/results/meta.RData")
meta <- left_join(data.frame(Bee_ID = ids, stringsAsFactors = F), 
                  read.table("../bee_samples_listed/all.meta", sep = "\t", header = T, stringsAsFactors = F),
                  by = "Bee_ID")


# read in haplotypes from ANGSD
mtdna <- read.table("results/haplotypes_mtdna.haplo.gz", stringsAsFactors = F, header = T, na.strings = "N") %>%
  data.table::setnames(c("chr", "pos", "major", ids))
haplos <- mtdna[ , ids]
mtdna_snps <- mtdna[ , c("chr", "pos", "major")]
apply(haplos, 1, table)


# what depth of coverage do I have, est across 1000 random mtdna sites?
cov_mtdna <- read.table("results/coverage_random_pos_1000_mtdna.txt") %>%
  data.table::setnames(c("chr", "start", "end", "V4", "V5", "strand", ids))
meta$cov_mtdna <- apply(cov_mtdna[ , ids], 2, mean)
cov_mtdna$mean <- apply(cov_mtdna[ , ids], 1, mean)
cov_mtdna %>%
  ggplot(., aes(x = start, y = mean)) +
  geom_point()
summary(cov_mtdna$mean)
group_names = data.frame(group = c("A", "C", "M", "AR_2018", "CA_2018", "N_CA", "S_CA"),
                         group_name = c("A", "C", "M", "Argentina 2018", "California 2018", "N. California 2014", "S. California 2014"),
                         stringsAsFactors = F) %>%
  mutate(group_name = factor(group_name, ordered = T, levels = .$group_name))
p_cov_line <- cov_mtdna %>%
  pivot_longer(cols = ids, names_to = "Bee_ID", values_to = "coverage") %>%
  left_join(., meta, by = "Bee_ID") %>%
  left_join(., group_names, by = "group") %>%
  mutate(source = ifelse(source == "Ramirez", "Cridland", source)) %>%
  ggplot(., aes(x = end, y = coverage, 
                group = group_name, 
                color = group_name)) +
  geom_point(size = 0.1) +
  facet_grid(source~., scales = "free_y") +
  theme_light() +
  labs(color = "", x = "Position (mtDNA)", y = "Total coverage")
p_cov_smooth <- cov_mtdna %>%
  pivot_longer(cols = ids, names_to = "Bee_ID", values_to = "coverage") %>%
  left_join(., meta, by = "Bee_ID") %>%
  left_join(., group_names, by = "group") %>%
  mutate(source = ifelse(source == "Ramirez", "Cridland", source)) %>%
  ggplot(., aes(x = end, y = coverage, 
                group = group_name, 
                color = group_name)) +
  geom_smooth(se = FALSE) +
  facet_grid(source~., scales = "fixed") +
  #facet_grid(source~., scales = "free_y") +
  theme_light() +
  labs(color = "", x = "Position (mtDNA)", y = "Total coverage")
p_cov_smooth
ggsave("plots/mtdna_coverage_est_lines.png", 
       plot = p_cov_line,
       device = "png", 
       width = 7.5, height = 4, units = "in", dpi = 600)
ggsave("plots/mtdna_coverage_est_smooth.png", 
       plot = p_cov_smooth,
       device = "png", 
       width = 7.5, height = 4, units = "in", dpi = 600)

# look at data and filter for quality of SNPs - how many are tri-allelic sites?
#sum(apply(haplos, 1, function(x) sum(!is.na(unique(x)))) >4)/nrow(haplos) # about 6%, filter out:
#haplos_no3 <- haplos[apply(haplos, 1, function(x) sum(!is.na(unique(x)))) <= 4, ] %>%
#  dplyr::select(., - c("chr", "pos", "major")) # keep only genotype columns
#colnames(haplos_no3) <- ids
#haplos_no3.snps <- haplos[apply(haplos, 1, function(x) sum(!is.na(unique(x)))) <= 4, ] %>%
#  dplyr::select(., c("chr", "pos", "major")) # keep only snp position columns
# how much information do I have about haplotypes?
#summary(apply(haplos_no3, 2, function(x) sum(!is.na(x)))) # not that much missing data

# how many tri-allelic sites? just 3. only 1 is convincing though..
mtdna_snps$triallelic <- apply(haplos, 1, function(x) sum(!is.na(unique(x))) > 2)
table(mtdna_snps$triallelic)
# identify major (most common) and minor (2nd most common) alleles
mtdna_snps$allele_A <- apply(haplos, 1, function(x) names(sort(table(x), decreasing=TRUE))[1])
mtdna_snps$allele_a <- apply(haplos, 1, function(x) names(sort(table(x), decreasing=TRUE))[2])
table(mtdna_snps$allele_A == mtdna_snps$major) # always true, great

# major allele = 1, minor allele = 0, any other or no data = NA
counts <- t(sapply(1:nrow(haplos), function(i) sapply(haplos[i, ], function(x)
                                                    ifelse(x == mtdna_snps$allele_A[i], 1, 
                                                           ifelse(x == mtdna_snps$allele_a[i], 0,
                                                                  NA))))) %>%
  as.data.frame(.)
# get allele freqs for ACM groups -- any fixed differences?
heatmap(as.matrix(counts))
heatmap(as.matrix(counts[ , meta$Bee_ID[meta$population %in% c("A", "C", "M")]]))
freqs <- bind_cols(mtdna_snps, counts) %>%
  pivot_longer(., cols = all_of(ids), names_to = "Bee_ID", values_to = "geno") %>%
  left_join(., meta, by = "Bee_ID") %>%
  group_by(chr, pos, population) %>%
  summarise(freq = mean(geno, na.rm = T),
            n = sum(!is.na(geno)))
freqs %>%
  filter(population %in% c("A", "C", "M")) %>%
  ggplot(., aes(x = pos, y = freq, size = n, color = population)) +
  geom_point()
acm_freqs <- freqs %>% # estimated freq of minor allele at snp?
  dplyr::select(-n) %>%
  filter(population %in% c("A", "C", "M")) %>%
  pivot_wider(., names_from = "population", values_from = "freq")
acm_ns <- freqs %>% # how many included samples?
  dplyr::select(-freq) %>%
  filter(population %in% c("A", "C", "M")) %>%
  pivot_wider(., names_from = "population", values_from = "n")
acm_freqs %>%
  filter(abs(A-C) > 0.9)
acm_freqs %>%
  filter(abs(A-M) > 0.9)
#acm_freqs %>%
#  filter(M %in% c(0,1)) %>%
#  View()
aims_mtdna <- acm_freqs %>%
  left_join(., mtdna_snps, by = c("chr", "pos")) %>%
  filter(abs(A-M) > 0.8 & abs(A-C) > 0.8) %>%
  left_join(., acm_ns, by = c("chr", "pos"), suffix = c("p", "n"))

# global ancestry:
load("../global_ancestry/results/NGSAdmix/ACM_K3_combined_sept19_chr_prunedBy250.rData")

# not perfect markers, but plot clines at these 2 snps:
p_mtdna <- filter(freqs, pos %in% aims_mtdna$pos) %>%
  mutate(marker = paste(chr, pos, sep = ":")) %>%
  left_join(., meta.pop, by = "population") %>%
  rename(continent = zone) %>%
  filter(!(population %in% c("A", "C", "M"))) %>%
  ggplot(.) +
  geom_point(aes(x = abs(lat), y = 1-freq, alpha = n, 
                 color = continent)) +
  facet_grid(marker~continent) +
  #theme_classic() +
  theme_light() +
  ylab("Allele frequency") +
  xlab("Degrees latitude from the equator") +
  geom_point(data = d_admix_ACM_labelled %>%
               filter(!is.na(continent)) %>%
               group_by(population, continent) %>%
               summarise(A = mean(A), lat = mean(lat)), 
             aes(x = abs(lat),
                 y = A), 
             color = "black",
             shape = 1) +
  geom_hline(data = aims_mtdna %>%
              dplyr::select(chr, pos) %>%
              left_join(., acm_freqs) %>%
              mutate(marker = paste(chr, pos, sep = ":")),
            aes(yintercept = 1 - A), color = col_ACM["A"],
            linetype = "dashed") +
  scale_color_manual(values = col_NA_SA_both) +
  labs(alpha = "Sample size", color = "Continent")
p_mtdna
ggsave("plots/mtdna_snp_clines.png", 
       plot = p_mtdna,
       device = "png", 
       width = 5.2, height = 4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/mtdna_snp_clines.png", 
       plot = p_mtdna,
       device = "png", 
       width = 5.2, height = 4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/mtdna_snp_clines.tiff", 
       plot = p_mtdna,
       device = "tiff", 
       width = 5.2, height = 4, units = "in", dpi = 600)

# convert to ape recognized format:
a <- as.DNAbin(t(as.matrix(haplos)), fill.with.gaps = T)
print(a) 
# just reference bees:
a_ACM <- as.DNAbin(t(as.matrix(haplos[ , meta$Bee_ID[meta$population %in% c("A", "C", "M")]])))

# alternatively, create a binary matrix where rows are loci and columns are individuals. 
# 1 means major allele 0 means minor allele
#a.01 <- haplos_no3
#for (i in 1:nrow(haplos_no3)){
#  for (j in 1:ncol(haplos_no3)){
#    a.01[i, j] <- as.integer(haplos_no3[i, j] != haplos_no3.snps[i, "major"]) # 0 for major allele, 1 for minor, NA for missing data
#  }
#}

# use NJ tree from ape package
# get distance under a particular model of molecular evolution
# Q: is it ok that I only include variant sites? I may need the whole fasta sequence or my tree is just wonky because of low coverage data
a.dist <- dist.dna(a, model = "K80", variance = FALSE,
         gamma = FALSE, pairwise.deletion = FALSE,
         base.freq = NULL, as.matrix = FALSE)
# plot neighbor-joining tree
a.nj <- njs(a.dist)
plot(a.nj, tip.color = rainbow(7)[as.factor(meta$group)])
print(a.nj)

# just ref bees:
a_ACM.dist <- dist.dna(a_ACM, model = "K80", variance = FALSE,
                   gamma = FALSE, pairwise.deletion = FALSE,
                   base.freq = NULL, as.matrix = FALSE)
# plot neighbor-joining tree
a_ACM.nj <- njs(a_ACM.dist)
plot(a_ACM.nj, 
     tip.color = 
       rainbow(3)[as.factor(meta$population[meta$population %in% c("A", "C", "M")])])



