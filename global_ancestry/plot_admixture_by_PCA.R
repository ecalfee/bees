library(ggplot2)
library(dplyr)
library(tidyr)
# load 
# get metadata for individuals included in analysis
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")

# get ID's for PCA/NGSadmix data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table("../bee_samples_listed/pass1.list", stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") 
# NOTE: bee ap50 is a duplicate entry I will fix later (!)

# get PCA data
cov_dir <- "results/PCA"
cov_file <- "ordered_scaffolds_prunedBy1000.cov"
cov_data <- read.table(file.path(cov_dir, cov_file),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
n = 2 # make small dataframe w/ only first n eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:n])
colnames(pca_small) = paste0("PC", 1:n)
eigen_small <- pca$values[1:n]

# get NGSadmix results
# starting with pass1 analysis from 1st round of sequencing
K = 3
name = paste0("K", K, "_ordered_scaffolds_prunedBy250")
file = paste0("results/NGSAdmix/", name, ".qopt")
admix <- read.table(file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3")
#colnames(admix) <- c("A", "M", "C")

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(bees, admix, pca_small)  %>%
  arrange(., lat) %>%
  arrange(., source) %>%
  arrange(., group)

# get in PC space the location of M C and A ancestry
ancestries = c("A", "M", "C")
anc_pos <- data.frame(sapply(ancestries, function(a) filter(d, group == a) %>%
  select(starts_with("PC")) %>%
  colMeans(.)))


# function to calculate the PC distance between any sample and the A/C/M groups
calc_PC_dist <- function(sample_PC_pos, ancestry_PC_pos = anc_pos){
 apply(ancestry_PC_pos, 2, function(p) sqrt(sum((sample_PC_pos-p)^2)))
  }
# test
#s <- filter(d, Bee_ID == "AR2512") %>%
#  select(., starts_with("PC"))
#calc_PC_dist(s)

# calculate PC distances for each sample and add to new dataframe 'dist' 
dist <- select(d, starts_with("PC")) %>%
  apply(., 1, function(x) calc_PC_dist(sample_PC_pos = x)) %>%
  t(.) %>%
  cbind(d, .) %>%
  gather(., "refAncestry", "PC_dist", ancestries) %>%
  gather(., "ngsAdmixAncestry", "proportion", c("anc1", "anc2", "anc3"))

# plot
p1 <- ggplot(dist) +
  geom_point(aes(x = proportion, y = PC_dist, color = refAncestry)) +
  facet_wrap(~ngsAdmixAncestry)
plot(p1)
ggsave("plots/PC_vs_NGSadmix_ancestry_all.png", 
       plot = p1, 
       device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)

p2A <- dist %>%
  filter(refAncestry == "A" & ngsAdmixAncestry == "anc1") %>%
  ggplot(.) +
  geom_point(aes(x = proportion, y = PC_dist, color = group)) +
  ggtitle("PC vs. A ancestry")

p2M <- dist %>%
  filter(refAncestry == "M" & ngsAdmixAncestry == "anc2") %>%
  ggplot(.) +
  geom_point(aes(x = proportion, y = PC_dist, color = group)) +
  ggtitle("PC vs. M ancestry")

p2C <- dist %>%
  filter(refAncestry == "C" & ngsAdmixAncestry == "anc3") %>%
  ggplot(.) +
  geom_point(aes(x = proportion, y = PC_dist, color = group)) +
  ggtitle("PC vs. C ancestry")

# ancestry aligns well with positions of groups in PCA space
plot(p2A)
plot(p2M)
plot(p2C)

ggsave("plots/PC_vs_NGSadmix_ancestry_A.png", 
       plot = p2A, 
       device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
ggsave("plots/PC_vs_NGSadmix_ancestry_M.png", 
       plot = p2M, 
       device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
ggsave("plots/PC_vs_NGSadmix_ancestry_C.png", 
       plot = p2C, 
       device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)


# now use ancestry proportions to predict a PC position for each sample
# how well does it do?
# function to calculate predicted PC coordinates given A/C/M admixture proportions
predict_PC_pos <- function(sample_admixture, ancestry_PC_pos = anc_pos){
  as.matrix(ancestry_PC_pos) %*% as.numeric(sample_admixture)
}
# test
#a <- filter(d, Bee_ID == "AR2512") %>%
#  select(., starts_with("anc"))
#predict_PC_pos(a)

# predict positions for all samples:
predict <- select(d, starts_with("anc")) %>%
  apply(., 1, function(x) predict_PC_pos(sample_admixture = x)) %>%
  t(.)
colnames(predict) <- paste0("PC", 1:n)

# create a new database with predicted, rather than actual, PC positions
d2 <- select(d, -starts_with("PC")) %>%
  cbind(., predict)

# make new column for if PC is a prediction or not
d2$predict <- "predicted_PCA"
d$predict <- "true_PCA"

# plot actual and predicted PC positions:
# ancestry predictes PC's really well 
# -- first two PCs are predominantly driven by ancestry
p3 <- rbind(d, d2) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = group)) +
  facet_wrap(~predict) + 
  ggtitle("True PCA vs. prediction from admixture %")
plot(p3) 
ggsave("plots/PC_true_vs_predicted_all.png", 
       plot = p3, 
       device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)


# in euclidean distance, how far off are our predicted points?
d_diff <- select(d, starts_with("PC")) - select(d2, starts_with("PC"))
p4 <- d_diff %>%
  cbind(select(d, -starts_with("PC")), .) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = group)) +
  ggtitle("Diff true PCA vs. prediction from admixture %")
# differences in ancestry-predicted and true PC position
# for admixed individuals are on the same order as
# variation around PC position within reference A/C/M populations
plot(p4)
ggsave("plots/PC_diff_PCA_true_vs_predicted.png", 
       plot = p4, 
       device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)


# calculate the euclidean distance between true and predicted PCs
d_diff_dist <- sqrt(rowMeans(d_diff^2))
p5 <- cbind(d, d_diff_dist) %>%
  ggplot(aes(d_diff_dist)) +
  geom_histogram(aes(fill = group))
plot(p5)
ggsave("plots/PC_diff_hist_true_vs_predicted.png", 
       plot = p5, 
       device = png(), 
       width = 8, height = 8, units = "in",
       dpi = 200)

# I conclude the error seems mostly random, with a few AR_2018 outliers
cbind(d, d_diff_dist) %>%
  filter(group == "AR_2018") %>%
  ggplot(aes(d_diff_dist)) +
  geom_histogram(aes(fill = population))
