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
n = 10 # make small dataframe w/ only first n eigenvectors
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
calc_PC_dist <- function(sample_PC_pos, ancestry_PC_pos = anc_pos, eigen_values = eigen_small){
 apply(ancestry_PC_pos, 2, function(p) sqrt(sum(((sample_PC_pos-p)*eigen_values)^2)))
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
ggplot(dist) +
  geom_point(aes(x = proportion, y = PC_dist, color = refAncestry)) +
  facet_wrap(~ngsAdmixAncestry)
dist %>%
  filter(refAncestry == "A" & ngsAdmixAncestry == "anc1") %>%
  ggplot(.) +
  geom_point(aes(x = proportion, y = PC_dist, color = group))
dist %>%
  filter(refAncestry == "M" & ngsAdmixAncestry == "anc2") %>%
  ggplot(.) +
  geom_point(aes(x = proportion, y = PC_dist, color = group))
dist %>%
  filter(refAncestry == "C" & ngsAdmixAncestry == "anc3") %>%
  ggplot(.) +
  geom_point(aes(x = proportion, y = PC_dist, color = group))
