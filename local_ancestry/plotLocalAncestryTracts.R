# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)

bees <- read.table("results/SNPs/thin1kb_common3/pass1_2018.ploidy", stringsAsFactors = F, 
                                     header = F, sep = "\t")$V1
prior <- c(0.4, 0.4, 0.2) # prior on global ancestry proportions
colnames(pop_admix_global) = c("popN", "popMaize", "popMex")
#dir_post = "results/ancestry_hmm/thin1kb_common3/pass1_2018_0.4_0.4_0.2/fixed_t_60_30"
#dir_post = "results/ancestry_hmm/thin1kb_common3/pass1_2018_0.05_0.05_0.9"
dir_post = "results/ancestry_hmm/thin1kb_common3/pass1_2018_0.4_0.4_0.2/reverse_order_CMA"
genotypes = c("AA", "AC", "AM", "CC", "CM", "MM")

# get posterior for individuals
getPost = function(id, dir = dir_post){ 
  post <- read.table(paste0(dir, "/", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., CC = X2.0.0) %>%
    rename(., CM = X1.1.0) %>%
    rename(., AC = X1.0.1) %>%
    rename(., MM = X0.2.0) %>%
    rename(., AA = X0.0.2) %>%
    rename(., AM = X0.1.1) %>%
    #rename(., AA = X2.0.0) %>%
    #rename(., AC = X1.1.0) %>%
    #rename(., AM = X1.0.1) %>%
    #rename(., CC = X0.2.0) %>%
    #rename(., MM = X0.0.2) %>%
    #rename(., CM = X0.1.1) %>%
    mutate(., Bee_ID = id)
  post$max_p = apply(post[ , genotypes], 1, max) # which ancestry state has the highest probability?
  small1 <- tidyr::gather(post, "anc", "p", genotypes) %>%
    filter(., p == max_p) %>% # only keep highest posterior prob. ancestry; and filter to one in every nSkip + 1
    arrange(., chrom, position)
  return(small1)
}

# ancestry proportions
mean_anc <- function(post){
  anc1 <- c(1, .5, .5, 0, 0, 0)
  anc2 <- c(0, .5, 0, 1, .5, 0)
  anc3 <- c(0, 0, .5, 0, .5, 1)
  anc_matrix <- cbind(anc1, anc2, anc3)
  anc <- t(apply(post[ , 3:8], 1, function(row) row %*% anc_matrix))
}

# moderate A ancestry
postCA0401 = read.table(paste0(dir_post, "/", "CA0401", ".posterior"), stringsAsFactors = F, header = T)
postSmallCA0401 = getPostSmall(id = "CA0401")

# high A ancestry
postCA0401 = read.table(paste0(dir_post, "/", "CA0401", ".posterior"), stringsAsFactors = F, header = T)
postSmallCA0401 = getPostSmall(id = "CA0401")
table(postSmallCA0401$anc)/nrow(postSmallCA0401)

# low A ancestry
postCA1410 = read.table(paste0(dir_post, "/", "CA1410", ".posterior"), stringsAsFactors = F, header = T)
postSmallCA1410 = getPostSmall(id = "CA1410")
table(postSmallCA1410$anc)/nrow(postSmallCA1410)
apply(mean_anc(postCA1410), 2, mean)

for (i in c("CA0401", "AR2711", "CA1410", "AR0302")){
  posti <- read.table(paste0(dir_post, "/", i, ".posterior"), 
                      stringsAsFactors = F, header = T)
  apply(mean_anc(posti), 2, mean)
}
bee_subset <- c("CA0401", "AR2711", "CA1410", "AR0302")
# for each bee, calculate mean ancestry
anc_all <- sapply(bees, function(id)
  apply(mean_anc( # get mean ancestry from posterior
    read.table(paste0(dir_post, "/", id, ".posterior"), 
                            stringsAsFactors = F, header = T)), 
        2, mean))
mean_anc_all <- as.data.frame(t(anc_all))
colnames(mean_anc_all) <- c("C", "M", "A")
mean_anc_all$Bee_ID <- rownames(mean_anc_all)

# plot ancestry one way
mean_anc_all %>%
  gather(., "anc", "mean", c("A", "C", "M")) %>%
  ggplot(., aes(x = Bee_ID, y = mean)) +
  geom_point(aes(color = anc))
# plot ancestry similar to NGSAdmix results to compare:
mean_anc_all %>%
  gather(., "anc", "mean", c("A", "C", "M")) %>%
  ggplot(., aes(fill=anc, y=mean, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + facet_wrap(~substr(Bee_ID, 1, 2))



# get all the posterior data
small_all = do.call(rbind,
                lapply(bees, function(i) 
  getPost(id = i)))

# plot all ind's ancestry (with uncertainty) over some region
small_all %>%
  filter(chrom == "Group1") %>%
  ggplot(aes(x=position, y = Bee_ID)) +
  geom_point(aes(color = anc, alpha = p), size = .5)

# plot all ind's ancestry for the whole genome, chromosome by chromosome
for (chrom in paste0("Group", 1:16)){
  small_all %>%
    filter(chrom == chrom) %>%
    ggplot(aes(x=position, y = Bee_ID)) +
    geom_point(aes(color = anc, alpha = p), size = .5) + 
    ggtitle(paste0("Ancestry 2018 bees -- ", chrom))
  ggsave(paste0("plots/local_ancestry_tracts_pass1_2018_", chrom, ".png"), 
         device = "png", 
         width = 20, height = 16, units = "in",
         dpi = 200)
}

# what confidence does the HMM have in the calls it makes?
hist(small_all$p)
summary(small_all$p)
# does the confidence vary by ancestry?
small_all %>%
  group_by(., anc) %>%
  summarize(mean_p = mean(p))
# very confident in all calls, but especially homozygous for ancestry calls.

