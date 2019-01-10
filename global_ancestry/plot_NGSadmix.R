# plot global ancestry NGSadmix results
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette

# get metadata for individuals included in NGSadmix analysis
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")

# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table("../bee_samples_listed/pass1.list", stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") 
# NOTE: bee ap50 is a duplicate entry I will fix later (!)

K=3 # 3 admixing populations
colorsK=c("red", "cornflowerblue", "navy")

# starting with pass1 analysis from 1st round of sequencing
name = paste0("K", K, "_ordered_scaffolds_prunedBy250")
#name = paste0("K", K, "_ordered_scaffolds_prunedBy1000") # minor differences
file = paste0("results/NGSAdmix/", name, ".qopt")
#admix <- read.table(file)[,3:1] # switched arbitrary order of ancestries to make visual comparison with pruneby250 easy
admix <- read.table(file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3)

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(bees, admix)  %>%
  arrange(., lat) %>%
  arrange(., source) %>%
  arrange(., group)
  
# plot 'STRUCTURE-like' ancestry plots
png(paste0("plots/", name, ".png"), # saves plot as ping in ../plots/
      height = 5, width = 8, units = "in", res = 150)
d %>%
  select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA,
          main = "all bees",
          ylab="admixture") 
title(xlab = " A  |  Argentina '18  |  C  |  California '18  |  M  |  N & S CA 1994-2015")
dev.off()

# ggplot2 need 'tidy' formatted data
p1 <- d %>% 
  tidyr::gather(., "ancestry", "p", c("anc1", "anc2", "anc3")) %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + facet_wrap(~group)
plot(p1)
ggsave("plots/NGS_admix_all.png", 
       plot = p1, 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

p2 <- d %>% tidyr::gather(., "ancestry", "p", c("anc1", "anc2", "anc3")) %>%
  filter(group == "CA_2018") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  ggtitle("California 2018 bee samples")
plot(p2)
ggsave("plots/NGS_admix_CA_2018.png", 
       plot = p2, 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

p3 <- d %>% tidyr::gather(., "ancestry", "p", c("anc1", "anc2", "anc3")) %>%
  filter(group == "AR_2018") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  ggtitle("Argentina 2018 bee samples")
plot(p3)
ggsave("plots/NGS_admix_AR_2018.png", 
       plot = p3, 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

p4 <- d %>% tidyr::gather(., "ancestry", "p", c("anc1", "anc2", "anc3")) %>%
  filter(group == "S_CA") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~population) +
  ggtitle("Southern CA 1994-2015 samples")
plot(p4)
ggsave("plots/NGS_admix_S_CA_1994-2015.png", 
       plot = p4, 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

p5 <- d %>% tidyr::gather(., "ancestry", "p", c("anc1", "anc2", "anc3")) %>%
  filter(group == "N_CA") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~population) +
  ggtitle("Northern CA 1994-2015 samples")
plot(p5)
ggsave("plots/NGS_admix_N_CA_1994-2015.png", 
       plot = p5, 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

p6 <- d %>% tidyr::gather(., "ancestry", "p", c("anc1", "anc2", "anc3")) %>%
  filter(group %in% c("A", "C", "M")) %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~group)
plot(p6)
ggsave("plots/NGS_admix_A_C_M_Ref.png", 
       plot = p6, 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

