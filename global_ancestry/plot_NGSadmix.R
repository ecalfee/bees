# plot global ancestry NGSadmix results
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette

# get metadata for individuals included in NGSadmix analysis
meta1 <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")
# NOTE: quick fix to view kohn bees (need meta):
kohn_meta <- data.frame(Bee_ID = c("SanDiego001", "SanDiego002", "Mexico001"),
                        source = "Kohn",
                        strain = "unknown",
                        year = 2015,
                        group = "Kohn",
                        geographic_location = c("San Diego", "San Diego", "Mexico"),
                        population = c("SanDiego_2015", "SanDiego_2015", "Mexico_2015"),
                        stringsAsFactors = F)
# get wallberg meta data
wallberg_ACMO = data.frame(geographic_location = c("Italy", "Austria", 
                             "Jordan", "Turkey",
                             "South Africa", "Nigeria",
                             "Norway", "Sweden", "Spain"),
                           group = c("C", "C",
                                     "O", "O",
                                     "A", "A",
                                     "M", "M", "M"),
                           stringsAsFactors = F)
wallberg_meta <- read.table("../bee_samples_listed/Wallberg_2014_all.meta",
                            header = T, sep = "\t", stringsAsFactors = F) %>%
  unique(.) %>%
  dplyr::select(c("SRA_Sample", "geo_loc_name")) %>%
  rename(Bee_ID = SRA_Sample,
         geographic_location = geo_loc_name) %>%
  mutate(source = "Wallberg",
         year = 2014) %>%
  left_join(., wallberg_ACMO, by = "geographic_location")

# merge all meta data
meta <- bind_rows(meta1, kohn_meta, wallberg_meta)

# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
#IDs <- read.table("../bee_samples_listed/with_duplicated_ap50/pass1.list", stringsAsFactors = F,
#                  header = F) # note ap50 is duplicated in the GL file output of ANGSD and NGSadmix results
#IDs <- read.table("../bee_samples_listed/pass1_plus_kohn.list", stringsAsFactors = F,
#                  header= F)
IDs <- read.table("../bee_samples_listed/pass1_plus_kohn_and_wallberg.list", stringsAsFactors = F)
colnames(IDs) <- c("Bee_ID")
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") 

# which bees are O from Jordan (admixed)?
bees_Jordan_O <- which(bees$geographic_location=="Jordan")
# I will cut out these rows from the GL file to re-run NGSAdmix.
rows_to_exclude_Jordan_O = sapply((which(bees$geographic_location=="Jordan"))*3, function(x) x + 1:3)
bees_noJordanO = bees[-bees_Jordan_O,]
bees <- bees_noJordanO

#K = 3 # 3 admixing populations
K = 4
colorsK=c("red", "cornflowerblue", "navy")

# starting with pass1 analysis from 1st round of sequencing
#prefix = "ordered_scaffolds_prunedBy250"
#prefix = "ordered_scaffolds_prunedBy1000" # minor differences
#prefix = "ordered_scaffolds_pass1_plus_kohn_prunedBy251"
#prefix = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_prunedBy1000"
#prefix = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_prunedBy251"
prefix = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_noJordanO_prunedBy251"

name = paste0("K", K, "_", prefix)
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
  dplyr::select(., colnames(admix)) %>%
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
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + facet_wrap(~group)
plot(p1)
ggsave(paste0("plots/NGS_admix_all_", name, ".png"), 
       plot = p1, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p2 <- d %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "CA_2018") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  ggtitle("California 2018 bee samples")
plot(p2)
ggsave(paste0("plots/NGS_admix_CA_2018_", name, ".png"), 
       plot = p2, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p3 <- d %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "AR_2018") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  ggtitle("Argentina 2018 bee samples")
plot(p3)
ggsave(paste0("plots/NGS_admix_AR_2018_", name, ".png"), 
       plot = p3, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p4 <- d %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "S_CA" | group == "Kohn") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~population) +
  ggtitle("Southern CA (and Mexico) 1994-2015 samples")
plot(p4)
ggsave(paste0("plots/NGS_admix_S_CA_Mex_1994-2015_", name, ".png"), 
       plot = p4, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p5 <- d %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "N_CA") %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~population) +
  ggtitle("Northern CA 1994-2015 samples")
plot(p5)
ggsave(paste0("plots/NGS_admix_N_CA_1994-2015_", name, ".png"), 
       plot = p5, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p6 <- d %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group %in% c("A", "C", "M")) %>%
  ggplot(., aes(fill=ancestry, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~group)
plot(p6)
ggsave(paste0("plots/NGS_admix_A_C_M_Ref_", name, ".png"), 
       plot = p6, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

