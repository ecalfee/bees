# plot global ancestry NGSadmix results
library(scales)
library(gtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggmap)
library(scatterpie)
library(OpenStreetMap)
library(maps)
library(mapproj)
library(geosphere)
library(mapdata)
library(viridis)
library(viridisLite)
library(RColorBrewer)
source("../colors.R") # for color palette
#Rio Claro, Sao Paulo Brazil: 22.4149° S, 47.5651° W (Google maps)
sao_paulo <- data.frame(long = -47.5651, lat = -22.4149)
# for maps:
my_api <- read.table("../maps/google_maps_api_EC2018.txt",
                     header = F, stringsAsFactors = F)$V1
ggmap::register_google(key = my_api)

# get metadata for all individuals included in NGSadmix analysis (plus extras)
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  mutate(continent = ifelse(group == "AR_2018", "S. America",
                            ifelse( group %in% c("CA_2018", "MX_2018",
                                                 "N_CA", "S_CA"),
                                    "N. America",
                            NA)))
meta$label <- meta$Bee_ID
pop2014_incl <- c("Stebbins_2014", "Stanislaus_2014", "Avalon_2014",
                  "Placerita_2014", "Riverside_2014", "Davis_2014")
for (i in pop2014_incl){ # make short labels for included bees from these ca_bee populations
  meta$label[meta$population == i & !is.na(meta$population)] <- 
    paste0(meta$geographic_location_short[meta$population == i & !is.na(meta$population)],
           "_",
           1:sum(meta$population == i, na.rm = T))
}
pop2014_2018_inc <- c(pop2014_incl, meta$population[meta$group %in% c("CA_2018", "AR_2018")])

#IDs <- read.table("../bee_samples_listed/with_duplicated_ap50/pass1.list", stringsAsFactors = F,
#                  header = F) # note ap50 is duplicated in the GL file output of ANGSD and NGSadmix results
#IDs <- read.table("../bee_samples_listed/pass1_plus_kohn.list", stringsAsFactors = F,
#                  header= F)
#IDs <- read.table("../bee_samples_listed/pass1_plus_kohn_and_wallberg.list", stringsAsFactors = F)
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

# which bees overlap with Dani's analysis?
bees_overlap_Dani <- c("CA0108", "CA0303", "AR1410", "SRCD9A", "SanDiego001", "SanDiego002", "Mexico001",
                       "SRR957075", "SRR957080", "SRR957061", "SRS549709") # plus one of each reference bee just for checking link to ancestry

# which bees are O from Jordan (admixed)?
bees_Jordan_O <- which(bees$geographic_location=="Jordan")
# I will cut out these rows from the GL file to re-run NGSAdmix.
rows_to_exclude_Jordan_O = sapply((which(bees$geographic_location=="Jordan"))*3, function(x) x + 1:3)
bees_noJordanO = bees[-bees_Jordan_O,]
#bees <- bees_noJordanO

# which A bees are not Wallberg A bees?
bees_A_not_wallberg <- which(bees$group=="A" & bees$source != "Wallberg")
# I will cut out these rows from the GL file to re-run NGSAdmix.
rows_to_exclude_not_Wallberg_A = sapply((bees_A_not_wallberg)*3, function(x) x + 1:3)
#bees_only_wallberg_A = bees[-bees_A_not_wallberg, ]
#bees <- bees_only_wallberg_A

# what if I only use the wallberg reference bees?
bees_ref_not_wallberg <- which(bees$group %in% c("A", "C", "M", "O") & bees$source != "Wallberg")
rows_to_exclude_not_Wallberg_ref = sapply((bees_ref_not_wallberg)*3, function(x) x + 1:3)
#bees_only_wallberg_ref = bees[-bees_ref_not_wallberg, ]
#bees <- bees_only_wallberg_ref


K = 3 # 3 admixing populations
#K = 4
#K = 5
#K = 6
colorsK=cbPalette[K]

# starting with pass1 analysis from 1st round of sequencing
#prefix1 = "ordered_scaffolds_prunedBy250"
#prefix1 = "ordered_scaffolds_prunedBy1000" # minor differences
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_prunedBy251"
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_prunedBy1000"
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_prunedBy251"
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_noJordanO_prunedBy1000"
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_noJordanO_prunedBy251"
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_onlyWallbergA_prunedBy251"
#prefix1 = "ordered_scaffolds_pass1_plus_kohn_and_wallberg_onlyWallbergREF_prunedBy251"
n = 250 # snps thinned to 1 every nth
prefix1 = paste0("ordered_scaffolds_", prefix, "_prunedBy", n)

name = paste0("K", K, "_", prefix1)
file = paste0("results/NGSAdmix/", name, ".qopt")
#admix <- read.table(file)[,3:1] # switched arbitrary order of ancestries to make visual comparison with pruneby250 easy
admix <- read.table(file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3)
# get allele freq. estimates for each ancestry (another output of NGSAdmix)
allele_freq_est <- read.table(paste0("results/NGSAdmix/", name, ".fopt.gz")) 
colnames(allele_freq_est) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3)

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d_admix <- bind_cols(bees, admix)  %>%
  arrange(., lat) %>%
  arrange(., source) %>%
  arrange(., group) %>%
  filter(., est_coverage > .05) # filters out one bee, AR1212, which had extremely low coverage -- I think it wasn't actually added to library pool
  
# plot 'STRUCTURE-like' ancestry plots
png(paste0("plots/", name, ".png"), # saves plot as ping in ../plots/
      height = 5, width = 8, units = "in", res = 150)
d_admix %>%
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

# what should the different ancestries be called? use the group with the highest frequency
anc <- data.frame(ancestry = colnames(admix),
                       ancestry_label = sapply(colnames(admix), function(x) names(which.max(tapply(d_admix[ , x], d_admix$population, sum)))),
                       stringsAsFactors = F) %>%
  mutate(ancestry_label = factor(ancestry_label, levels = c("C", "M", "A")))
# NOTE: now above code is specific to 3 ancestries because I wanted to order them A/M/C
# for plotting


# ggplot2 need 'tidy' formatted data
p1 <- d_admix %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~group) +
  scale_fill_manual(values = col_ACM, name = "Ancestry")
plot(p1)
ggsave(paste0("plots/NGS_admix_all_", name, ".png"), 
       plot = p1, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p2 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "CA_2018") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, y=p, x=reorder(Bee_ID, lat))) +
  geom_bar(stat = "identity", position = "fill") + 
  ggtitle("California 2018 bee samples") +
  xlab("Individual bee samples (ordered by latitude)") +
  ylab("Ancestry fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_ACM, name = "Ancestry")
plot(p2)
ggsave(paste0("plots/NGS_admix_CA_2018_", name, ".png"), 
       plot = p2, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)
# mexico
p2m <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "MX_2019") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, y=p, x=reorder(Bee_ID, lat))) +
  geom_bar(stat = "identity", position = "fill") + 
  ggtitle("Mexico 2019 bee samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_ACM, name = "Ancestry")
plot(p2m)
ggsave(paste0("plots/NGS_admix_MX_2019_", name, ".png"), 
       plot = p2m, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

#--------------------------------------------------------------------------------
# plots for manuscript:

# CA 2014 and 2018
p2_2014 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "CA_2018" | population %in% pop2014_incl) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, alpha = as.factor(year), y=p, x=reorder(label, lat))) +
  geom_bar(stat = "identity", position = "fill", width = 0.99) + # width=1 gets rid of lines, but some are still visible in pdf. better to make them uniform
  ggtitle("California 2014 & 2018 bee samples") +
  xlab("Individual bee samples (ordered by latitude)") +
  ylab("Ancestry fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_discrete(range = c(0.6, 1)) +
  labs(alpha = "Collection") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") + 
  theme_classic()
plot(p2_2014 )
ggsave(paste0("plots/NGS_admix_CA_2014_and_2018_", name, ".png"), 
       plot = p2_2014, 
       device = "png", 
       width = 17, height = 8, units = "in",
       dpi = 200)
ggsave(paste0("../../bee_manuscript/figures/NGS_admix_CA_2014_and_2018_", name, ".pdf"), 
       plot = p2_2014 +
         coord_flip() +
         theme(axis.ticks.y = element_blank(), # I can figure out how to add ticks at pop divisions later
               axis.text.y = element_blank()), 
       device = "pdf", 
       width = 4, height = 8, units = "in",
       dpi = 200)



p3 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "AR_2018") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill = ancestry_label, alpha = as.factor(year), y = p, x = reorder(Bee_ID, lat))) +
  geom_bar(stat = "identity", position = "fill", width = 0.99) + 
  ggtitle("Argentina 2018 bee samples") + 
  xlab("Individual bee samples (ordered by latitude)") +
  ylab("Ancestry fraction") +
  scale_alpha_discrete(range = c(1, 1)) +
  labs(alpha = "Collection") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  theme_classic()
plot(p3)
ggsave(paste0("plots/NGS_admix_AR_2018_", name, ".png"), 
       plot = p3, 
       device = "png", 
       width = 17, height = 8, units = "in",
       dpi = 200)
ggsave(paste0("../../bee_manuscript/figures/NGS_admix_AR_2018_", name, ".pdf"), 
       plot = p3 +
         coord_flip() +
         theme(axis.ticks.y = element_blank(), # I can figure out how to add ticks at pop divisions later
               axis.text.y = element_blank()), 
       device = "pdf", 
       width = 4, height = 8, units = "in",
       dpi = 200)

# make map with pie charts
admix.pops <- d_admix %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(population %in% pop2014_2018_inc) %>%
  left_join(., anc, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  tidyr::spread(., ancestry_label, p) %>%
  group_by(population) %>%
  summarise(A = mean(A),
            M = mean(M),
            C = mean(C),
            long = mean(long),
            lat = mean(lat),
            year = mean(year),
            n = n()) %>%
  mutate(continent = ifelse(lat < 0, "S. America", "N. America")) %>%
  arrange(lat) %>%
  mutate(population_factor = factor(population, ordered = T, levels = .$population)) %>% # get ordered by lat
  mutate(shape = c(rep(1:6, 6), 1:3))
  

#autoplot.osmtime()
tile_america <- openmap(upperLeft = c(47.52, -126.75),
                        lowerRight = c(-44.47, -33.05),
                        type = "osm",
                        mergeTiles = T)
tile_america2 <- openmap(upperLeft = c(45.0, -127.0),
                        lowerRight = c(-45.0, -33.0),
                        type = "bing",
                        mergeTiles = T)
tile_america3 <- openmap(upperLeft = c(47, -126.75),
                         lowerRight = c(-46, -33.05),
                         type = "esri",
                         mergeTiles = T)
tile_america4 <- openmap(upperLeft = c(48, -145),
                         lowerRight = c(-45, -32),
                         type = "esri",
                         mergeTiles = T)
autoplot(tile_america3) +
  geom_point(aes(x = long*10^5,
                 y = lat*10^5), #,
             #col = popN), 
             data = filter(meta, population %in% pop2014_2018_inc),
             cex = .25,
             col = "black",
             alpha = .5)
p_world <- autoplot(openproj(tile_america3)) + # squishes map horizontally
  geom_point(aes(x = long,
                 y = lat), #,
             #col = popN), 
             data = filter(meta, population %in% pop2014_2018_inc),
             cex = .25,
             col = "black",
             alpha = .5)
#world2 <- map_data("world2Hires")
world2 <- map_data("world2")
world <- map_data("world")
states <- map_data("state")
dim(states)
ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, fill = region, group = group), color = "white") + 
  coord_fixed(1.3) +
  guides(fill=FALSE)

p_world <- world %>%
  filter(., ! region %in% c("Greenland", "Antarctica", "Canada")) %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               alpha = .3,
               #alpha = .7,
               fill = dark2[8]
               ) +
  xlim(c(-140, -20)) +
  coord_fixed(1.1, 
              xlim = c(-126, -23), 
              ylim = c(-39, 39.5)) +
  theme_classic() +
  geom_point(data = sao_paulo, 
             aes(x = long, y = lat),
             color = col_ACM["A"],
             shape = 8,
             size = 2) +
  geom_point(aes(x = long,
                 y = lat,
                  col = continent), #,
             #col = popN), 
             data = admix.pops,
             cex = .25,
             alpha = .5) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme(legend.position = "None")
plot(p_world)
ggsave("../../bee_manuscript/figures/world_map_samples.pdf", # ugh. blurry (!)
       plot = p_world, 
       device = "pdf", 
       width = 6, height = 6, units = "in")
# just SA samples:
SA_pie <- world %>%
  filter(., ! region %in% c("Greenland", "Antarctica", "Canada")) %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               #alpha = .3,
               alpha = .7,
               fill = dark2[8]
  ) +
  xlim(c(-140, -20)) +
  coord_fixed(1.1, 
              xlim = c(-62, -57), 
              ylim = c(-36.5, -28.5)) +
  theme_classic() +
  geom_point(data = sao_paulo, 
             aes(x = long, y = lat),
             color = col_ACM["A"],
             shape = 8,
             size = 2) +

  geom_scatterpie(data = admix.pops, 
                  aes(long, lat, r = .15),
                  cols = c("A", "C", "M"), 
                  alpha = 1,
                  lwd = 0) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_manual(values = col_ACM, name = NULL) +
  theme(legend.position = "None") +
  ggtitle("Argentina")
plot(SA_pie)
ggsave("../../bee_manuscript/figures/SA_pie_map_ancestry.pdf",
       plot = SA_pie, 
       device = "pdf", 
       width = 3, height = 6, units = "in")
# N. America only
NA_pie <- world %>%
  filter(., ! region %in% c("Greenland", "Antarctica", "Canada")) %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               #alpha = .3,
               alpha = .7,
               fill = dark2[8]
  ) +
  xlim(c(-140, -20)) +
  coord_fixed(1.3, 
              xlim = c(-122, -116), 
              ylim = c(32, 38.5)) +
  theme_classic() +
  geom_scatterpie(data = admix.pops, 
                  aes(long, lat, r = .15),
                  cols = c("A", "C", "M"), 
                  alpha = 1,
                  lwd = .1,
                  color = "white") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_manual(values = col_ACM, name = NULL) +
  theme(legend.position = "None") +
  ggtitle("California")
plot(NA_pie)
ggsave("../../bee_manuscript/figures/NA_pie_map_ancestry.pdf",
       plot = NA_pie, 
       device = "pdf", 
       width = 3, height = 6, units = "in")

# SA points, not pie chart
SA_points <- world %>%
  filter(., ! region %in% c("Greenland", "Antarctica", "Canada")) %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               #alpha = .3,
               #alpha = .7,
               fill = dark2[8]
  ) +
  xlim(c(-140, -20)) +
  theme_classic() +
  geom_point(aes(x = long,
                 y = lat,
                 color = population,
                 shape = factor(shape)), 
             data = filter(admix.pops, continent == "S. America"),
             cex = 2) +
             #alpha = .5) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis_d(option = "magma", direction = -1) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  theme(legend.position = "None") +
  ggtitle("Argentina")
plot(SA_points +
       coord_fixed(1.1, 
                   xlim = c(-62, -57),
                   ylim = c(-36.5, -28.5)))
plot(SA_points +
       coord_fixed(1.1, 
                   xlim = c(-62, -55),
                   ylim = c(-36.5, -27)))
ggsave("../../bee_manuscript/figures/SA_point_map_samples_zoom_in.pdf",
       plot = SA_points + 
         coord_fixed(1.1, 
                     xlim = c(-62, -57),
                     ylim = c(-36.5, -28.5)), 
       device = "pdf", 
       width = 3, height = 6, units = "in")
ggsave("../../bee_manuscript/figures/SA_point_map_samples_zoom_out.pdf",
       plot = SA_points +
         coord_fixed(1.1, 
                     xlim = c(-62, -55),
                     ylim = c(-36.5, -27)), 
       device = "pdf", 
       width = 3, height = 6, units = "in")

# NA just colored points, not pie charts:
NA_points <- world %>%
  filter(., ! region %in% c("Greenland", "Antarctica", "Canada")) %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               #alpha = .3,
               #alpha = .7,
               fill = dark2[8]
  ) +
  xlim(c(-140, -20)) +
  theme_classic() +
  geom_point(aes(x = long,
                 y = lat,
                 color = population_factor,
                 shape = factor(shape)), 
             data = filter(admix.pops, continent == "N. America"),
             cex = 2) +
  #alpha = .5) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis_d(option = "viridis") +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  theme(legend.position = "None") +
  ggtitle("California")
plot(NA_points +
       coord_fixed(1.3, 
                   xlim = c(-122, -116), 
                   ylim = c(32, 38.5)))
plot(NA_points +
       coord_fixed(1.3, 
                   xlim = c(-124, -114.5), 
                   ylim = c(32, 39.5)))
ggsave("../../bee_manuscript/figures/NA_point_map_samples_zoom_in.pdf",
       plot = NA_points +
                     coord_fixed(1.3, 
                                 xlim = c(-122, -116), 
                                 ylim = c(32, 38.5)),
       device = "pdf", 
       width = 3, height = 6, units = "in")
ggsave("../../bee_manuscript/figures/NA_point_map_samples_zoom_out.pdf",
       NA_points +
         coord_fixed(1.3, 
                     xlim = c(-124, -114.5), 
                     ylim = c(32, 39.5)), 
       device = "pdf", 
       width = 3, height = 6, units = "in")
AR_end1 <- filter(admix.pops,
                  continent == "S. America") %>%
  filter(., lat == min(lat))
AR_end2 <- filter(admix.pops,
                  continent == "S. America") %>%
  filter(., lat == max(lat))
distm(AR_end1[ , c("long", "lat")], AR_end2[ , c("long", "lat")], fun = distGeo)

# --------------------------------------------------------------------------

p4 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "S_CA" | group == "Kohn") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~population) +
  ggtitle("Southern CA (and Mexico) 1994-2015 samples") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  scale_fill_manual(values = col_ACM, name = "Ancestry")
plot(p4)
ggsave(paste0("plots/NGS_admix_S_CA_Mex_1994-2015_", name, ".png"), 
       plot = p4, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p5 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "N_CA") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~population) +
  ggtitle("Northern CA 1994-2015 samples") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  scale_fill_manual(values = col_ACM, name = "Ancestry")
plot(p5)
ggsave(paste0("plots/NGS_admix_N_CA_1994-2015_", name, ".png"), 
       plot = p5, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

p6 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group %in% c("A", "C", "M", "O")) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, y=p, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_ACM, name = "Ancestry")
plot(p6 + facet_wrap(~group))
ggsave(paste0("plots/NGS_admix_Ref_", name, ".png"), 
       plot = p6 + facet_wrap(~group), 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

plot(p6 + facet_wrap(~source))
ggsave(paste0("plots/NGS_admix_Ref_bySource", name, ".png"), 
       plot = p6 + facet_wrap(~source), 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)

# print file with population frequencies of each admixed population 
# to use as priors in local ancestry inference
d_ACM <- d_admix %>%
  filter(source %in% c("Calfee", "Ramirez")) %>%
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  left_join(., anc, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  spread(ancestry_label, p)
d_pop_ACM <- d_ACM %>%
  group_by(population) %>%
  summarise(A = mean(A),
            C = mean(C),
            M = mean(M),
            n = n())
d_pop_ACM %>%
  filter(A > .001) %>% # do not run ancestry inference on populations with fewer than 4 individuals or extremely low A ancestry
  write.table(., paste0("results/NGSAdmix/", name, ".pop.anc"),
              quote = F, col.names = T, row.names = F, sep = "\t")
# write out individual ancestry too
d_ACM %>%
  filter(., population %in% unique(d_pop_ACM$population[d_pop_ACM$A > .001])) %>% # do not run ancestry inference on populations with fewer than 4 individuals or extremely low A ancestry
  dplyr::select(Bee_ID, population, A, C, M) %>%
  write.table(., paste0("results/NGSAdmix/", name, ".ind.anc"),
              quote = F, col.names = T, row.names = F, sep = "\t")

# What is estimated Fst between these ancestries according to estimated ancestry allele frequencies?
allele_freq_est1 <- allele_freq_est %>%
  mutate(anc1_anc2 = (anc1+anc2)/2) %>%
  mutate(anc1_anc3 = (anc1+anc3)/2) %>%
  mutate(anc2_anc3 = (anc2+anc3)/2) %>%
  mutate(total = (anc1+anc2+anc3)/3)
allele_het_est1 <- 2*allele_freq_est1*(1-allele_freq_est1)
allele_het_est <- apply(allele_het_est, 2, mean)
1-mean(allele_het_est[c("anc1", "anc2")])/allele_het_est["anc1_anc2"]
1-mean(allele_het_est[c("anc1", "anc3")])/allele_het_est["anc1_anc3"]
1-mean(allele_het_est[c("anc2", "anc3")])/allele_het_est["anc2_anc3"]
# fairly high Fst between groups (as expected)
1-allele_het_est/allele_het_est["total"]

# only makes sense for K4 with O group -- how much are O and C ancestry correlated? 
# possibly (?) positive correlation O-C for low C and neg correlation for high C .. but maybe should be neg by design for high C b/c its a percent
ggplot(d, aes(anc3, anc2, color = group)) + geom_point() + facet_wrap(~source)

# look at ancestry estimates for a few bees
# to send to Dani for comparison with her results
# note: ancestry translations only work for 251, k=4, w/ Kohn & Wallberg:
d_admix %>%
  mutate(A = round(anc1, 2), C = round(anc2, 2),
         O = round(anc3, 2), M = round(anc4, 2)) %>%
  filter(Bee_ID %in% bees_overlap_Dani) %>%
  arrange(source) %>%
  select(c("Bee_ID", "geographic_location", "year", "group", "source", "A", "C", "M", "O")) %>%
  write.table(., paste0("plots/NGS_admix_results_subset_", name, ".txt"),
              quote = F, row.names = F, sep = "\t")
# note: ancestry translations only work for 251, k=3, w/ Kohn & Wallberg:
d_admix %>%
  mutate(A = round(anc3, 2), C = round(anc2, 2),
         M = round(anc1, 2)) %>%
  filter(Bee_ID %in% bees_overlap_Dani) %>%
  arrange(source) %>%
  select(c("Bee_ID", "geographic_location", "year", "group", "source", "A", "C", "M")) %>%
  write.table(., paste0("plots/NGS_admix_results_subset_", name, ".txt"),
              quote = F, row.names = F, sep = "\t")
# what does the distribution of SNPs look like?
sites <- read.table(paste0("results/input/", prefix, ".var.sites")) %>%
  tidyr::separate(., V1, sep = "_", into = c("scaffold", "pos")) %>%
  dplyr::mutate(pos = as.numeric(pos))
sites$diff <- diff(c(0,sites$pos))
summary(sites$diff[sites$diff >= 0])
summary(sites$diff[sites$diff >= 250])
hist(sites$diff[sites$diff >= 0 & sites$diff < 251])


# Ancestry translation true if using 
# "K3_ordered_scaffolds_pass1_plus_kohn_prunedBy251" data
filter(d, source == "Calfee") %>%
  mutate(M_ancestry = anc1) %>%
  mutate(C_ancestry = anc2) %>%
  mutate(A_ancestry = anc3) %>%
  select(Bee_ID, geographic_location, population, year, lat, long, 
         popN, indN, M_ancestry, C_ancestry, A_ancestry) %>%
  write.table(., paste0("results/A_ancestry_estimates_for_Jodie_", name, ".txt"),
              quote = F, col.names = T, row.names = F, sep = "\t")
lm.a.lat <- with(filter(d, source == "Calfee"), lm(anc3 ~ abs(lat)*geographic_location))
lm.logita.lat <- with(filter(d, source == "Calfee"), lm(gtools::logit(anc3) ~ abs(lat)*geographic_location))
summary(lm.a.lat)
summary(lm.logita.lat)

# quickly plot latitude trends
filter(d, source == "Calfee") %>%
  ggplot(., aes(x = abs(lat), y = anc3, color = geographic_location)) +
  geom_point() +
  xlab("Latitude (abs value)") +
  ylab("Percent African ancestry") +
  facet_wrap(~geographic_location) +
  ggtitle("A ancestry ~ latitude across 2 Africanized honey bee clines")
ggsave(paste0("plots/A_ancestry_by_latitude_", name, ".png"),
       device = "png",
       width = 10, height = 8, units = "in",
       dpi = 200)


# testing effect of recomb. bin on global ancestry estimate
#(1) get ancestry estimates for A ancestry across recomb bins
byR_list <- vector("list", 5)  # empty list w/ 5 bins
for (i in 1:5){
  p5 <- paste0("recomb_10kb_5bins_", i, "_ordered_scaffolds_CA_AR_MX_harpur_sheppard_kohn_wallberg_NoGroupUn_prunedBy10")
  name5 <- paste0("K", K, "_", p5)
  file5 <- paste0("results/NGSAdmix/", name5, ".qopt")
  admix5 <- read.table(file5)
  colnames(admix5) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3)
  # join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
  d5 <- bind_cols(bees, admix5)  %>%
    arrange(., lat) %>%
    arrange(., source) %>%
    arrange(., group) %>%
    filter(., est_coverage > .05) # filters out one bee, AR1212, which had extremely low coverage -- I think it wasn't actually added to library pool
  # what should the different ancestries be called? use the group with the highest frequency
  anc5 <- data.frame(ancestry = colnames(admix5),
                    ancestry_label = sapply(colnames(admix5), function(x) names(which.max(tapply(d5[ , x], d5$population, mean)))),
                    stringsAsFactors = F)
  byR_list[[i]] <- d5 %>%
    mutate(., rbin = i) %>%
    tidyr::gather(., "ancestry", "p", colnames(admix5)) %>%
    left_join(., anc5, by = "ancestry") # save in my list
}
byR <- do.call(rbind, byR_list)
# plot mean ancestry by r bin across CA/AR clines
byR %>%
  filter(source == "Calfee" & geographic_location %in% c("Argentina", "California")) %>% 
  #filter(geographic_location == "Argentina") %>%
  filter(ancestry_label == "A") %>%
  ggplot(., aes(x = abs(lat), y = p, color = rbin, group = rbin)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~geographic_location) +
  ggtitle("no clear differences in ancestry slope for high (5) vs low (1) recomb. bins")
ggsave("plots/slope_anc_lat_by_rbin5.png", device = "png",
       width = 8, height = 6, units = "in")

# plot mean ancestry for populations across r bins
byR %>%
  filter(source == "Calfee" & geographic_location %in% c("Argentina", "California")) %>% 
  filter(geographic_location == "Argentina") %>%
  filter(ancestry_label == "A") %>%
  mutate(rbin = as.factor(rbin)) %>%
  ggplot(., aes(x = rbin, y = p, fill = rbin)) +
  geom_boxplot() +
  facet_wrap(~population, scales = "free_y") +
  ggtitle("Argentina: no clear patterns in A ancestry for high (5) vs low (1) recomb. bins")
ggsave("plots/boxplot_A_anc_Argentina_pop_by_rbin5.png", device = "png",
       width = 8, height = 6, units = "in")

byR %>%
  filter((source == "Calfee" & geographic_location == "California") |
  source == "Ramirez") %>% 
  filter(ancestry_label == "A") %>%
  mutate(rbin = as.factor(rbin)) %>%
  ggplot(., aes(x = rbin, y = p, fill = rbin)) +
  geom_boxplot() +
  facet_wrap(~population, scales = "free_y") +
  ggtitle("California: no clear patterns in A ancestry for high (5) vs low (1) recomb. bins")
ggsave("plots/boxplot_A_anc_California_pop_by_rbin5.png", device = "png",
       width = 8, height = 6, units = "in")


# plot reference bees ancestry across recomb. bins
byR %>%
  filter(source == "Harpur" | source == "Sheppard") %>% 
  mutate(rbin = as.factor(rbin)) %>%
  ggplot(., aes(x = rbin, y = p, fill = ancestry_label)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free_y") +
  ggtitle("reference A/C/M ancestries clearly distinguished across high (5) vs low (1) recomb. bins")
ggsave("plots/ref_bees_ACM_well_identifiable_across_rbin5.png", device = "png",
       width = 8, height = 6, units = "in")
