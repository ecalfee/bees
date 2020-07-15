# plot maps and global ancestry NGSadmix results
library(dplyr)
library(tidyr)
library(ggplot2)
library(scatterpie)
library(geosphere)
library(viridisLite)
library(gridExtra)
library(cowplot)

source("../colors.R") # for color palette

#Rio Claro, Sao Paulo Brazil: 22.4149° S, 47.5651° W (Google maps)
sao_paulo <- data.frame(long = -47.5651, lat = -22.4149)

# get metadata for all individuals included in NGSadmix analysis (plus extras)
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  mutate(continent = ifelse(group == "AR_2018", "S. America",
                            ifelse( group %in% c("CA_2018", "N_CA", "S_CA"),
                                    "N. America",
                            NA))) %>%
  dplyr::select(-c(Date, Time))

meta$label <- meta$Bee_ID

load("../local_ancestry/results/meta.RData") # meta.ind and meta.pop

pop2014_incl <- c("Stebbins_2014", "Stanislaus_2014", "Avalon_2014",
                  "Placerita_2014", "Riverside_2014", "Davis_2014")
for (i in pop2014_incl){ # make short labels for included bees from these ca_bee populations
  meta$label[meta$population == i & !is.na(meta$population)] <- 
    paste0(meta$geographic_location_short[meta$population == i & !is.na(meta$population)],
           "_",
           1:sum(meta$population == i, na.rm = T))
}
pop2014_2018_inc <- c(pop2014_incl, unique(meta$population[meta$group %in% c("CA_2018", "AR_2018") & meta$toSequence]))
write.table(pop2014_2018_inc, "../bee_samples_listed/byPop/combined_sept19_pops.list",
            col.names = F, row.names = F, quote = F, sep = "\t")

prefix <- "combined_sept19"
# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table(paste0("../bee_samples_listed/", prefix, ".list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")

# get coverage estimates
coverage <- read.table(paste0("../geno_lik_and_SNPs/results/", prefix, "/coverage/mean_ind_coverage.chr.random_pos_1000.txt"),
                       header = T, stringsAsFactors = F, sep = "\t")
# mean coverage all samples:
left_join(IDs, meta, by = "Bee_ID") %>%
  left_join(., coverage, by = "Bee_ID") %>% group_by(source) %>% summarise(mean = mean(est_coverage))


# join all data together
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") %>%
  dplyr::left_join(., coverage, by = "Bee_ID")

K = 3 # 3 admixing populations

# starting with pass1 analysis from 1st round of sequencing
n = 250 # snps thinned to 1 every nth
prefix1 = paste0(prefix, "_chr_prunedBy", n)

name = paste0("K", K, "_", prefix1)
file = paste0("results/NGSAdmix/", name, ".qopt")
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
  

# what should the different ancestries be called? use the group with the highest frequency
anc <- data.frame(ancestry = colnames(admix),
                       ancestry_label = sapply(colnames(admix), function(x) names(which.max(tapply(d_admix[ , x], d_admix$population, sum)))),
                       stringsAsFactors = F) %>%
  mutate(ancestry_label = factor(ancestry_label, levels = c("C", "M", "A")))
# NOTE: now above code is specific to 3 ancestries because I wanted to order them A/M/C
# for plotting
d_admix_ACM_labelled <- d_admix %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  left_join(., anc, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  tidyr::spread(., ancestry_label, p)
save(d_admix_ACM_labelled, file = paste0("results/NGSAdmix/ACM_", name, ".rData"))

# ------------------------------------------
# plot ngsadmix results and maps for sampling and historical invasion

# NGS Admix results for reference pops
d_ref_places <- d_admix %>% 
  filter(population %in% c("C", "M", "A")) %>%
  mutate(geographic_location_short = ifelse(geographic_location_short == "Pag Island, Croatia", "Croatia", geographic_location_short)) %>%
  mutate(place = paste(population, geographic_location_short)) %>%
  arrange(place) %>%
  mutate(order = 1:nrow(.) + as.integer(factor(place))) 
p_admix_ref_places <- d_ref_places %>%
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(y = p, x = order, labels = place)) +
  geom_bar(aes(fill = ancestry_label), stat = "identity", position = "fill", width = 0.95) + # width=1 gets rid of lines, but some are still visible in pdf. better to make them uniform
  ylab("Ancestry fraction") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") + 
  theme_classic() +
  scale_x_continuous(breaks = (group_by(d_ref_places, place) %>% summarise(pos = mean(order)))$pos, 
                     labels = unique(d_ref_places$place), name = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/NGS_admix_refACM.png", 
       plot = p_admix_ref_places,
       device = "png", 
       width = 5.2, height = 4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/NGS_admix_refACM.png", 
       plot = p_admix_ref_places,
       device = "png", 
       width = 5.2, height = 4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/NGS_admix_refACM.tiff", 
       plot = p_admix_ref_places,
       device = "tiff", 
       width = 5.2, height = 4, units = "in", dpi = 600)

# NGSAdmix for samples from hybrid zones
# CA 2014 and 2018
p2_2014 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "CA_2018" | population %in% pop2014_incl) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, alpha = as.factor(year), y=p, x=reorder(label, lat))) +
  geom_bar(stat = "identity", position = "fill", width = 0.99) + # width=1 gets rid of lines, but some are still visible in pdf. better to make them uniform
  
  ylab("Ancestry fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_discrete(range = c(0.6, 1)) +
  labs(alpha = "Collection") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") + 
  xlab("Individual bee samples (ordered by latitude)") +
  theme_classic() 
plot(p2_2014 )
ggsave(paste0("plots/NGS_admix_CA_2014_and_2018_", name, ".png"), 
       plot = p2_2014 +
         ggtitle("California 2014 & 2018 bee samples"), 
       device = "png", 
       width = 17, height = 8, units = "in",
       dpi = 200)
ggsave(paste0("../../bee_manuscript/figures/NGS_admix_CA_2014_and_2018_", name, ".pdf"), 
       plot = p2_2014 +
         coord_flip() +
         theme(axis.title.y = element_blank(),
               axis.ticks.y = element_blank(), # I can figure out how to add ticks at pop divisions later
               axis.text.y = element_blank()), 
       device = "pdf", 
       width = 4, height = 8, units = "in",
       dpi = 200)

p3 <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "AR_2018") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill = ancestry_label, alpha = as.factor(year), y = p, x = reorder(Bee_ID, lat))) +
  geom_bar(stat = "identity", position = "fill", width = 0.99) + 
  xlab("Individual bee samples (ordered by latitude)") +
  ylab("Ancestry fraction") +
  scale_alpha_discrete(range = c(1, 1)) +
  labs(alpha = "Collection") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  theme_classic()
plot(p3)
ggsave(paste0("plots/NGS_admix_AR_2018_", name, ".png"), 
       plot = p3 + ggtitle("Argentina 2018 bee samples"), 
       device = "png", 
       width = 17, height = 8, units = "in",
       dpi = 200)
ggsave(paste0("../../bee_manuscript/figures/NGS_admix_AR_2018_", name, ".pdf"), 
       plot = p3 +
         coord_flip() +
         theme(axis.title.y = element_blank(),
               axis.ticks.y = element_blank(), # I can figure out how to add ticks at pop divisions later
               axis.text.y = element_blank()), 
       device = "pdf", 
       width = 4, height = 8, units = "in",
       dpi = 200)

# summarise admixture by pop (e.g. pie chart)
admix.ind <- d_admix %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(population %in% pop2014_2018_inc) %>%
  left_join(., anc, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  tidyr::spread(., ancestry_label, p) 
admix.pops <- admix.ind %>%
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
  
# get maps
world <- map_data("world") # world data, country outlines
states <- map_data("state") # united states data

p_world <- world %>%
  filter(., ! region %in% c("Greenland", "Antarctica", "Canada")) %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               alpha = .3,
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
                  col = continent),
             data = admix.pops,
             cex = .25,
             alpha = 1) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme(legend.position = "None")
buffer = 1
arrow_NA <- data.frame(x1 = sao_paulo$long - buffer, 
                       y1 = sao_paulo$lat + buffer*1.25, 
                       x2 = -114.5 + buffer*2, 
                       y2 = 32)
arrow_SA <- data.frame(x1 = sao_paulo$long - buffer*1.5, 
                       y1 = sao_paulo$lat - buffer*.5, 
                       x2 = mean(c(5-55, -62)), 
                       y2 = -27 + buffer)
p_world_arrows <- p_world + 
       # draw SA rectangle
       geom_polygon(data = data.frame(X = c(-62, -55, -55, -62), 
                                      Y = c(-37, -37, -27, -27)), 
                    aes(x = X, y = Y), fill = NA, lwd = 0.75,
                    color = col_NA_SA_both["S. America"]) +
       # draw NA rectangle
       geom_polygon(data = data.frame(X = c(-124, -114.5, -114.5, -124), 
                                      Y = c(32, 32, 39.5, 39.5)),
                    aes(x = X, y = Y), fill = NA, lwd = 0.75,
                    color = col_NA_SA_both["N. America"]) +
       # draw curved arrows for routes of expansion
       geom_curve(aes(x = x1, y = y1, 
                      xend = x2, yend = y2),
                      color = col_ACM["A"],
                      lwd = 0.75,
                      data = arrow_NA, 
                      curvature = 0.1,
                      arrow = arrow(length = unit(0.1, "inches"))) +
       geom_curve(aes(x = x1, y = y1, 
                 xend = x2, yend = y2),
                 color = col_ACM["A"],
                 lwd = 0.75,
                 data = arrow_SA, 
                 curvature = 0.05,
                 arrow = arrow(length = unit(0.1, "inches")))
plot(p_world_arrows)
ggsave("plots/world_map_samples.png",
       plot = p_world_arrows, 
       device = "png", 
       width = 6, height = 6, units = "in")
ggsave("../../bee_manuscript/figures/world_map_samples.pdf",
       plot = p_world_arrows, 
       device = "pdf", 
       width = 6, height = 6, units = "in")

# add dates to world map
dates_AHB <- read.table("../maps/spread_africanized_bees.csv",
                        sep = ",", header = T)
omit_dates <- c("Brazil", "Venezuela", "Uruguay", "Costa Rica")
left_dates <- c("Mexico", "Peru", "Bolivia", "Paraguay", "Panama")
right_dates <- c("California", "Texas", "Guyana", "Colombia", "Argentina")

bind_rows(data.frame(name = "Brazil", 
                     date = 1957,
                     location = "Rio Claro, Sao Paulo",
                     lat = -22.4149, lon = -47.5651,
                     source = "Kent 1988 The introduction and diffusion of the African honeybee in South America",
                     stringsAsFactors = F),
                     filter(dates_AHB, name %in% c(left_dates, right_dates))) %>%
  dplyr::select(., date, name, location, lat, lon, source) %>%
  arrange(date) %>%
  mutate(source = ifelse(source == "Kent 1988 The introduction and diffusion of the African honeybee in South America",
                         "Kent 1988, The introduction and diffusion of the African honeybee in South America",
                         source)) %>%
  write.table("../../bee_manuscript/files_supp/S6_Table_Africanized_honey_bee_invasion_historical_dates_and_locations.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")

p_world_labels <- p_world_arrows +
  geom_text(data = mutate(sao_paulo, label = "Brazil\n1957"), aes(x = long, y = lat, label = label), size = 4, 
            hjust = 0.2, vjust = -0.3, 
            color = col_ACM["A"]) +
  geom_text(data = filter(dates_AHB, name %in% right_dates), aes(x = lon, y = lat, label = paste(name, date)), size = 3, 
            hjust = -0.2, vjust = 0.5) +
  geom_text(data = filter(dates_AHB, name %in% left_dates), aes(x = lon, y = lat, label = paste(name, date)), size = 3, 
            hjust = 1.2, vjust = 0.5) +
  geom_point(data = dates_AHB, aes(x = lon, y = lat), size = 0.5)
p_world_labels
ggsave("plots/world_map_dates.png",
       plot = p_world_labels, 
       device = "png", 
       width = 6, height = 5, dpi = 600, units = "in")
ggsave("../../bee_manuscript/figures/world_map_dates.png",
       plot = p_world_labels, 
       device = "png", 
       width = 6, height = 5, dpi = 600, units = "in")

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
               alpha = .3,
               fill = dark2[8],
               color = "white"
  ) +
  xlim(c(-140, -20)) +
  theme_classic() +
  geom_point(aes(x = long,
                 y = lat,
                 color = population,
                 shape = factor(shape)),
             data = filter(admix.pops, continent == "S. America"),
             size = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis_d(option = "viridis", direction = 1,
                        begin = 0, end = .8) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  theme(legend.position = "None")

dist_100km_lat <- 100000/111319.5 # in degrees latitude

SA_scale <- data.frame(long = c(-55.5 - 0.5, -55.5 - 0.5),
                       lat = c(-35.5, -35.5-dist_100km_lat))

distm(c(SA_scale$long[1], SA_scale$lat[1]), c(SA_scale$long[2], SA_scale$lat[2]), fun = distHaversine)

SA_points_zoom <- SA_points + 
  coord_quickmap(# approximate mercator projection where lat = long distances (ok for zoomed in maps)
    xlim = c(-62, -55),
    ylim = c(-37, -27)) +
  geom_line(data = SA_scale, aes(x = long, y = lat), color = "black") +
  geom_text(data = data.frame(lat = rep(mean(SA_scale$lat), 2), 
                              long = SA_scale$long + c(-.4, 0.4), dist = c("100", "km")),
             aes(x = long, y = lat, label = dist), size = 2, angle = 90)
SA_points_zoom

ggsave("plots/SA_point_map_samples_zoom_out.png",
       plot = SA_points_zoom +
         ggtitle("Argentina"), 
       device = "png", 
       width = 3, height = 6, units = "in")
ggsave("../../bee_manuscript/figures/SA_point_map_samples_zoom_out.png",
       plot = SA_points_zoom, 
       device = "png", 
       width = 3, height = 6, units = "in", dpi = 600)


# NA just colored points, not pie charts:
NA_points <- world %>%
  ggplot(data = .) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               alpha = .3,
               color = "white",
               fill = dark2[8]
  ) +
  geom_polygon(aes(x = long, y = lat, group = group),
              alpha = 0, 
              color = "white",
              fill = dark2[8],
              data = filter(states, region %in% c("nevada", "arizona"))) +
  geom_point(aes(x = long,
                 y = lat,
                 color = population_factor,
                 shape = factor(shape)), 
             data = filter(admix.pops, continent == "N. America"),
             cex = 2) +
  xlim(c(-140, -20)) +
  theme_classic() +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis_d(option = "plasma", direction = -1,
                        begin = 0.1, end = 0.85) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  theme(legend.position = "None")

# smaller points
NA_scale <- data.frame(long = c(-122.5, -122.5),
                       lat = c(33.5, 33.5-dist_100km_lat))
NA_points_zoom <- NA_points +
  #coord_fixed(1.1,
  coord_quickmap(
              xlim = c(-124, -114.5), 
              ylim = c(32, 39.5)) +
  geom_line(data = NA_scale, aes(x = long, y = lat), color = "black") +

  geom_text(data = data.frame(lat = rep(mean(NA_scale$lat), 2), 
                              long = NA_scale$long + c(-.4, 0.4), dist = c("100", "km")),
            aes(x = long, y = lat, label = dist), size = 2, angle = 90)
NA_points_zoom
ggsave("plots/NA_point_map_samples_zoom_out.png",
       NA_points_zoom +
         ggtitle("California"), 
       device = "png", 
       width = 3, height = 6, units = "in")
ggsave("../../bee_manuscript/figures/NA_point_map_samples_zoom_out.png",
       NA_points_zoom, 
       device = "png", 
       width = 3, height = 6, units = "in", dpi = 600)

# make a scale bar with the symbols for my structure-like plot:
CA_bar <- admix.pops %>%
  filter(continent == "N. America") %>%
  arrange(desc(lat)) %>%
  mutate(end = cumsum(n)) %>%
  mutate(start = end - n) %>%
  mutate(center = start + n/2)
CA_bar_symbols <- ggplot() +
  geom_polygon(aes(x = X, y = Y, group = Group),
               alpha = 0.3,
               fill = dark2[8],
               data = data.frame(X = c(-1, 1, 1, -1),
                                 Y = c(-1, -1, 150, 150),
                                 Group = rep(1, 4))) +
  geom_point(data = CA_bar,
             aes(x = 0, y = center, 
                 color = population_factor,
                 shape = factor(shape)), 
             size = 3) +
  scale_color_viridis_d(option = "plasma", direction = -1,
                        begin = 0.1, end = 0.85) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  #theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "None") +
  xlim(-.1, .1) +
  ylim(0, 135) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_hline(data = CA_bar, aes(yintercept = start),
             color = "white") +
  geom_hline(data = CA_bar, aes(yintercept = end),
             color = "white")


# make the same symbol set for AR Argentina samples:
AR_bar <- admix.pops %>%
  filter(continent == "S. America") %>%
  arrange(lat) %>%
  mutate(end = cumsum(n)) %>%
  mutate(start = end - n) %>%
  mutate(center = start + n/2)
AR_bar_symbols <- ggplot() +
  geom_polygon(aes(x = X, y = Y, group = Group),
               alpha = 0.3,
               fill = dark2[8],
               data = data.frame(X = c(-1, 1, 1, -1),
                                 Y = c(-1, -1, 150, 150),
                                 Group = rep(1, 4))) +
  geom_point(data = AR_bar,
             aes(x = 0, y = center, 
                 color = population_factor,
                 shape = factor(shape)), 
             size = 3) +
  scale_color_viridis_d(option = "viridis", direction = 1,
                        begin = 0, end = .8) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  #theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "None") +
  xlim(-.1, .1) +
  ylim(0, 178) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_hline(data = AR_bar, aes(yintercept = start),
             color = "white") +
  geom_hline(data = AR_bar, aes(yintercept = end),
             color = "white")


# putting plots together
p_world_together <- p_world_labels +
  annotation_custom(grob = ggplotGrob(NA_points_zoom +
                                        ggtitle("California") +
                                        theme( # get rid of axes
                                          axis.line = element_blank(),
                                          axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.title = element_blank(),
                                          panel.background = element_rect(fill = "transparent", colour = NA),
                                          plot.margin = unit(c(0,0,0,0), "null"),
                                          #panel.spacing = unit(c(0,0,0,0), "null"),
                                          panel.border = element_rect(colour = col_NA_SA_both["N. America"], 
                                                                      fill = NA, 
                                                                      size = 2),
                                          plot.title = element_text(hjust = 0.5, size = 10,
                                                                    margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = "in")))), 
                    xmin = -50,
                    xmax = -20,
                    ymin = 10) +
  annotation_custom(grob = ggplotGrob(SA_points_zoom +
                                        
                                        ggtitle("Argentina") +
                                        theme( # get rid of axes
                                          axis.line = element_blank(),
                                          axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.title = element_blank(),
                                          panel.background = element_rect(fill = "transparent", colour = NA),
                                          panel.border = element_rect(colour = col_NA_SA_both["S. America"], 
                                                                      fill = NA, 
                                                                      size = 2),
                                          plot.title = element_text(hjust = 0.5, size = 10,
                                                                    margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = "in")))), 
                    xmin = -125, 
                    xmax = -90, 
                    ymax = 0)
p_world_together

# redo components so they fit together better in the plot:
NA_plot <- d_admix %>% 
  left_join(., meta.pop[ , c("population", "lat")], by = "population") %>%
  arrange(desc(lat.y), desc(lat.x)) %>% # lat within populations, but then group by pop (b/c some 2014 and 2018 pops overlap slightly in CA)
  mutate(order = 1:nrow(.)) %>%
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "CA_2018" | population %in% pop2014_incl) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, alpha = as.factor(year), y=p, x=reorder(label, order))) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  scale_alpha_discrete(range = c(0.6, 1)) +
  labs(alpha = "Collection") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") + 
  theme_classic() +
  scale_y_reverse(name = "California", breaks = c(0, 0.5, 1), labels = c("1", "0.5", "0")) +
  ggtitle("<- Brazil") +
  theme( # get rid of axes
    plot.title = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  guides(fill = "none",
         alpha = "none")
NA_plot
NA_plot_byindlat <- d_admix %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "CA_2018" | population %in% pop2014_incl) %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill=ancestry_label, alpha = as.factor(year), y=p, x=reorder(label, -lat))) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  scale_alpha_discrete(range = c(0.6, 1)) +
  labs(alpha = "Collection") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") + 
  theme_classic() +
  scale_y_reverse(name = "California", breaks = c(0, 0.5, 1), labels = c("1", "0.5", "0")) +
  ggtitle("<- Brazil") +
  theme( # get rid of axes
    plot.title = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  guides(fill = "none",
         alpha = "none")
NA_w_shapes <- NA_plot +
  geom_point(data = CA_bar,
             aes(x = center + 0.5, 
                 y = 1.1, 
                 color = population_factor,
                 shape = factor(shape)), 
             fill = "white",
             size = 1,
             alpha = 1) +
  scale_color_viridis_d(option = "plasma", direction = -1,
                        begin = 0.1, end = 0.85) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  guides(color = "none", shape = "none") +
  geom_segment(data = CA_bar %>% mutate(ancestry_label = "A") %>% bind_rows(., data.frame(end = 0)), 
               aes(x = end + 0.5, xend = end + 0.5, 
                                  y = 1, yend = 1.15,
                   lwd = 0.1),
               color = "black", size = 0.1, alpha = 1)

NA_w_shapes
SA_plot <- d_admix %>% tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  filter(group == "AR_2018") %>%
  left_join(., anc, by = "ancestry") %>%
  ggplot(., aes(fill = ancestry_label, y = p, x = reorder(Bee_ID, lat))) +
  geom_bar(stat = "identity", position = "fill", width = 1) + 
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  theme_classic() + 
  scale_y_reverse(name = "Argentina", breaks = c(0, 0.5, 1), labels = c("1", "0.5", "0")) +
  theme( # get rid of axes
    plot.title = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  guides(color = "none", fill = "none")
SA_w_shapes <- SA_plot + 
  geom_point(data = AR_bar,
             aes(x = center + 0.5, 
                 y = 1.1, 
                 color = population_factor,
                 shape = factor(shape)), 
             fill = "white",
             size = 1,
             alpha = 1) +
  scale_color_viridis_d(option = "viridis", direction = 1,
                        begin = 0, end = .8) +
  scale_shape_manual(values = c(15, 1, 17, 5, 19, 6)) +
  guides(color = "none", shape = "none") +
  geom_segment(data = AR_bar %>% mutate(ancestry_label = "A") %>% bind_rows(., data.frame(end = 0)), 
               aes(x = end + 0.5, xend = end + 0.5, 
                   y = 1, yend = 1.15,
                   lwd = 0.1),
               color = "black", size = 0.1, alpha = 1)
AMC_legend = ggplot(data = data.frame(Ancestry = factor(ACM, levels = c("A", "M", "C"), ordered = T), x = 1:3, y = 1:3), aes(x = x, y = y, color = Ancestry)) +
  geom_point() +
  scale_color_manual(values = col_ACM, 
                     name = element_blank()) +
                     #name = expression(""*symbol('\254')* " Brazil")) +
  theme_classic() +
  theme(legend.key = element_rect(size = 0.01, color = "white"),
        legend.key.size = unit(0.5, units = "cm"),
        legend.spacing = unit(0, units = "cm"))
p_world_admix_tall <- grid.arrange(grobs = list(ggplotGrob(p_world_together),
                                                cowplot::get_legend(AMC_legend),                                                              ggplotGrob(NA_w_shapes),
                                                ggplotGrob(SA_w_shapes)),
                                   layout_matrix = rbind(c(1,1),
                                                         c(NA, 2),
                                                         c(3,2),
                                                         c(4,4)),
                                   heights = c(7, 0.25, 1.1, 1.1),
                                   widths = c(8, 1))


plot(p_world_admix_tall)
ggsave("plots/world_map_ngsadmix_tall.png",
       plot = p_world_admix_tall, 
       device = "png", 
       width = 7.5, height = 6.75, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/world_map_ngsadmix_tall.png",
       plot = p_world_admix_tall, 
       device = "png", 
       width = 7.5, height = 6.75, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/Fig1.tiff",
       plot = p_world_admix_tall, 
       device = "tiff", 
       width = 7.5, height = 6.75, units = "in", dpi = 600,
       compression = "lzw", type = "cairo")
ggsave("../../bee_manuscript/figures_main/world_map_ngsadmix_tall_300dpi.tiff",
       plot = p_world_admix_tall, 
       device = "tiff", 
       width = 7.5, height = 6.75, units = "in", dpi = 300,
       compression = "lzw", type = "cairo")



# ----------------------------------------------
# print file with population frequencies of each admixed population 
# to use as priors in local ancestry inference

admix.pops %>%
  dplyr::select(c("population", "A", "C", "M", "n")) %>%
  write.table(., paste0("results/NGSAdmix/", name, ".pop.anc"),
              quote = F, col.names = T, row.names = F, sep = "\t")
# write out individual ancestry too
admix.ind %>%
  dplyr::select(Bee_ID, population, A, C, M) %>%
  write.table(., paste0("results/NGSAdmix/", name, ".ind.anc"),
              quote = F, col.names = T, row.names = F, sep = "\t")



#---------------------------------------------------------------------------------------------
# Comparison of ancestry_hmm and NGSAdmix global ancestry estimates:
# first read in data from ancestry_hmm, summarised by individual:
# first read in individual alpha estimates for mean A ancestry
indAlpha <- do.call(rbind,
                    lapply(meta.pop$population, function(p)
                      read.table(paste0("../local_ancestry/results/ancestry_hmm/combined_sept19/posterior/anc/", p, ".alpha.anc"),
                                 stringsAsFactors = F, header = T)))


# compare the priors with the mean inferred ancestry from ancestry_hmm
# make individual ancestry plots too:
compare_anc_ind <- indAlpha %>% # ancestry_hmm genomewide estimates for individuals
  tidyr::gather(., "ancestry", "ancestry_hmm", c("A", "C", "M")) %>%
  left_join(., tidyr::gather(admix.ind[ , c("Bee_ID", "A", "C", "M")], "ancestry", "NGSAdmix", c("A", "C", "M")), # NGSAdmix genomewide estimates for individuals
            by = c("ancestry"="ancestry", "ID"="Bee_ID"))
p_ind <- compare_anc_ind %>%
  ggplot(., aes(NGSAdmix, ancestry_hmm, color = ancestry)) +
  geom_point(alpha = .5, pch = 1) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
  xlab("NGSAdmix") +
  ylab("ancestry_hmm") +
  scale_color_manual(values = col_ACM, name = "Ancestry") +
  coord_fixed()
p_pop <- compare_anc_ind %>%
  left_join(., meta.ind, by = c("ID"="Bee_ID")) %>%
  group_by(population, ancestry) %>%
  summarise(NGSAdmix = mean(NGSAdmix),
            ancestry_hmm = mean(ancestry_hmm)) %>%
  ggplot(., aes(NGSAdmix, ancestry_hmm, color = ancestry)) +
  geom_point(alpha = .75) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
  xlab("NGSAdmix") +
  ylab("ancestry_hmm") +
  scale_color_manual(values = col_ACM, name = "Ancestry") +
  coord_fixed()
p_pop
p_pop_ind <- grid.arrange(p_ind + ggtitle("A") + theme(legend.position = "none"), 
                          p_pop + ggtitle("B") + theme(legend.position = "none"), 
                          get_legend(p_pop)$grobs[[1]],
                          widths = c(5, 5, 2),
                          nrow = 1, ncol = 3, right = 2)
p_pop_ind
# save plot locally and in bee_manuscript figures folder.
ggsave("plots/mean_ancestry_prior_posterior_ancestry_hmm.png",
       plot = p_pop_ind,
       device = "png",
       width = 7.5, height = 4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/mean_ancestry_prior_posterior_ancestry_hmm.png",
       plot = p_pop_ind,
       device = "png",
       width = 7.5, height = 4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/mean_ancestry_prior_posterior_ancestry_hmm.tiff",
       plot = p_pop_ind,
       device = "tiff",
       width = 7.5, height = 4, units = "in", dpi = 600)
compare_anc_ind %>%
  group_by(ancestry) %>%
  summarise(corr = cor(ancestry_hmm, NGSAdmix, method = "pearson"))
compare_anc_ind %>%
  group_by(ancestry) %>%
  summarise(corr = cor(ancestry_hmm, NGSAdmix, method = "pearson"))
