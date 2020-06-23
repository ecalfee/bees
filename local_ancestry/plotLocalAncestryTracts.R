# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(reshape2)
library(viridis)
library(reshape2) # melt function
library(LaplacesDemon)
library(emdbook)
library(betareg)
library(grid)
library(gridExtra)
library(MASS) # for mvrnorm
library(plotly)
library(ggpointdensity)
library(ggExtra)
library(gtable)
library(egg) # for plot layouts
library(coda) # for hpdi calc
library(sp) # for polygons
library(GISTools)
#library(hex) # for hex plot
#library(ggridges) # to get density plot without bottom line on x axis
source("../colors.R") # for color palette
#source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions
#source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions
source("k_matrix.R") # made its own local copy for this git repository
source("calc_FDRs.R") # scripts to calculate false discovery rates

#************************************************************************************************************#

# Newest bee sequences
# load metadata
# load population ancestry frequencies:
#pops <- read.table("../bee_samples_listed/byPop/pops_included.list", stringsAsFactors = F)$V1
pops <- read.table("../bee_samples_listed/byPop/combined_sept19_pops.list", stringsAsFactors = F)$V1
bees <- do.call(rbind, 
                lapply(pops, function(p) data.frame(Bee_ID = read.table(paste0("../bee_samples_listed/byPop/", p, ".list"),
                                                    stringsAsFactors = F)$V1, population = p, stringsAsFactors = F)))

# get metadata
meta.ind0 <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  left_join(bees, ., by = c("Bee_ID", "population"))
dates <- rbind(filter(meta.ind0, geographic_location == "Argentina") %>%
  mutate(date = as.Date(Date, "%d.%m.%y")) %>%
  dplyr::select(c("Bee_ID", "Date", "date")),
  filter(meta.ind0, geographic_location == "California") %>%
    mutate(date = as.Date(Date, "%m/%d/%y")) %>%
    dplyr::select(c("Bee_ID", "Date", "date")))
meta.ind <- left_join(meta.ind0, dates, by = c("Bee_ID", "Date")) %>%
  dplyr::select(-Date) %>%
  rename(time = Time)
meta.pop <- meta.ind %>%
  dplyr::select(c("population", "source", "year", "group", "lat", "long")) %>%
  dplyr::group_by(population, source, year, group) %>%
  dplyr::summarise(n_bees = n(),
                   lat = mean(lat),
                   long = mean(long)) %>%
  left_join(data.frame(population = pops, stringsAsFactors = F), # limits to included populations
            ., by = "population") %>%
  # not relevant anymore, not using Riverside 1999 historical samples:
  mutate(lat = ifelse(population == "Riverside_1999", .[.$population == "Riverside_2014", "lat"], lat)) %>%
  mutate(zone = ifelse(group == "AR_2018", "S. America", "N. America")) %>%
  arrange(lat) %>%
  mutate(region = ifelse(zone == "S. America",
                         ifelse(lat < -32.26, "Low A", "High A"),
                         ifelse(lat > 32.72, "Low A", "High A")))

# included populations by latitude:
pops_by_lat <- meta.pop$population[order(meta.pop$lat)]
meta.AR.order.by.lat <- data.frame(population = pops_by_lat, stringsAsFactors = F) %>%
  left_join(., meta.pop, by = "population") %>%
  filter(zone == "S. America") %>%
  mutate(abs_lat_SA_c = abs(lat) - mean(abs(lat))) # absolute latitude centered for SA

# regions defined based on cline center for 'low' and 'high' sides of the cline
AR_pops_S <- meta.pop$population[meta.pop$zone == "S. America" & meta.pop$region == "Low A"]
AR_pops_N <- meta.pop$population[meta.pop$zone == "S. America" & meta.pop$region == "High A"]
CA_pops <- meta.pop$population[meta.pop$zone == "N. America"]
# now with legend for 'low A' and 'high A'
NS_segments = data.frame(Region = c("Low A", "High A", "Low A"),
                         starts = c(0.5, 
                                    length(AR_pops_S) + 0.5, 
                                    length(c(AR_pops_S, AR_pops_N)) + 0.5),
                         ends = c(length(AR_pops_S) + 0.5, 
                                  length(c(AR_pops_S, AR_pops_N)) + 0.5, 
                                  length(c(AR_pops_S, AR_pops_N, CA_pops)) + 0.5))

save(file = "results/meta.RData", list = c("meta.ind", "meta.pop", "pops_by_lat", "meta.AR.order.by.lat", "NS_segments"))

# get SNP sites where ancestry was called
sites0 <- read.table("results/SNPs/combined_sept19/chr.var.sites", stringsAsFactors = F,
                     sep = "\t", header = F)[ , 1:2]
colnames(sites0) <- c("scaffold", "pos")
chr_lengths <- cbind(read.table("../data/honeybee_genome/chr.names", stringsAsFactors = F),
                     read.table("../data/honeybee_genome/chr.lengths", stringsAsFactors = F)) %>%
  data.table::setnames(c("chr", "scaffold", "chr_length")) %>%
  mutate(chr_n = 1:16) %>%
  mutate(chr_end = cumsum(chr_length)) %>%
  mutate(chr_start = chr_end - chr_length) %>%
  mutate(chr_mid = (chr_start + chr_end)/2)

sites <- left_join(sites0, chr_lengths[ , c("chr", "scaffold", "chr_n", "chr_start")], by = "scaffold") %>%
  mutate(cum_pos = pos + chr_start)

# get ancestry frequencies for each population across the genome
dir_results <- "results/ancestry_hmm/combined_sept19/posterior"

popA <- lapply(pops_by_lat, function(p) read.table(paste0(dir_results, "/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A <- do.call(cbind, popA)
rm(popA)
colnames(A) <- pops_by_lat

popM <- lapply(pops_by_lat, function(p) read.table(paste0(dir_results, "/anc/", p, ".M.anc"),
                                            stringsAsFactors = F))
M <- do.call(cbind, popM)
colnames(M) <- pops_by_lat
rm(popM)
popC <- lapply(pops_by_lat, function(p) read.table(paste0(dir_results, "/anc/", p, ".C.anc"),
                                            stringsAsFactors = F))
C <- do.call(cbind, popC)
colnames(C) <- pops_by_lat
rm(popC)

# mean ancestry across populations
meanA0 <- apply(A, 1, mean)
meanC0 <- apply(C, 1, mean)
meanM0 <- apply(M, 1, mean)

# mean ancestry across the genome for each pop
admix_proportions <- data.frame(population = pops_by_lat,
                                A = apply(A, 2, mean),
                                C = apply(C, 2, mean),
                                M = apply(M, 2, mean),
                                stringsAsFactors = F)


# get time of admixture estimates
time_pops <- read.table(paste0(dir_results, "/", "time_pops.txt"), 
                        stringsAsFactors = F, header = F) %>%
  data.table::setnames("population")
admix_times <- do.call(rbind, lapply(c("A", "C", "M"), function(anc) read.table(paste0(dir_results, "/", "time_", anc, ".txt"), 
                     stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("ancestry_n", "time", "proportion")) %>%
    mutate(ancestry = anc) %>% 
    cbind(time_pops, .)))

# the ancestry proportion values in the log rom ancestry_hmm are just the priors from NGSAdmix:
# admix.pops and admix.ind are loaded from plot_NGSadmix.R script
left_join(admix.pops, filter(admix_times, ancestry == "A"), by = "population") %>%
  with(., plot(A, proportion, xlim = 0:1, ylim = 0:1, col = "blue", main = "should be perfect match -- showing prior"))
abline(0, 1)


# plot ancestry times vs. mean ancestry proportion:
admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  filter(., ancestry != "C") %>% # no time, first ancestry
  ggplot(., aes(x = abs(lat), y = time, color = zone, shape = factor(year))) +
  geom_point(alpha = .75) +
  facet_grid(. ~ ancestry) +
  scale_color_manual(values = col_NA_SA_both, name = "Continent") +
  #ggtitle("Inferred time of admixture pulses from HMM") +
  ylab("Time (generations)") +
  xlab("Degrees latitude from the equator") +
  labs(shape = "Sample") +
  theme_classic() +
  scale_shape_manual(values = c(17, 19))
ggsave("plots/time_of_admixture_vs_latitude.png",
       height = 3, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures/time_of_admixture_vs_latitude.png",
       height = 3, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/time_of_admixture_vs_latitude.tiff",
       height = 3, width = 5.2, units = "in", dpi = 600)
summary(filter(admix_times, ancestry == "A")$time)
hist(filter(admix_times, ancestry == "A")$time)
filter(admix_times, ancestry == "A") %>%
  arrange(time) %>%
  View(.)
admix_times %>%
  group_by(ancestry) %>%
  summarise(mean = mean(time),
            median = median(time),
            min = min(time),
            max = max(time))
range(filter(admix_times, ancestry == "A")$time)

# compare California and Argentina: 
# For the same ancestry proportions, do they have similar inferred times of admixture?
# This would be evidence for different amounts of ongoing gene flow in the two zones.
admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  #filter(., ancestry == "A") %>%
  filter(., ancestry != "C") %>%
  mutate(ancestry_labels = paste(ancestry, "ancestry admixture pulse")) %>%
  ggplot(., aes(x = proportion, y = time, color = zone, shape = factor(year))) +
  geom_point(alpha = .75) +
  #ggtitle("Time since admixture estimated from ancestry block lengths") +
  facet_grid(. ~ ancestry) +
  scale_color_manual(values = col_NA_SA_both, name = "Hybrid zone") +
  ylab("Time (generations in the past)") +
  xlab("Admixture proportion") +
  labs(color = "Hybrid Zone", shape = "Year") +
  theme_classic() +
  scale_shape_manual(values = c(17, 19))
ggsave("plots/California_has_shorter_ancestry_blocks.png",
       height = 3, width = 6, units = "in")


#CA_pops_included <- meta.pop[meta.pop$group %in% c("N_CA", "S_CA", "CA_2018") & meta.pop$year >= 2014, ]
CA_pops_included <- meta.pop[meta.pop$zone == "N. America", ]
CA_A <- A[, CA_pops_included$population]
CA_A_earlier <- A[, meta.pop$population[meta.pop$group %in% c("N_CA", "S_CA") & meta.pop$year < 2014]]
CA_A_2014 <- A[, meta.pop$population[meta.pop$group %in% c("N_CA", "S_CA") & meta.pop$year == 2014]]
CA_A_2018 <- A[, meta.pop$population[meta.pop$group %in% c("CA_2018")]]
plot(apply(CA_A_2014, 1, mean), apply(CA_A_2018, 1, mean))
plot(apply(CA_A_earlier, 1, mean), apply(CA_A_2018, 1, mean))
AR_pops_included <- meta.pop[meta.pop$zone == "S. America", ]
AR_A <- A[, AR_pops_included$population]

# take mean across individuals, not across populations:
meanA_CA0 <- apply(CA_A, 1, mean)
meanA_CA <- apply(CA_A, 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanA_AR0 <- apply(AR_A, 1, mean)
meanA_AR <- apply(AR_A, 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))

# with only included pops, mean of all individuals:
meanA <- apply(A, 1, function(x) sum(x*meta.pop$n_bees)/sum(meta.pop$n_bees))

# C and M ancestries:
meanC_CA <- apply(C[ , CA_pops_included$population], 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanC_AR <- apply(C[ , AR_pops_included$population], 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanC <- apply(C, 1, function(x) sum(x*meta.pop$n_bees)/sum(meta.pop$n_bees))
meanM_CA <- apply(M[ , CA_pops_included$population], 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanM_AR <- apply(M[ , AR_pops_included$population], 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanM <- apply(M, 1, function(x) sum(x*meta.pop$n_bees)/sum(meta.pop$n_bees))

# save ancestry data for later access:
save(file = "results/A.RData", list = c("A", "sites", "meanA", "meanA_CA", "meanA_AR"))
save(file = "results/C.RData", list = c("C", "sites", "meanC", "meanC_CA", "meanC_AR"))
save(file = "results/M.RData", list = c("M", "sites", "meanM", "meanM_CA", "meanM_AR"))
#load("results/A.RData")
#load("results/C.RData")
#load("results/M.RData")

# make plots - all loci, mean ancestry across all pops
png("plots/histogram_A_ancestry_all_loci.png", height = 6, width = 8, units = "in", res = 300)
hist(meanA, main = "all loci: mean ancestry across populations")
abline(v = max(meanA), col = "blue")
abline(v = min(meanA), col = "blue")
abline(v = quantile(meanA, .99), col = "orange")
abline(v = quantile(meanA, .01), col = "orange")
dev.off()

sites %>%
  filter(meanA_CA > quantile(meanA_CA, .99) & 
            meanA_AR > quantile(meanA_AR, .99)) %>%
  dplyr::select(scaffold) %>%
  table() # outliers are on just a few scaffolds
# plot
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter(meanA_CA > quantile(meanA_CA, .99) & 
           meanA_AR > quantile(meanA_AR, .99)) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) # really just a few peaks
sites %>%
  filter(meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  dplyr::select(scaffold) %>%
  table() # low A outliers still only on 14 scaffolds
#plot
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter(meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) # really just a few peaks
sites %>%
  filter(meanA_CA > quantile(meanA_CA, .25) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  dplyr::select(scaffold) %>%
  table() # several regions, but at least one big hit on chr 11 - Group11.18
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>% # plot
  filter(meanA_CA > quantile(meanA_CA, .25) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) # several regions, but at least one big hit on chr 11 - Group11.18
sites %>%
  filter(meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR > quantile(meanA_AR, .25)) %>%
  dplyr::select(scaffold) %>%
  table() # more spread out but a few scaffolds with higher concentrations
sites %>%
  filter(meanA_CA < quantile(meanA_CA, .75) & 
           meanA_AR > quantile(meanA_AR, .99)) %>%
  dplyr::select(scaffold) %>%
  table()
sites %>%
  filter(meanA_CA > quantile(meanA_CA, .99) & 
           meanA_AR < quantile(meanA_AR, .75)) %>%
  dplyr::select(scaffold) %>%
  table()
# plot shared high and shared low on same axes:
#plot
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter((meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR < quantile(meanA_AR, .01)) |
           (meanA_CA > quantile(meanA_CA, .99) & 
              meanA_AR > quantile(meanA_AR, .99))) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) + # really just a few peaks
  ggtitle("shared peaks high and low A ancestry")


# plot all data. oops I need to use position relative to the chromosome instead.
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>% # plot
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point(size = .1) +
  facet_wrap(~chr, scales = "free_x")

png("plots/scatterplot_CA_vs_AR_A_ancestry_all_loci.png", height = 6, width = 8, units = "in", res = 300)
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(meanA_CA, .99) & 
                                                   meanA_AR > quantile(meanA_AR, .99)) |
                                                  (meanA_CA < quantile(meanA_CA, .01) & 
                                                     meanA_AR < quantile(meanA_AR, .01)), 
          alpha("orange", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, orange = top 1% shared outliers")
abline(lm(meanA_AR~meanA_CA), col = "blue")
dev.off()
cor(meanA_CA, meanA_AR)

# color outliers by scaffold/chromosome:
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>% # plot
  ggplot(aes(x = CA, y = AR, color = scaffold_n)) +
  geom_point(size = .1) +
  facet_wrap(~chr) + 
  ggtitle("African ancestry frequencies in CA vs. AR, all loci")
ggsave("plots/mean_A_ancestry_CA_vs_AR_rainbow_by_chr.png", 
       device = "png", height = 10, width = 12, units = "in")
# zoom in on Group1.23 because it has both high and low outliers:
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  #filter(chr == "Group1" & scaffold_n %in% 20:23) %>%
  filter(chr == "Group1") %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point(size = .1) +
  #facet_wrap(~scaffold_n) + 
  ggtitle("African ancestry frequencies in CA and AR")
ggsave("plots/mean_A_ancestry_CA_and_AR_Group1_scaffolds_20-23.png", 
       device = "png", height = 10, width = 12, units = "in")


# plot K matrix (! order by latitude!)
zAnc_bees = make_K_calcs(t(A[ , pops_by_lat]))

save(file = "results/zAnc.RData", list = "zAnc_bees")
#load("results/zAnc.RData")
# what are mean correlations?
# within CA
# within AR S
# within AR N
# CA -> AR S
# CA -> AR N
# AR S -> AR N


mean_corr_k <- get_mean_from_K(cov2cor(zAnc_bees$K))
mean_corr_k
mean_corr_k %>%
  write.table(., "results/mean_anc_corr_grouped.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
mean_cov_k <- get_mean_from_K(zAnc_bees$K)
mean_cov_k
mean_cov_k %>%
  write.table(., "results/mean_anc_cov_grouped.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
# plot the lower triangle of K
lower_tri <- zAnc_bees$K
lower_tri[lower.tri(lower_tri, diag = T)] <- NA
melt(lower_tri) %>% 
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

k_plot_all <- melt(cov2cor(zAnc_bees$K)) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     #limits = c(-.12, .32),
                     name = "Ancestry\ncorrelation") +
  ggtitle("K correlation matrix - bees") +
  xlab("") +
  ylab("")

k_plot_all

segment_buffer = 0.15
k_plot_fancy <- k_plot_all +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments,
               aes(x = starts + segment_buffer, xend = ends - segment_buffer,
               y = -0.5, yend = -0.5,
               color = Region), size = 2, inherit.aes = F) +
  geom_segment(data = NS_segments,
               aes(x = -0.25, xend = -0.25,
                   y = starts + segment_buffer, yend = ends - segment_buffer,
                   color = Region), size = 2, inherit.aes = F) +
  #add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5)) +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5)) +
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), 
                   #expand = expand_scale(add = 2),
                   labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"),
                   #expand = expand_scale(add = 2),
                   labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A)

k_plot_fancy
ggsave("plots/k_correlation_matrix_all_pops.png", 
       plot = k_plot_fancy,
       height = 3, width = 4, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures/k_correlation_matrix_all_pops.png", 
       plot = k_plot_fancy,
       height = 3, width = 4, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures_main/k_correlation_matrix_all_pops.tif", 
       plot = k_plot_fancy,
       height = 3, width = 4, 
       units = "in", device = "tiff", dpi = 600,
       compression = "lzw", type = "cairo")

# make covariance matrix for the supplement:
k_cov_plot_offdiag <- melt(zAnc_bees$K) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     option = "viridis",
                     name = "Ancestry\ncovariance") +
  ggtitle("K covariance matrix - bees") +
  xlab("") +
  ylab("") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments,
               aes(x = starts + segment_buffer, xend = ends - segment_buffer,
                   y = -0.5, yend = -0.5,
                   color = Region), size = 2, inherit.aes = F) +
  geom_segment(data = NS_segments,
               aes(x = -0.25, xend = -0.25,
                   y = starts + segment_buffer, yend = ends - segment_buffer,
                   color = Region), size = 2, inherit.aes = F) +
  # add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5), color = "white") +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5), color = "white") +
  #geom_text(colour = "darkgray", aes(y = -3, label = zone1),  position = position_dodge(width=0.9))
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A)

k_cov_plot_offdiag

k_cov_plot_diag <- melt(zAnc_bees$K) %>%
  #filter(Var1 == Var2) %>% # diagonal only
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     option = "inferno",
                     discrete = F,
                     name = "Ancestry\ncovariance") +
  
  ggtitle("K covariance matrix - bees") +
  xlab("") +
  ylab("") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments,
               aes(x = starts + segment_buffer, xend = ends - segment_buffer,
                   y = -0.5, yend = -0.5,
                   color = Region), size = 2, inherit.aes = F) +
  geom_segment(data = NS_segments,
               aes(x = -0.25, xend = -0.25,
                   y = starts + segment_buffer, yend = ends - segment_buffer,
                   color = Region), size = 2, inherit.aes = F) +
  # add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5), color = "white") +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5), color = "white") +
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A) #+
  #guides(fill = guide_legend(order = 2), color = guide_legend(order = 1))

k_cov_plot_diag
k_cov_plot_fancy <- arrangeGrob(grobs = list(k_cov_plot_diag + guides(color = F),
                                             k_cov_plot_offdiag + guides(color = F),
                                             ggpubr::get_legend(k_cov_plot_diag + guides(fill = F) + 
                                                                  theme(legend.position = "bottom"))),
                                ncol = 2, nrow = 2, 
                                layout_matrix = rbind(c(1,2),c(3,3)),
                                heights = c(10,1), widths = c(1,1))
plot(k_cov_plot_fancy)

ggsave("plots/k_covariance_matrix_all_pops.png", 
       plot = k_cov_plot_fancy,
       height = 3, width = 7.5, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures/k_covariance_matrix_all_pops.png", 
       plot = k_cov_plot_fancy,
       height = 3, width = 7.5, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/k_covariance_matrix_all_pops.tif", 
       plot = k_cov_plot_fancy,
       height = 3, width = 7.5, 
       units = "in", device = "tiff", dpi = 600,
       compression = "lzw", type = "cairo")

# look at covariance within European ancestry
# proportion M conditional on European ancestry
M_within_Euro = M/(M + C) # proportion European ancestry that's M
zAnc_M_Euro = make_K_calcs(t(M_within_Euro[ , pops_by_lat]))
k_M_plot_all <- #melt(cov2cor(zAnc_M_Euro$K)) %>%
  melt(zAnc_M_Euro$K) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     name = "Ancestry\ncorrelation") +
  ggtitle("K correlation matrix - w/in European ancestry") +
  xlab("") +
  ylab("")

k_M_plot_all
ggsave("../../bee_manuscript/figures/K_m_plot_all.png",
       plot = k_M_plot_all,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, device = "png")

k_plot_fancy <- k_plot_all +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments, 
               aes(x = starts, xend = ends,
                   y = 0.5, yend = 0.5,
                   color = Region), inherit.aes = F) +
  geom_segment(data = NS_segments, 
               aes(x = 0.5, xend = 0.5,
                   y = starts, yend = ends,
                   color = Region), inherit.aes = F) +  
  # add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5)) +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5)) +
  #geom_text(colour = "darkgray", aes(y = -3, label = zone1),  position = position_dodge(width=0.9))
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A)

k_plot_fancy



# make a new K matrix but omit outlier points:
ZAnc_bees_noOutliers = make_K_calcs(t(A[!(meanA > quantile(meanA, .9) | 
                                            meanA < quantile(meanA, .1)), 
                                        pops_by_lat]))
melt(cov2cor(ZAnc_bees_noOutliers$K)) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     #limits = c(-.35, .55)) +
                     limits = c(-.12, .32)) +
  xlab("") +
  ylab("") +
  ggtitle("K correlation matrix - bees")
ggsave("plots/k_NO_OUTLIERS_correlation_matrix_all_pops.png", 
       height = 7, width = 8, 
       units = "in", device = "png")
ggsave("../../bee_manuscript/figures/k_NO_OUTLIERS_correlation_matrix_all_pops.pdf", 
       height = 7, width = 8, 
       units = "in", device = "pdf")


summary(meanA)
sd(meanA)

# simulate from MVN distribution:
#n_sim <- 10^6 # slow
n_sim <- 10^5
#n_sim = length(meanA) # simulate data of same length
set.seed(101)
MVNsim <- mvrnorm(n = n_sim, 
                  mu = zAnc_bees$alpha, 
                  Sigma = zAnc_bees$K, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)



# sets bounds at 0 and 1
MVNsim_bounded <- data.frame(MVNsim, stringsAsFactors = F) 
MVNsim_bounded[MVNsim < 0] <- 0
MVNsim_bounded[MVNsim > 1] <- 1

MVNsim_AR <- MVNsim[ , AR_pops_included$population]
MVNsim_AR_bounded <- MVNsim_bounded[ , AR_pops_included$population]
MVNsim_CA <- MVNsim[ , CA_pops_included$population]
MVNsim_CA_bounded <- MVNsim_bounded[ , CA_pops_included$population]
# mean across populations
meanA_MVNsim_AR_bounded0 <- apply(MVNsim_AR_bounded, 1, mean)
meanA_MVNsim_CA_bounded0 <- apply(MVNsim_CA_bounded, 1, mean)
# mean across individuals
meanA_MVNsim_AR_bounded <- apply(MVNsim_AR_bounded[ , AR_pops_included$population], 1, 
                                  function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_bounded <- apply(MVNsim_CA_bounded[ , CA_pops_included$population], 1, 
                                  function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
# before bounding:
# mean across individuals
meanA_MVNsim_AR_unbounded <- apply(MVNsim_AR[ , AR_pops_included$population], 1, 
                                 function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_unbounded <- apply(MVNsim_CA[ , CA_pops_included$population], 1, 
                                 function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))


# combined mean across individuals
meanA_MVNsim_bounded <- apply(cbind(MVNsim_AR_bounded, MVNsim_CA_bounded)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))
meanA_MVNsim <- apply(cbind(MVNsim_AR, MVNsim_CA)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))

# save MVN simulation as data objects
save(MVNsim_bounded, meanA_MVNsim_bounded, meanA_MVNsim_AR_bounded, meanA_MVNsim_CA_bounded,
     file = "results/MVNsim_bounded.RData")
#load("results/MVNsim_bounded.RData")

# Add in constraint that all covariances must be 0 or positive
# (effectively I zero out negative covariances):
summary(zAnc_bees$K[zAnc_bees$K < 0])
sum(zAnc_bees$K < 0)/prod(dim(zAnc_bees$K)) # ~13% of covariances < 0
summary(zAnc_bees$K[zAnc_bees$K > 0])
# conclusion: negative covariances are on the same order as positive covariances
summary(zAnc_bees$K[T])
summary(as.matrix(zAnc_bees$K)[lower.tri(as.matrix(zAnc_bees$K), diag = F)]) # don't include diagonals
# simulate from MVN distribution:
K_zero <- zAnc_bees$K
K_zero[zAnc_bees$K < 0] <- 0
set.seed(101)
MVNsim_zero <- mvrnorm(n = n_sim, 
                  mu = zAnc_bees$alpha, 
                  Sigma = K_zero, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)

# how many simulations exceed the bounds?
# overall few, but more than 25% for some populations
sum(MVNsim < 0)/sum(table(MVNsim<0)) # low outliers need to be set to bound
sum(MVNsim > 1)/sum(table(MVNsim>1))
sum(MVNsim_zero < 0)/sum(table(MVNsim_zero<0)) # high outliers need to be set to bound
sum(MVNsim_zero > 1)/sum(table(MVNsim_zero>1))
hist(MVNsim[MVNsim_zero<0])
table(MVNsim_bounded[MVNsim_zero < 1])
apply(MVNsim, 2, function(x) sum(x > 1))/(nrow(MVNsim_zero))
summary(apply(MVNsim, 2, function(x) sum(x <0))/(nrow(MVNsim_zero)))
perc_sim_freq_out_of_bounds = data.frame(population = colnames(MVNsim),
           Lower = apply(MVNsim, 2, function(x) sum(x < 0))/nrow(MVNsim),
           Upper = apply(MVNsim, 2, function(x) sum(x > 1))/nrow(MVNsim),
           stringsAsFactors = F) %>%
  pivot_longer(data = ., cols = c("Lower", "Upper"), names_to = "bound", values_to = "p") %>% 
  left_join(., admix_proportions, by = "population") %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = A, y = p, color = zone)) +
  geom_point(alpha = .75) +
  facet_wrap(~bound) +
  xlab("Mean A ancestry proportion (ancestry_hmm)") +
  ylab("Percent sims exceeding bound") +
  #ggtitle("Percent simulated population ancestry frequencies < 0 (MVN)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic()
perc_sim_freq_out_of_bounds
ggsave("plots/percents_MVN_sims_need_truncation_low.png",
       plot = perc_sim_freq_out_of_bounds,
       height = 3, width = 6, units = "in")
perc_sim_zero_freq_out_of_bounds = data.frame(population = colnames(MVNsim),
                                         Lower = apply(MVNsim_zero, 2, function(x) sum(x < 0))/nrow(MVNsim_zero),
                                         Upper = apply(MVNsim_zero, 2, function(x) sum(x > 1))/nrow(MVNsim_zero),
                                         stringsAsFactors = F) %>%
  pivot_longer(data = ., cols = c("Lower", "Upper"), names_to = "bound", values_to = "p") %>% 
  left_join(., admix_proportions, by = "population") %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = A, y = p, color = zone)) +
  geom_point(alpha = .75) +
  facet_wrap(~bound) +
  xlab("Mean A ancestry proportion (ancestry_hmm)") +
  ylab("Percent sims exceeding bound") +
  #ggtitle("Percent simulated population ancestry frequencies < 0 (MVN)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic()





MVNsim_zero_bounded <- data.frame(MVNsim_zero, stringsAsFactors = F) # sets bounds at 0 and 1
MVNsim_zero_bounded[MVNsim_zero < 0] <- 0
MVNsim_zero_bounded[MVNsim_zero > 1] <- 1
MVNsim_AR_zero <- MVNsim_zero[ , AR_pops_included$population]
MVNsim_AR_zero_bounded <- MVNsim_zero_bounded[ , AR_pops_included$population]
MVNsim_CA_zero <- MVNsim_zero[ , CA_pops_included$population]
MVNsim_CA_zero_bounded <- MVNsim_zero_bounded[ , CA_pops_included$population]
# mean across individuals
meanA_MVNsim_AR_zero_bounded <- apply(MVNsim_AR_zero_bounded[ , AR_pops_included$population], 1, 
                                 function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_zero_bounded <- apply(MVNsim_CA_zero_bounded[ , CA_pops_included$population], 1, 
                                 function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
# combined mean across individuals
meanA_MVNsim_zero_bounded <- apply(cbind(MVNsim_AR_zero_bounded, MVNsim_CA_zero_bounded)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                              function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))
meanA_MVNsim_zero <- apply(cbind(MVNsim_AR_zero, MVNsim_CA_zero)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                      function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))

# save data objects from this simulation:
save(MVNsim_zero_bounded, meanA_MVNsim_zero_bounded, meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_CA_zero_bounded,
     file = "results/MVNsim_zero_bounded.RData")

# other simulations without covariances:
# compare distribution between binomial and MVN for 10^5 simulations
# first read in individual alpha estimates for mean A ancestry
CA_indAlpha <- do.call(rbind,
                       lapply(CA_pops_included$population, function(p) 
                         read.table(paste0("results/ancestry_hmm/combined_sept19/posterior/anc/", p, ".alpha.anc"),
                                    stringsAsFactors = F, header = T)))
AR_indAlpha <- do.call(rbind,
                       lapply(AR_pops_included$population, function(p) 
                         read.table(paste0("results/ancestry_hmm/combined_sept19/posterior/anc/", p, ".alpha.anc"),
                                    stringsAsFactors = F, header = T)))

set.seed(101) # use same seed
PoiBinsim_CA <- apply(do.call(rbind,
                              lapply(CA_indAlpha$A, 
                                     function(alpha)
                                       rbinom(n = n_sim, size = 2, prob = alpha)))/2, 2, mean)
PoiBinsim_AR <- apply(do.call(rbind,
                              lapply(AR_indAlpha$A, 
                                     function(alpha)
                                       rbinom(n = n_sim, size = 2, prob = alpha)))/2, 2, mean)
PoiBinsim_combined <- (PoiBinsim_CA*sum(CA_pops_included$n_bees) + PoiBinsim_AR*sum(AR_pops_included$n_bees))/sum(c(CA_pops_included$n_bees, AR_pops_included$n_bees))

# MVN simulation bounded with no covariances
set.seed(101)
MVNsim_no_cov <- mvrnorm(n = n_sim, 
                         mu = zAnc_bees$alpha, 
                         Sigma = diag(diag(zAnc_bees$K), length(zAnc_bees$alpha), length(zAnc_bees$alpha)), 
                         tol = 1e-6, 
                         empirical = FALSE, 
                         EISPACK = FALSE)
MVNsim_no_cov_bounded <- MVNsim_no_cov
MVNsim_no_cov_bounded[MVNsim_no_cov < 0] <- 0
MVNsim_no_cov_bounded[MVNsim_no_cov > 1] <- 1
MVNsim_AR_no_cov_bounded <- MVNsim_no_cov_bounded[ , AR_pops_included$population]
MVNsim_CA_no_cov_bounded <- MVNsim_no_cov_bounded[ , CA_pops_included$population]
# mean across individuals
meanA_MVNsim_AR_no_cov_bounded <- apply(MVNsim_no_cov_bounded[ , AR_pops_included$population], 1, 
                                        function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_no_cov_bounded <- apply(MVNsim_no_cov_bounded[ , CA_pops_included$population], 1, 
                                        function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
# combined mean across individuals
meanA_MVNsim_no_cov_bounded <- apply(MVNsim_no_cov_bounded[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                                     function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))

#summarise
summary(meanA_MVNsim)
summary(meanA_MVNsim_bounded) # makes little difference
summary(meanA_MVNsim_zero)
summary(meanA_MVNsim_zero_bounded)
# plots
hist(apply(MVNsim, 1, mean))
hist(meanA_MVNsim)
hist(meanA_MVNsim_bounded)
plot(meanA_MVNsim_AR_bounded, meanA_MVNsim_CA_bounded, xlim = c(0, 1), ylim = c(0, 1),
     main = "simulation CA vs. AR frequencies (zero-d K matrix MVN sim in blue)")
points(meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_CA_zero_bounded, col = "blue")
plot(meanA_MVNsim_bounded, meanA_MVNsim_zero_bounded, col = "orange", 
     xlim = c(0.1, 0.7), ylim = c(0.1, 0.7),
     main = "effect of zero-ing out small negative covariances")
abline(a = 0, b = 1)
var(meanA_MVNsim_bounded)
var(meanA_MVNsim_zero_bounded) # very slightly higher variance -- could just be simulation noise
# this makes sense because the covariances I zero-d out were very slight
summary(diag(zAnc_bees$K))
summary(zAnc_bees$K[lower.tri(zAnc_bees$K, diag = F)])

#QQ-plots:
png("plots/QQ_plots_against_MVN.png",
    height = 10, width = 15, res = 300, units = "in")
par(mfrow=c(2,2))
qqplot(meanA_MVNsim_bounded, meanA_MVNsim,
       main = "Effect of 0-1 bounds on MVN distribution")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_bounded, meanA,
       main = "QQ-plot: MVN dist. vs. data, combined mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_CA_bounded, meanA_CA,
       main = "QQ-plot: MVN dist. vs. data, CA mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_AR_bounded, meanA_AR,
       main = "QQ-plot: MVN dist. vs. data, AR mean")
abline(0, 1, col = "blue")
dev.off()

png("plots/QQ_plots_against_MVN_zero-ed_out_negK.png",
    height = 10, width = 15, res = 300, units = "in")
par(mfrow=c(2,2))
qqplot(meanA_MVNsim_zero_bounded, meanA_MVNsim_bounded,
       main = "Effect of non-neg K on MVN distribution")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_zero_bounded, meanA,
       main = "QQ-plot: MVN dist. non-neg K vs. data, combined mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_CA_bounded,
       main = "QQ-plot: Effect of non-neg K on MVN CA mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_AR_bounded,
       main = "QQ-plot: Effect of non-neg K on MVN AR mean")
abline(0, 1, col = "blue")
dev.off()



                  
plot_QQ <- function(NA1, NA2, SA1, SA2, axis1, axis2, ps_qq = seq(0, 1, length.out = 10000)){
  d_qq <- bind_rows(data.frame(p = ps_qq,
                                           zone = "N. America", stringsAsFactors = F) %>%
                                  mutate(q1 = quantile(NA1, p), 
                                         q2 = quantile(NA2, p)),
                                data.frame(p = ps_qq,
                                           zone = "S. America", stringsAsFactors = F) %>%
                                  mutate(q1 = quantile(SA1, p), 
                                         q2 = quantile(SA2, p)))
  ggplot(data = d_qq) + 
    geom_point(aes(x = q1, y = q2, color = p)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab(axis1) +
    ylab(axis2) +
    scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
    theme_classic() +
    facet_grid(.~zone)#, scales = "free_x")
}
# make QQ plots for paper:
# effect of truncating at [0,1] is small:
MVN_vs_01_truncated_MVN_qq = plot_QQ(NA1 = meanA_MVNsim_CA_bounded, NA2 = meanA_MVNsim_CA_unbounded, 
        SA1 = meanA_MVNsim_AR_bounded, SA2 = meanA_MVNsim_AR_unbounded, 
        axis1 = "MVN simulation", axis2 = "MVN simulation unbounded", 
        ps_qq = seq(0, 1, length.out = 10000))
MVN_vs_01_truncated_MVN_qq
ggsave("plots/qq_effect_truncation_on_MVN.png",
       plot = MVN_vs_01_truncated_MVN_qq,
       height = 3, width = 6, 
       units = "in", device = "png")
effect_truncation_on_MVN <- arrangeGrob(perc_sim_freq_out_of_bounds + ggtitle("A"),
                                        MVN_vs_01_truncated_MVN_qq + ggtitle("B"),
                                         nrow = 2,
                                         ncol = 1)
effect_truncation_on_MVN
ggsave("../../bee_manuscript/figures/effect_truncation_on_MVN.png", 
       plot = effect_truncation_on_MVN,
       height = 6, width = 5.2, dpi = 600,
       units = "in", device = "png")
ggsave("../../bee_manuscript/figures_supp/effect_truncation_on_MVN.tiff", 
       plot = effect_truncation_on_MVN,
       height = 6, width = 5.2, dpi = 600,
       units = "in", device = "tiff")
# overall good fit of data to MVN:
mvn_vs_data_qq = plot_QQ(NA1 = meanA_MVNsim_CA_bounded, NA2 = meanA_CA, 
        SA1 = meanA_MVNsim_AR_bounded, SA2 = meanA_AR, 
        axis1 = "MVN simulation", axis2 = "Observed A frequencies", 
        ps_qq = seq(0, 1, length.out = 10000))
mvn_vs_data_qq
# very high expected false-positive rate with poisson-binomial null:
poibin_vs_data_qq = plot_QQ(NA1 = PoiBinsim_CA, NA2 = meanA_CA, 
        SA1 = PoiBinsim_AR, SA2 = meanA_AR, 
        axis1 = "Poisson binomial simulation", axis2 = "Observed A frequencies", 
        ps_qq = seq(0, 1, length.out = 10000))
# and high expected false-positive rate with no-covariance MVN:
mvn_no_cov_vs_data_qq = plot_QQ(NA1 = meanA_MVNsim_CA_no_cov_bounded, NA2 = meanA_CA, 
        SA1 = meanA_MVNsim_AR_no_cov_bounded, SA2 = meanA_AR, 
        axis1 = "MVN simulation with zero covariances", axis2 = "Observed A frequencies", 
        ps_qq = seq(0, 1, length.out = 10000))
# make joint plot for these QQs:
qq_vs_data_plots_combined <- arrangeGrob(mvn_vs_data_qq + ggtitle("A"),
                                         mvn_no_cov_vs_data_qq + ggtitle("B"),
                                         poibin_vs_data_qq + ggtitle("C"),
                                     nrow = 3,
                                     ncol = 1)
# plot them all together as a facet wrap instead:
sim_NA = list(meanA_MVNsim_CA_bounded, meanA_MVNsim_CA_no_cov_bounded, PoiBinsim_CA)
sim_SA = list(meanA_MVNsim_AR_bounded, meanA_MVNsim_AR_no_cov_bounded, PoiBinsim_AR)
sim_names = c("MVN", "MVN zero covariances", "Poisson binomial")
ps_qq = seq(0, 1, length.out = 10000)

d_qq <- do.call(rbind, lapply(1:3, function(i)
  bind_rows(data.frame(p = ps_qq,
             zone = "N. America",
             sim = sim_names[i],
             stringsAsFactors = F) %>%
    mutate(q2 = quantile(meanA_CA, p),
      q1 = quantile(sim_NA[[i]], p)),
    data.frame(p = ps_qq,
             zone = "S. America",
             sim = sim_names[i],
             stringsAsFactors = F) %>%
      mutate(q2 = quantile(meanA_AR, p), 
             q1 = quantile(sim_SA[[i]], p)))))
ggplot(data = d_qq) + 
    geom_point(aes(x = q1, y = q2, color = p), size = 0.5) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Simulated A frequencies") +
    ylab("Observed A frequencies") +
    scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
    theme_classic() +
    facet_grid(zone~sim)#, scales = "free_x")
ggsave("../../bee_manuscript/figures/qq_vs_data_sim_comparison.png",
       width = 7.5, height = 6, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/qq_vs_data_sim_comparison.tiff",
       width = 7.5, height = 6, 
       units = "in", dpi = 600, device = "tiff")
ggsave("plots/qq_vs_data_sim_comparison.png",
       width = 7.5, height = 6, 
       units = "in", dpi = 600, device = "png")



# Compute kde2d
kd_data <- with(data = data.frame(x = meanA_CA, y = meanA_AR), 
           MASS::kde2d(y = y, x = x, n = 100,
                       lims = c(c(0,.6),c(0,.6))))
kd_mvn <- with(data = data.frame(x = meanA_MVNsim_CA_unbounded, 
                                 y = meanA_MVNsim_AR_unbounded), 
               MASS::kde2d(y = y, x = x, n = 100,
                           lims = c(c(0,.6),c(0,.6))))
ggplot(data.frame(kd_mvn)) +
  geom_density(aes(x = x, y = y, z = z))
# Plot with plotly
plot_ly(x = kd_data$x, y = kd_data$y, z = kd_data$z) %>% add_surface()
plot_ly(x = kd_mvn$x, y = kd_mvn$y, z = kd_mvn$z) %>% add_surface()
plot_ly(x = kd_data$x, 
        y = kd_data$y, 
        z = kd_data$z - kd_mvn$z,
        xlab = "N. America",
        ylab = "S. America") %>% 
  add_surface() %>% 
  layout(
    title = "2D Density Observed A frequency minus MVN",
    scene = list(
      xaxis = list(title = "NA", range = c(1, 0)),
      yaxis = list(title = "SA", range = c(1, 0)),
      zaxis = list(title = "Diff.")
    ))


# confirm identical 2d grid points for each density
identical(kd_data$x, kd_mvn$x) == T
identical(kd_data$y, kd_mvn$y) == T

# calc difference between density estimates
kd_diff = kd_data 
kd_diff$z = kd_data$z - kd_mvn$z

image(kd_diff, col = viridis(30))

# melt z data to long format from 100 x 100 matrix
rownames(kd_diff$z) = kd_diff$x
colnames(kd_diff$z) = kd_diff$y
kd_diff$z %>% # melt
  melt(., id.var = rownames(kd_diff)) %>%
  data.table::setnames(c("A_NA", "A_SA", "A_DIFF")) %>%
  ggplot(., aes(x = A_NA, A_SA, z = A_DIFF, fill = A_DIFF)) +
  stat_density2d(aes(alpha = A_DIFF),
                 geom="polygon", bins = 100)
  
  geom_tile() +
  stat_contour(aes(colour = ..level..), binwidth = .01) +
  scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
  scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"), midpoint=0) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  guides(colour=FALSE)

plot_2d_density <- function(){
  contour(kd_mvn, col = viridis(10)[1], nlevels = 7, lwd = 2,
          xlim = c(0.12, 0.35),
          ylim = c(0.32, 0.5),
          xlab = "N. America",
          ylab = "S. America")
  contour(kd_data, add = T, col = viridis(10)[9], cex = 2, nlevels = 7)
  legend("bottomright", legend = c("Observed A frequency", "Simulated A frequency (MVN)"),
         lwd = 1, lty = 1, col = viridis(10)[c(9,1)], cex = 0.5)
}
  
png(file = "../../bee_manuscript/figures/comparison_2d_density_data_mvn.png",
    height = 5, width = 5.2, res = 600, units = "in")
plot_2d_density()
dev.off()
tiff(file = "../../bee_manuscript/figures_supp/comparison_2d_density_data_mvn.tiff",
    height = 5, width = 5.2, res = 600, units = "in")
plot_2d_density()
dev.off()
png(file = "plots/comparison_2d_density_data_mvn.png",
    height = 5, width = 5.2, res = 600, units = "in")
plot_2d_density()
dev.off()

# for how long do my quantiles match up?
qs <- seq(0, 1, by = .01)
qs_meanA_MVNsim_bounded <- quantile(meanA_MVNsim_bounded, probs = qs) 
qs_meanA <- quantile(meanA, probs = qs) 
qs_meanA_AR <- quantile(meanA_AR, probs = qs) 
qs_meanA_CA <- quantile(meanA_CA, probs = qs) 
qs_meanA_MVNsim_CA_bounded <- quantile(meanA_MVNsim_CA_bounded, probs = qs)
qs_meanA_MVNsim_AR_bounded <- quantile(meanA_MVNsim_AR_bounded, probs = qs)
qqplot(meanA_MVNsim_bounded, meanA)
points(qs_meanA_MVNsim_bounded, qs_meanA, col = "blue")
tail(qs_meanA_MVNsim_bounded)
tail(qs_meanA)
head(qs_meanA_MVNsim_bounded)
head(qs_meanA)
lapply(list(qs_meanA_AR, qs_meanA_MVNsim_AR_bounded, 
            qs_meanA_CA, qs_meanA_MVNsim_CA_bounded, 
            qs_meanA, qs_meanA_MVNsim_bounded), tail)
lapply(list(qs_meanA_AR, qs_meanA_MVNsim_AR_bounded, 
            qs_meanA_CA, qs_meanA_MVNsim_CA_bounded, 
            qs_meanA, qs_meanA_MVNsim_bounded), head)
# probably just removing top 1% is fine, but I could be safe and remove top 3% outliers
high = "97%"
low = "3%"
is_outlier <- (meanA >= qs_meanA[high] | meanA_AR >=  qs_meanA_AR[high] | meanA_CA >= qs_meanA_CA[high] |
                 meanA <= qs_meanA[low] | meanA_AR <=  qs_meanA_AR[low] | meanA_CA <= qs_meanA_CA[low])
# what % are outliers high or low? About 12% of my data. I'll exclude this set for demographic inference
table(is_outlier)/length(is_outlier)

# plot 50%, 95%, 99% & .999% credible intervals for my empirical observations
LaplacesDemon::joint.pr.plot(meanA_CA, 
                             meanA_AR, # top 1% may not be strict enough to ID outliers
                             quantiles=c(0.5, 0.95, 0.99, .999))


png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_outliers_blue.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(apply(MVNsim_CA_bounded, 1, mean), .9999) | 
                                                   meanA_AR > quantile(apply(MVNsim_AR_bounded, 1, mean), .9999)) |
                                                  (meanA_CA < quantile(apply(MVNsim_CA_bounded, 1, mean), .0001) | 
                                                     meanA_AR < quantile(apply(MVNsim_AR_bounded, 1, mean), .0001)), 
                                                alpha("blue", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, blue = top 0.01% sim")
dev.off()


png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_results_orange.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = alpha("grey", .1),
     main = "comparison of A ancestry in California vs. Argentina zone, orange = 100k MVN simulations")
points(apply(MVNsim_CA_bounded, 1, mean), apply(MVNsim_AR_bounded, 1, mean), pch = 20, col = alpha("orange", .02))
dev.off()
#points(apply(MVNsim_CA, 1, mean), apply(MVNsim_AR, 1, mean), pch = 20, col = alpha("green", .1))

# plot with both outliers in blue and MVN sims in orange
png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_results_orange_outliers_blue.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(apply(MVNsim_CA_bounded, 1, mean), .9999) | 
                                                   meanA_AR > quantile(apply(MVNsim_AR_bounded, 1, mean), .9999)) |
                                                  (meanA_CA < quantile(apply(MVNsim_CA_bounded, 1, mean), .0001) | 
                                                     meanA_AR < quantile(apply(MVNsim_AR_bounded, 1, mean), .0001)), 
                                                alpha("blue", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, blue = top 0.01% sim")
points(apply(MVNsim_CA_bounded, 1, mean), apply(MVNsim_AR_bounded, 1, mean), pch = 20, col = alpha("orange", .02))
dev.off()

# how likely are these shared outliers in my simulation?
MVNsim_AR_bounded_quantile <- quantile(meanA_MVNsim_AR_bounded, .99)
MVNsim_CA_bounded_quantile <- quantile(meanA_MVNsim_CA_bounded, .99)
table(meanA_MVNsim_AR_bounded > MVNsim_AR_bounded_quantile & 
        meanA_MVNsim_CA_bounded > MVNsim_CA_bounded_quantile)/n_sim
.01^2 # it's about 5-6 times as likely under the MVN to get a double outlier at top .01% than it would be under true independence between N and S America
# but still same order of magnitude. Also I maybe need to simulate a larger # of trials to get accuracy in the tails (only expect 1000 outliers each out of 100k, very small # overlap)
# when I simulated 1 million MVN's I got .000521 as the probability of outliers in both CA and AR at .01%
MVNsim_AR_bounded_quantile_low <- quantile(meanA_MVNsim_AR_bounded, .01)
MVNsim_CA_bounded_quantile_low <- quantile(meanA_MVNsim_CA_bounded, .01)
table(meanA_MVNsim_AR_bounded < MVNsim_AR_bounded_quantile_low & 
        meanA_MVNsim_CA_bounded < MVNsim_CA_bounded_quantile_low)/n_sim
# on the low side it's similar. some variance from low # expected hits in my simulations


# these outliers are pretty surprising under a MVN framework
# how about under a poisson binomial framework? even more surprising

# compare simulations
# plot comparison
sim_compare <- 
  data.frame(poisson_binomial = PoiBinsim_combined,
           MVN_no_covariance = meanA_MVNsim_no_cov_bounded,
           #observed_data = sample(meanA, n_sim, replace = F), # downsample data to match length of simulations
           MVN_with_covariance = meanA_MVNsim_bounded) %>%
  tidyr::gather(., "distribution", "A") %>%
  bind_rows(., data.frame(distribution = "observed_data", A = meanA, stringsAsFactors = F))
           
# plot
p_sim_compare <- sim_compare %>%
  ggplot(., aes(x = A, color = distribution)) +
  #geom_density_line(lwd = 1) +
  geom_line(stat = "density", #alpha = .75, 
            aes(linetype = distribution)) +
  theme_classic() +
  xlab("Combined sample mean African ancestry") +
  ylab("Density") +
  scale_color_viridis_d(option = "viridis", name = NULL, 
                        limits = c("observed_data", "MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_linetype_manual(name = NULL, values = c(1,2,3,4),
                        limits = c("observed_data", "MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
#scale_color_discrete(name = "Distribution", labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
p_sim_compare
p_sim_compare2 <- sim_compare %>%
  ggplot(., aes(x = A, color = distribution)) +
  #geom_density_line(lwd = 1) +
  geom_line(stat = "density", alpha = .75, lwd = 1) +
  theme_classic() +
  xlab("Mean African ancestry") +
  ylab("Density") +
  scale_color_viridis_d(option = "viridis", name = NULL, 
                        limits = c("observed_data", "MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
# use alpha instead of different line types
p_sim_compare2

ggsave("plots/distribution_data_vs_poibin_vs_MVN_sim.png",
       plot = p_sim_compare2,
       device = "png",
       width = 5.2, height = 3, units = "in")
ggsave("../../bee_manuscript/figures/distribution_data_vs_poibin_vs_MVN_sim.png",
       plot = p_sim_compare2,
       device = "png",
       width = 5.2, height = 3, 
       units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/distribution_data_vs_poibin_vs_MVN_sim.tiff",
       plot = p_sim_compare2,
       device = "tiff",
       width = 5.2, height = 3, 
       units = "in", dpi = 600)

# make joint plot of kinship matrix and these distribution comparisons:

dist_k_plots_combined <- arrangeGrob(k_plot_fancy + ggtitle("A"),
                                      p_sim_compare2 + ggtitle("B"),
                                      ncol = 2,
                                      widths = c(3, 3))

ggsave("../../bee_manuscript/figures/k_matrix_and_poi_bin_mvn_dist_comparison.png",
       plot = dist_k_plots_combined,
       device = "png",
       width = 7.5, 
       height = 3, 
       units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/k_matrix_and_poi_bin_mvn_dist_comparison.tiff",
       plot = dist_k_plots_combined,
       device = "tiff",
       width = 7.5, 
       height = 3, 
       units = "in", dpi = 600)
ggsave("plots/k_matrix_and_poi_bin_mvn_dist_comparison.png",
       plot = dist_k_plots_combined,
       device = "png",
       width = 7.5, 
       height = 3, 
       units = "in", dpi = 600)


png(paste0("plots/QQ_plot_data_against_PoiBin.png"),
    height = 8, width = 8, res = 300, units = "in")
qqplot(PoiBinsim_combined, meanA,
       main = "QQ plot fit - A Ancestry Combined Sample",
       xlab = "Poisson Binomial Simulation",
       ylab = "Observed",
       col = "grey")
abline(0, 1, col = "black")
dev.off()
#make_qqplot_lines(meanA_MVNsim_zero_bounded, legend = T)

# plot QQ-plots:
#QQ-plots:
png("plots/QQ_plots_against_PoiBin.png",
    height = 10, width = 15, res = 300, units = "in")
par(mfrow=c(2,2))
qqplot(PoiBinsim_combined, meanA_MVNsim_bounded,
       main = "Poisson Binomial vs. MVN distribution null")
abline(0, 1, col = "blue")
qqplot(PoiBinsim_combined, meanA,
       main = "QQ-plot: Poisson Binomial dist. vs. data, combined mean")
abline(0, 1, col = "blue")
qqplot(PoiBinsim_CA, meanA_CA,
       main = "QQ-plot: Poisson Binomial dist. vs. data, CA mean")
abline(0, 1, col = "blue")
qqplot(PoiBinsim_AR, meanA_AR,
       main = "QQ-plot: Poisson Binomial dist. vs. data, AR mean")
abline(0, 1, col = "blue")
dev.off()


# What are my false discovery rates?
# get a 1% and 5% FDR for jointly shared outliers under the MVN (based on simulations):

# walk along sd's above mean ancestry for each zone
sd_CA <- sd(meanA_CA)
mu_CA <- mean(meanA_CA)
sd_AR <- sd(meanA_AR)
mu_AR <- mean(meanA_AR)
sd_range <- seq(0, 5, by = .005) # I can make this more precise later if I want

# get standard deviation cutoffs for FDR shared high loci
test_fdr_shared_high_sds <- sapply(sd_range, function(x)
  fdr_shared_high(a1 = mu_AR+x*sd_AR, a2 = mu_CA+x*sd_CA, pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
FDRs_shared_high_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_shared_high_sds<p], na.rm = T))

# shared low # standard deviations below mean
test_fdr_shared_low_sds <- sapply(sd_range, function(x)
  fdr_shared_low(a1 = mu_AR-x*sd_AR, a2 = mu_CA-x*sd_CA, pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
FDRs_shared_low_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_shared_low_sds<p], na.rm = T))



range_q_high <- seq(0.75, 1, by = 0.0001)
range_q_low <- seq(0, 0.25, by = 0.0001)

# high
test_fdr_CA_high2 <- sapply(range, function(x) 
  fdr_1pop_high(a = x, pop = meanA_CA, 
                sims = meanA_MVNsim_CA_bounded))
FDRs_CA_high2 <- sapply(FDR_values, function(p) min(range_q_high[test_fdr_CA_high2<p], na.rm = T))
# get standard deviation cutoffs for FDR shared high loci -- using MVN quantiles
test_fdr_shared_high_quantile <- sapply(range_q_high, function(x)
  fdr_shared_high(a1 = quantile(meanA_MVNsim_AR_bounded, x), a2 = quantile(meanA_MVNsim_CA_bounded, x), 
                  pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
FDRs_shared_high_quantile <- sapply(FDR_values, function(q) min(range_q_high[test_fdr_shared_high_quantile<q], 
                                    na.rm = T))
# about how many outliers?
table(meanA_CA > quantile(meanA_MVNsim_CA_bounded, FDRs_shared_high_quantile[2]) & 
        meanA_AR > quantile(meanA_MVNsim_AR_bounded, FDRs_shared_high_quantile[2]))/length(meanA_CA)
A_AR_CA %>%
  filter(., CA > quantile(meanA_MVNsim_CA_bounded, FDRs_shared_high_quantile[2]) & 
  AR > quantile(meanA_MVNsim_AR_bounded, FDRs_shared_high_quantile[2])) %>%
  group_by(chr) %>%
  summarise(n = n())


# shared low # standard deviations below mean
test_fdr_shared_low_quantile <- sapply(range_q_low, function(x)
  fdr_shared_low(a1 = quantile(meanA_MVNsim_AR_bounded, x), a2 = quantile(meanA_MVNsim_CA_bounded, x), 
                 pop1 = meanA_AR, pop2 = meanA_CA, 
                 sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
FDRs_shared_low_quantile <- sapply(FDR_values, function(q) min(range_q_low[test_fdr_shared_low_quantile<q], na.rm = T))
# no outliers


# high CA and high AR separately (less powerful):
test_fdr_CA_high_sds <- sapply(sd_range, function(x) 
  fdr_1pop_high(a = mu_CA+x*sd_CA, pop = meanA_CA, 
                  sims = meanA_MVNsim_CA_bounded))
FDRs_CA_high_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_CA_high_sds<p], na.rm = T))

test_fdr_AR_high_sds <- sapply(sd_range, function(x)
  fdr_1pop_high(a = mu_AR+x*sd_AR, pop = meanA_AR, 
                sims = meanA_MVNsim_AR_bounded))
FDRs_AR_high_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_AR_high_sds<p], na.rm = T))

# low CA and high AR separately (less powerful):
test_fdr_CA_low_sds <- sapply(sd_range, function(x) 
  fdr_1pop_low(a = mu_CA-x*sd_CA, pop = meanA_CA, 
                sims = meanA_MVNsim_CA_bounded))
FDRs_CA_low_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_CA_low_sds<p], na.rm = T))
# we do not find any low outliers in CA based on a .1, .05 or .01 FDR (underpowered)
test_fdr_AR_low_sds <- sapply(sd_range, function(x)
  fdr_1pop_low(a = mu_AR-x*sd_AR, pop = meanA_AR, 
               sims = meanA_MVNsim_AR_bounded))
FDRs_AR_low_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_AR_low_sds<p], na.rm = T))
# get standard deviations for false-discovery rates
FDRs_sds = data.frame(FDR_values = FDR_values,
                      shared_high = FDRs_shared_high_sds,
                      shared_low = FDRs_shared_low_sds,
                      CA_high = FDRs_CA_high_sds,
                      CA_low = FDRs_CA_low_sds,
                      AR_high = FDRs_AR_high_sds,
                      AR_low = FDRs_AR_low_sds,
                      stringsAsFactors = F)
write.table(FDRs_sds, "results/FDRs_MVN_01_high_low_A_SDs.txt", quote = F, col.names = T, row.names = F, sep = "\t")
data.frame(mu = c(mu_CA, mu_AR), sd = c(sd_CA, sd_AR), 
           zone = c("N. America", "S. America"), short_name = c("CA", "AR"),
           stringsAsFactors = F) %>% 
  write.table(., "results/mu_sd_CA_AR.txt", quote = F, col.names = T, row.names = F, sep = "\t")


# get FDR's more directly, without using mean and sd:
range01 <- seq(0, 1, by = .0001)

# high
test_fdr_CA_high2 <- sapply(range01, function(x) 
  fdr_1pop_high(a = x, pop = meanA_CA, 
                sims = meanA_MVNsim_CA_bounded))
FDRs_CA_high2 <- sapply(FDR_values, function(p) min(range01[test_fdr_CA_high2<p], na.rm = T))
test_fdr_AR_high2 <- sapply(range01, function(x) 
  fdr_1pop_high(a = x, pop = meanA_AR, 
                sims = meanA_MVNsim_AR_bounded))
FDRs_AR_high2 <- sapply(FDR_values, function(p) min(range01[test_fdr_AR_high2<p], na.rm = T))

# low
test_fdr_CA_low2 <- sapply(range01, function(x) 
  fdr_1pop_low(a = x, pop = meanA_CA, 
               sims = meanA_MVNsim_CA_bounded))
FDRs_CA_low2 <- sapply(FDR_values, function(p) max(range01[test_fdr_CA_low2<p], na.rm = T))


test_fdr_AR_low2 <- sapply(range01, function(x) 
  fdr_1pop_low(a = x, pop = meanA_AR, 
               sims = meanA_MVNsim_AR_bounded))
FDRs_AR_low2 <- sapply(FDR_values, function(p) max(range01[test_fdr_AR_low2<p], na.rm = T))


# translate SD cutoffs to % A ancestry cutoffs:
FDRs = data.frame(FDR_values = FDR_values,
                  #shared_high_CA = mu_CA + FDRs_sds$shared_high*sd_CA,
                  #shared_high_AR = mu_AR + FDRs_sds$shared_high*sd_AR,
                  #shared_low_CA = mu_CA - FDRs_sds$shared_low*sd_CA,
                  #shared_low_AR = mu_AR - FDRs_sds$shared_low*sd_AR,
                  CA_high = FDRs_CA_high2, #mu_CA + FDRs_sds$CA_high*sd_CA,
                  CA_low = FDRs_CA_low2, #mu_CA - FDRs_sds$CA_low*sd_CA,
                  AR_high = FDRs_AR_high2, #mu_AR + FDRs_sds$AR_high*sd_AR,
                  AR_low = FDRs_AR_low2, #mu_AR - FDRs_sds$AR_low*sd_AR,
                  stringsAsFactors = F)
write.table(FDRs, "results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", quote = F, col.names = T, row.names = F, sep = "\t")
FDRs <- read.table("results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", 
                   stringsAsFactors = F, header = T, sep = "\t")


# what percent of the genome matches these cutoffs?
table(meanA_CA > FDRs$CA_high[FDRs$FDR_values==.05])/length(meanA_CA) # 0.26% of loci
table(meanA_CA > FDRs$CA_high[FDRs$FDR_values==.1])/length(meanA_CA) # 0.34% of loci
table(sites$scaffold[meanA_CA > FDRs$CA_high[FDRs$FDR_values==.05]]) # where are they?
table(meanA_MVNsim_CA_bounded > FDRs$CA_high[FDRs$FDR_values==.1])/n_sim # quantile from sims

# high AR
table(meanA_AR > FDRs$AR_high[FDRs$FDR_values==.05])/length(meanA_AR) # 0.06% of loci (very few!)
table(meanA_AR > FDRs$AR_high[FDRs$FDR_values==.1])/length(meanA_AR) # 0.13% of loci (very few!)
table(sites$scaffold[meanA_AR > FDRs$AR_high[FDRs$FDR_values==.05]]) # where are they?
table(meanA_MVNsim_AR_bounded > FDRs$AR_high[FDRs$FDR_values==.05])/n_sim # quantile from sims
quantile(meanA_MVNsim_AR_bounded, .999)

# about .3% of the genome is found in high shared sites at 5% FDR
table(meanA_CA > FDRs$CA_high[FDRs$FDR_values == .05] & meanA_AR > FDRs$AR_high[FDRs$FDR_values == .05])/length(meanA_CA)
table(sites$scaffold[meanA_CA > FDRs$CA_high[FDRs$FDR_values == .05] & meanA_AR > FDRs$AR_high[FDRs$FDR_values == .05]])
table(sites$scaffold[meanA_CA > FDRs$CA_high[FDRs$FDR_values == .1] & meanA_AR > FDRs$AR_high[FDRs$FDR_values == .1]])
# I get a similar answer if I use the simulation quantiles to find outliers (few regions on a few chr)
table(sites$scaffold[meanA_CA > quantile(meanA_MVNsim_CA_bounded, 0.99) & meanA_AR > quantile(meanA_MVNsim_AR_bounded, 0.99)])
table(sites$scaffold[meanA_CA > quantile(meanA_MVNsim_CA_bounded, 0.999) & meanA_AR > quantile(meanA_MVNsim_AR_bounded, 0.999)])



# make a grid where you calculate quantile for joint distribution based off of MVN simulation
grid_n <- 10 # smoother if you use more like 100
range_CA <- seq(from = range(meanA_CA)[1], to = range(meanA_CA)[2], length.out = grid_n)
range_AR <- seq(from = range(meanA_AR)[1], to = range(meanA_AR)[2], length.out = grid_n)
quantile_MVN <- matrix(0, grid_n, grid_n)
for (i in 1:grid_n){
  for (j in 1:grid_n){
    quantile_MVN[i, j] <- mean(abs(meanA_MVNsim_CA_zero_bounded - mu_CA) <= abs(range_CA[i] - mu_CA) & 
                                 abs(meanA_MVNsim_AR_zero_bounded - mu_AR) <= abs(range_AR[j] - mu_AR)) 
  }
}
png("plots/contour_plot_joint_outliers_CA_AR_zero_bounded_MVN.png",
    height = 6, width = 8, units = "in", res = 300)
contour(x = range_CA, y = range_AR, z = quantile_MVN, 
        levels = c(0.25, 0.5, .95, .99, 1),
        xlab = "mean A ancestry in CA hybrid zone",
        ylab = "mean A ancestry in AR hybrid zone",       
        main = "contour plot for joint outliers, MVN zero bounded simulation")
points(meanA_CA, meanA_AR, col = alpha("blue", .05), pch = 20)
dev.off()

# actually maybe what I want is simpler -- just a density plot for the MVN sim underneath my data points
ggplot() +
  geom_density_2d(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
             aes(x = meanA_CA, y = meanA_AR), color = "blue") +
  geom_density_2d(data = data.frame(MVN_CA = meanA_MVNsim_CA_zero_bounded,
                                    MVN_AR = meanA_MVNsim_AR_zero_bounded), 
                  aes(x=MVN_CA, y=MVN_AR), color = "orange")
ggplot() +
  stat_density_2d(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
                  aes(x = meanA_CA, y = meanA_AR, fill = stat(level)), geom = "hex")

cor(meanA_CA, meanA_AR)
cor(meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_CA_zero_bounded)
# I can use stat_contour if I have a z gridded for what I want to display


# make a new kind of density plot:
p_density <- ggplot(data = data.frame(CA = meanA_CA, AR = meanA_AR), #[c(T, rep(F, 10)), ],
       aes(x = CA, y = AR)) +
  #geom_point()
  #geom_bin2d(bins = 100) +
  #stat_density2d() +
  ggpointdensity::geom_pointdensity(adjust = 2, size = .5) + 
  # change method away from kde2d for original k neighbor plot, but may have to thin
  # takes way longer -- also need to rasterize
  scale_color_viridis(option = "magma") +
  scale_fill_viridis(option = "magma") +
  theme_classic() +
  xlab("N. America") +
  ylab("S. America") #+
  #labs(color = "No. Nearby Points")
plot(p_density +
       geom_density_2d(data = data.frame(MVN_CA = meanA_MVNsim_CA_bounded,
                                         MVN_AR = meanA_MVNsim_AR_bounded), 
                       aes(x=MVN_CA, y=MVN_AR), color = "orange"))

b <- as.mcmc(cbind(meanA_MVNsim_CA_bounded, meanA_MVNsim_AR_bounded))
a <- emdbook::HPDregionplot(b, vars = c("meanA_MVNsim_CA_bounded", "meanA_MVNsim_AR_bounded"), 
                            prob = c(0.99, .9, .75), 
                            n = 100, # number of grid points
                            #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
                            #col=c("lightgreen", "salmon", "lightblue", "yellow"), 
                            lwd = 3, 
                            h = .05,
                            add = FALSE)
b_data <- as.mcmc(cbind(meanA_CA, meanA_AR))
a_data <- emdbook::HPDregionplot(b_data, vars = c("meanA_CA", "meanA_AR"), 
                                 prob = c(0.99, .9, .75), 
                                 n = 100, # number of grid points
                                 #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
                                 #col=c("lightgreen", "salmon", "lightblue", "yellow"), 
                                 lwd = 3, 
                                 h = .05,
                                 add = TRUE)
b_MVNnoCov <- as.mcmc(cbind(meanA_MVNsim_CA_no_cov_bounded, meanA_MVNsim_AR_no_cov_bounded))
a_MVNnoCov <- emdbook::HPDregionplot(b_MVNnoCov, vars = c("meanA_MVNsim_CA_no_cov_bounded", "meanA_MVNsim_AR_no_cov_bounded"), 
                            prob = c(0.99, .9, .75), 
                            n = 100,
                            lwd = 3, 
                            h = .05,
                            add = FALSE)
b_poiBin <- as.mcmc(cbind(PoiBinsim_CA, PoiBinsim_AR))
a_poiBin <- emdbook::HPDregionplot(b_poiBin, vars = c("PoiBinsim_CA", "PoiBinsim_AR"), 
                                     prob = c(0.99, .9, .75), 
                                     n = 100,
                                     lwd = 3, 
                                     h = .05,
                                     add = FALSE)


p_density2 <- ggplot() +
  geom_hex(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
           aes(x = meanA_CA, y = meanA_AR, color = ..count.., fill = ..count..), bins = 200) +
  scale_fill_viridis(option = "magma", end = 1, begin = 0.1, name = "SNP count") +
  scale_color_viridis(option = "magma", end = 1, begin = 0.1, name = "SNP count") +
  #theme_classic() +
  theme_light() +
  coord_fixed() +
  theme(panel.grid.minor = element_blank()) +
  xlab("African ancestry in North America") +
  ylab("African ancestry in South America") +
  geom_polygon(data = data.frame(a[[1]]), 
               aes(x = x, y = y),
               fill = NA, color = "orange", lwd = 1)
plot(p_density2)

ggsave("plots/A_CA_vs_AR_density.png",
       plot = p_density2,
       height = 4.3, width = 5.2, 
       units = "in", dpi = 600,
       device = "png")
ggsave("../../bee_manuscript/figures/A_CA_vs_AR_density.png",
       plot = p_density2,
       height = 4.3, width = 5.2, 
       units = "in", dpi = 600, 
       device = "png")
ggsave("../../bee_manuscript/figures_main/A_CA_vs_AR_density.tiff",
       plot = p_density2,
       height = 4.3, width = 5.2, 
       units = "in", dpi = 600, 
       device = "tiff")


sim_compare_CA <- 
  data.frame(poisson_binomial = PoiBinsim_CA,
             MVN_no_covariance = meanA_MVNsim_CA_no_cov_bounded,
             MVN_with_covariance = meanA_MVNsim_CA_bounded) %>%
  tidyr::gather(., "distribution", "A") %>%
  bind_rows(., data.frame(distribution = "observed_data", A = meanA_CA, stringsAsFactors = F)) %>%
  mutate(distribution = factor(distribution, levels = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial", "observed_data"), ordered = T))

p_sim_compare_CA_with_legend <- sim_compare_CA %>%
  filter(., distribution != "observed_data") %>%
  ggplot(., aes(x = A, color = distribution)) +
  geom_histogram(data = filter(sim_compare_CA, distribution == "observed_data"),
                 aes(y = ..density..), bins = 40, fill = "white", color = "black", lwd = 0.3) +
  geom_line(stat = "density", alpha = 0.75, lwd = 1, aes(linetype = distribution)) +
  theme_classic() +
  xlab("Mean African ancestry") +
  ylab("Density") +
  #viridis(4)[c(1,4,3,2)]
  scale_color_manual(values = c("orange", viridis(4)[2:3]), name = "Model", 
                        limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_linetype_manual(values=c(1, 2, 3), name = "Model",
                        limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  #xlim(c(0.05, 0.6)) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     labels = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     limits = c(0.05, 0.6)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "in"),
        panel.border = element_blank(), 
        panel.spacing = unit(0, "in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p_sim_compare_CA_with_legend
p_sim_compare_CA <- p_sim_compare_CA_with_legend + guides(color = F, linetype = F)


sim_compare_AR <- 
  data.frame(poisson_binomial = PoiBinsim_AR,
             MVN_no_covariance = meanA_MVNsim_AR_no_cov_bounded,
             MVN_with_covariance = meanA_MVNsim_AR_bounded) %>%
  tidyr::gather(., "distribution", "A") %>%
  bind_rows(., data.frame(distribution = "observed_data", A = meanA_AR, stringsAsFactors = F)) %>%
  mutate(distribution = factor(distribution, levels = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial", "observed_data"), ordered = T))

p_sim_compare_AR <- sim_compare_AR %>%
  filter(., distribution != "observed_data") %>%
  ggplot(., aes(x = A, color = distribution)) +
  geom_histogram(data = filter(sim_compare_AR, distribution == "observed_data"),
                 aes(y = ..density..), bins = 40, fill = "white", color = "black", lwd = 0.3) +
  geom_line(stat = "density", alpha = 0.75, lwd = 1, aes(linetype = distribution)) +
  theme_classic() +
  xlab("Mean African ancestry") +
  ylab("Density") +
  scale_color_manual(values = c("orange", viridis(4)[2:3]), name = NULL, 
                     limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                     labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_linetype_manual(values=c(1, 2, 3), name = "Model",
                        limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  guides(color = F, linetype = F) +
  xlim(c(0.15, 0.7)) +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        axis.line.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "in"),
        panel.border = element_blank(), 
        panel.spacing = unit(0, "in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p_sim_compare_AR

# extract legends
legend_p_sim_compare_bottom <- get_legend(p_sim_compare_CA_with_legend +
                                     theme(plot.margin = unit(c(1,1,1,1), "lines"),
                                           legend.box.margin = unit(c(1,1,1,1), "lines"),
                                           legend.position = "bottom"))
legend_p_sim_compare <- get_legend(p_sim_compare_CA_with_legend +
                                     theme(plot.margin = unit(c(1,1,1,1), "lines"),
                                           legend.box.margin = unit(c(1,1,1,1), "lines")))
legend_p_density3 <- get_legend(p_density2 +
                                  theme(legend.position = "right",
                                    plot.margin = unit(c(1,1,1,1), "lines"),
                                        legend.box.margin = unit(c(1,1,1,1), "lines")))

p_density3 <- p_density2 + 
  #guides(color = F, fill = F) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0, unit='cm'),
        legend.title = element_text(size = 10)#,
        #axis.text = element_blank(),
        #axis.ticks = element_blank()
        ) +
  guides(fill = guide_colorbar(title.vjust = 0.75)) +
  #theme(legend.position = c(0.9, 0.25)) +
  #xlim(c(0.05, 0.6)) +
  ylim(c(0.15, 0.7)) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     labels = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     limits = c(0.05, 0.6))#+
  #geom_vline(xintercept = c(0.2, 0.4, 0.6)) +
  #geom_hline(yintercept = c(0.2, 0.4, 0.6))
plot(p_density3)


g0 <- ggplotGrob(p_density3)
#plot(g)
panel_id <- g0$layout[g0$layout$name == "panel", c("t","l")]
#g <- grid.arrange(nrow = 2, grobs = list(g0, legend_p_sim_compare_bottom$grobs[[1]]))
g <- g0
g <- gtable_add_cols(g, unit(1.5, "in"))
g <- gtable_add_grob(g, ggplotGrob(p_sim_compare_AR + 
                                     theme(axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
                                           axis.line = element_line(colour = "grey87", size = rel(0.5)))),
                     t = panel_id$t, l = ncol(g))
g <- gtable_add_rows(g, unit(1.5, "in"), 0) # add row on top
g <- gtable_add_grob(g, ggplotGrob(p_sim_compare_CA + 
                                     theme(axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
                                           axis.line = element_line(colour = "grey87", size = rel(0.5)))),
                     t = 1, l = panel_id$l)
g <- gtable_add_rows(g, unit(0.4, "in"), -1)
g <- gtable_add_grob(g, legend_p_sim_compare_bottom$grobs[[1]],
                     t = -1, l = panel_id$l) # add row to bottom?
grid.newpage()
grid.draw(g)
ggsave("plots/A_CA_vs_AR_density_marginal_hist.png",
       plot = g,
       height = 5.4, width = 5.2, 
       units = "in", dpi = 600,
       device = "png")
ggsave("../../bee_manuscript/figures/A_CA_vs_AR_density_marginal_hist.png",
       plot = g,
       height = 5.4, width = 5.2, 
       #height = 6.4, width = 6, 
       units = "in", dpi = 600, 
       device = "png")
ggsave("../../bee_manuscript/figures_main/A_CA_vs_AR_density_marginal_hist.tif",
       plot = g,
       height = 5.4, width = 5.2, 
       units = "in", dpi = 600, 
       device = "tiff",
       compression = "lzw", type = "cairo"
       )

# what % of data points fall within approximate 99% HPDI for 
# each model?

# plot
p_density2 +
  geom_polygon(data = data.frame(a_poiBin[[1]]), 
               aes(x = x, y = y),
               fill = NA, color = "blue", lwd = 1) +
  geom_polygon(data = data.frame(a_MVNnoCov[[1]]), 
               aes(x = x, y = y),
               fill = NA, color = "green", lwd = 1) +
  geom_polygon(data = data.frame(a_data[[1]]), 
               aes(x = x, y = y),
               fill = NA, color = "red", lwd = 1) 

poly_MVN <- Polygons(srl = list(Polygon(coords = with(a[[1]], cbind(x, y)), hole = F)),
                     ID = "MVN")
poly_MVNNoCov <- Polygons(srl = list(Polygon(coords = with(a_MVNnoCov[[1]], cbind(x, y)), hole = F)),
                          ID = "NoCov")
poly_PoiBin <- Polygons(srl = list(Polygon(coords = with(a_poiBin[[1]], cbind(x, y)), hole = F)),
                        ID = "Poibin")
polygons <- SpatialPolygons(Srl = list(poly_MVN, poly_MVNNoCov, poly_PoiBin))
plot(polygons)
some_points <- A_AR_CA %>%
  dplyr::select(CA, AR) %>%
  head() %>%
  SpatialPoints(coords = .)
# % of points. note: it's approximate in part because polygons are approximately defined with finite # of points
# and also because I use the sample of 100k points from these distributions to approximate the density rather than the true density.
GISTools::poly.counts(pts = some_points, polys = polygons)/nrow(A_AR_CA)


# # thin points and plot nearest neighbors:
# p_neighbors <- ggplot(data = data.frame(CA = meanA_CA, AR = meanA_AR)[c(T, rep(F, 100)), ],
#                     aes(x = CA, y = AR)) +
#   #geom_point()
#   #geom_bin2d(bins = 100) +
#   #stat_density2d() +
#   ggpointdensity::geom_pointdensity(adjust = .1, size = .5) + 
#   # change method away from kde2d for original k neighbor plot, but may have to thin
#   # takes way longer -- also need to rasterize
#   scale_color_viridis(option = "magma") +
#   scale_fill_viridis(option = "magma") +
#   theme_classic() +
#   xlab("N. America") +
#   ylab("S. America") #+
# #labs(color = "No. Nearby Points")
# plot(p_neighbors)
# ggsave("plots/A_CA_vs_AR_neighbors.png",
#        plot = p_density,
#        height = 5, width = 6, units = "in", device = "png")
# ggsave("../../bee_manuscript/figures/A_CA_vs_AR_neighbors.png",
#        height = 5, width = 6, units = "in", device = "png")



# plot(meanA_CA, meanA_AR, col = alpha("grey", alpha = .3), pch = 20,
#      main = "Ancestry frequencies in both hybrid zones across SNPs",
#      xlab = "Mean African ancestry in California bees",
#      ylab = "Mean African ancestry in Argentina bees")


# 
# png("plots/density_MVN_overlay_on_scatterplot_w_ind_zone_quantiles.png", 
#     height = 8, width = 8, units = "in", res = 300)
# plot(meanA_CA, meanA_AR, col = alpha("grey", alpha = .3), pch = 20,
#      main = "Ancestry frequencies in both hybrid zones across SNPs",
#      xlab = "Mean African ancestry in California bees",
#      ylab = "Mean African ancestry in Argentina bees")
# #plot(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded, col = "grey", pch = 20)
# emdbook::HPDregionplot(cbind(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded), 
#                        prob = c(0.99, 0.95, 0.75, 0.5), 
#                        n = 100, # number of grid points
#                        #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
#                        col=c("lightgreen", "salmon", "lightblue", "yellow"), 
#                        lwd = 3, 
#                        add = TRUE)
# 
# abline(v = quantile(meanA_MVNsim_CA_zero_bounded, c(0.01, 0.99)), col = "lightgreen")
# abline(h = quantile(meanA_MVNsim_AR_zero_bounded, c(0.01, 0.99)), col = "lightgreen")
# legend("bottomright", legend = c(0.99, 0.95, 0.75, 0.5), 
#        title = "Neutral simulation (MVN) density",
#        lty = 1,
#        lwd = 4,
#        col=c("lightgreen", "salmon", "lightblue", "yellow"))
# dev.off()


# any genes in outlier regions?

# write a bed file with outlier regions (to compare with honeybee genes):
# ancestry calls are extended to halfway between any two calls and end at the position of the last/first call for a scaffold
# note: this does not extend all the way to the ends of the chromosomes (so some genes may not have ancestry calls)
sites_bed <- read.table("results/SNPs/combined_sept19/chr.var.sites.bed",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "start", "end", "snp_id"))
#table(paste(sites$scaffold, sites$pos, sep = "_") == sites_bed$snp) # match up, good
A_AR_CA <- sites %>%
  dplyr::select(chr, pos, cum_pos) %>%
  bind_cols(sites_bed, .) %>%
  mutate(CA = meanA_CA, AR = meanA_AR, combined = meanA) %>%
  # add FDRs
  mutate(FDR_CA_high = ifelse(CA >= FDRs[FDRs$FDR_values == .01, "CA_high"],
                              .01, ifelse(CA >= FDRs[FDRs$FDR_values == .05, "CA_high"],
                                          .05, ifelse(CA >= FDRs[FDRs$FDR_values == .1, "CA_high"],
                                                      .1, NA)))) %>%
  mutate(FDR_AR_high = ifelse(AR >= FDRs[FDRs$FDR_values == .01, "AR_high"],
                              .01, ifelse(AR >= FDRs[FDRs$FDR_values == .05, "AR_high"],
                                          .05, ifelse(AR >= FDRs[FDRs$FDR_values == .1, "AR_high"],
                                                      .1, NA)))) %>%
  mutate(FDR_CA_low = ifelse(CA <= FDRs[FDRs$FDR_values == .01, "CA_low"],
                             .01, ifelse(CA <= FDRs[FDRs$FDR_values == .05, "CA_low"],
                                         .05, ifelse(CA <= FDRs[FDRs$FDR_values == .1, "CA_low"],
                                                     .1, NA)))) %>%
  mutate(FDR_AR_low = ifelse(AR <= FDRs[FDRs$FDR_values == .01, "AR_low"],
                             .01, ifelse(AR <= FDRs[FDRs$FDR_values == .05, "AR_low"],
                                         .05, ifelse(AR <= FDRs[FDRs$FDR_values == .1, "AR_low"],
                                                     .1, NA)))) %>%
  mutate(FDR_shared_high = sapply(1:nrow(.), function(i) max(as.numeric(FDR_CA_high[i]), as.numeric(FDR_AR_high[i])))) %>%
  mutate(FDR_shared_low = sapply(1:nrow(.), function(i) max(as.numeric(FDR_CA_low[i]), as.numeric(FDR_AR_low[i])))) %>%
  dplyr::select(scaffold, start, end, snp_id, AR, CA, FDR_shared_high, FDR_AR_high, FDR_CA_high, 
                FDR_shared_low, FDR_AR_low, FDR_CA_low, combined, chr, pos, cum_pos)
  

# write bed file with mean ancestry for Argentina and California included bees
# also include whether a site meets a FDR threshold for selection
A_AR_CA %>%
  write.table(., 
              "results/mean_ancestry_AR_CA.bed", 
              sep = "\t", quote = F, col.names = F, row.names = F)
save(A_AR_CA, file = "results/mean_ancestry_AR_CA.RData")
#load("results/mean_ancestry_AR_CA.RData")


# # take away diagonal contribution of binomial sampling variance
# # the diagonal elements have variance due to drift + binomial sampling variance
# # = drift var + (n/n-1)*n*p*q/(n^2) where n/n-1 is the effect of small sample size n and npq is the standard binomial variance
# # p is the African ancestry allele frequency, q = 1-p and n = number of haplotypes = 2*n_bees in a population
# # but then because I'm using variance in frequencies, not counts, I divide by n^2
# # except this isn't quite right because it should be observed allele freq at a locus, not genomewide mean
# n_bees_per_pop <- sapply(names(zAnc_bees$alpha), function(p) meta.pop[meta.pop$population == p, "n_bees"])
# sampling_var_diag <- ((2*n_bees_per_pop/(2*n_bees_per_pop - 1))*n_bees_per_pop*2*zAnc_bees$alpha*(1 - zAnc_bees$alpha))/(2*n_bees_per_pop)^2 
# K_bees_minus_sampling <- zAnc_bees$K
# for (i in 1:length(sampling_var_diag)){
#   K_bees_minus_sampling[i,i] <- K_bees_minus_sampling[i,i] - sampling_var_diag[i]
# }
# K_bees_minus_sampling %>% # not quite
#   melt() %>%
#   dplyr::filter(!(Var1 == "Avalon_2014" | Var2 == "Avalon_2014")) %>%
#   ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_viridis(begin = 0, end = 1, direction = 1) +
#   #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
#   #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
#   ggtitle("K matrix - bees AR & CA")
# # simulate the variance:
# a = .5
# n_binom = 10
# 
# # simulate 100 loci, just binomial variance
# sim_pop1 <- rbinom(size = n_binom, n = n_sim, p = a)/n_binom  
# var(sim_pop1)
# mean(sim_pop1)
# cov(sim_pop1, sim_pop1)
# mean((sim_pop1 - mean(sim_pop1))^2) # observed variance
# mean(sim_pop1)*(1-mean(sim_pop1))/n_binom # theoretical no small sample size correction - looks correct
# mean(sim_pop1)*(1-mean(sim_pop1))/(n_binom-1) # theoretical w/ small sample size correction - incorrect
# # now add normal noise around allele freq due to drift before binomial sampling:
# sim_drift2 <- rnorm(n = n_sim, mean = a, sd = .1)
# sim_pop2 <- rbinom(size = n_binom, n = sim_drift2, p = a)/n_binom
# var(sim_drift2)
# mean(sim_pop2)
# var(sim_pop2) # observed variance
# mean((sim_pop2 - mean(sim_pop2))^2)
# mean(sim_pop2)*(1-mean(sim_pop2))/n_binom # mean theoretical
# 
# 
# log(1/(1 - diag(zAnc_bees$K)/(zAnc_bees$alpha*(1 - zAnc_bees$alpha))))
# log(diag(zAnc_bees$K)/(zAnc_bees$alpha*(1 - zAnc_bees$alpha)))


# r <- read.table("results/SNPs/thin1kb_common3/included_pos_on_Wallberg_Amel4.5_chr.r.bed",
#                 sep = "\t", stringsAsFactors = F) %>%
#   mutate(snp_id = read.table("results/SNPs/thin1kb_common3/included_pos_on_Wallberg_Amel4.5_chr.map",
#                              stringsAsFactors = F, header = F, sep = "\t")$V2) %>%
#   rename(chr = V1) %>%
#   rename(cM_Mb = V4) %>%
#   mutate(pos_chr = V2 + 1) %>% # because .bed files are 0 indexed
#   dplyr::select(c("snp_id", "chr", "pos_chr", "cM_Mb"))
# 
# # **** check on the recombination rate files!
# 
# d <- bind_cols(sites, A) %>%
#   left_join(., r, by = c("chr", "snp_id"))
# d$r_bin10 = with(d, cut(cM_Mb,  # note need to load map from scaffolds_to_chr.R
#                              breaks = unique(quantile(c(0, r$cM_Mb, 100), # I extend bounds so that every value gets in a bin
#                                                       p = seq(0, 1, by = .1))),
#                              right = T,
#                              include.lowest = T))
# table(d$r_bin10) # fewer SNPs represented in lowest recomb. bins
# head(d[is.na(d$r_bin10),]) # good, should be none
# 
# d$r_bin5 = with(d, cut(cM_Mb,  # note need to load map from scaffolds_to_chr.R
#                         breaks = unique(quantile(c(0, r$cM_Mb, 100), # I extend bounds so that every value gets in a bin
#                                                  p = seq(0, 1, by = .2))),
#                         right = T,
#                         include.lowest = T))
# table(d$r_bin5) # fewer SNPs represented in lowest recomb. bins
# table(d$r_bin5, is_outlier)
# head(d[is.na(d$r_bin5),]) # good, should be none
# 
# # make K matrix again with only high recombination rate region of genome
# 
# #a <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "(38.8,100]", ]))
# #a2 <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 %in% c("[0,8.49]","(8.39,15.6]"),]))
# a <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "(39.8,100]", ]))
# a2 <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "[0,14]", ]))
# 
# #a <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "(39.8,100]" & !is_outlier, ]))
# #a2 <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "[0,14]" & !is_outlier, ]))
# 
# # note: I need some way of summarising these K matrices statistically 
# # to test hypotheses about low vs. high r regions.
# 
# melt(a$K) %>%
#   mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
#   ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_viridis(begin = 0, end = 1, direction = 1) +
#   ggtitle("K matrix covariances r > 38cM/Mb - bees AR & CA no within-pop var")
# ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_highr.png", 
#        height = 6, width = 8, 
#        units = "in", device = "png")
# melt(cov2cor(a$K)) %>%
#   mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
#   ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_viridis(begin = 0, end = 1, direction = 1) +
#   ggtitle("K matrix correlations r > 38cM/Mb - bees AR & CA no within-pop var")
# ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_highr_correlations.png", 
#        height = 6, width = 8, 
#        units = "in", device = "png")
# melt(a2$K) %>%
#   mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
#   ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_viridis(begin = 0, end = 1, direction = 1) +
#   ggtitle("K matrix covariances r < 14cM/Mb - bees AR & CA no within-pop var")
# ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_lowr.png", 
#        height = 6, width = 8, 
#        units = "in", device = "png")
# melt(cov2cor(a2$K)) %>%
#   mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
#   ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_viridis(begin = 0, end = 1, direction = 1) +
#   ggtitle("K matrix correlations r < 14cM/Mb - bees AR & CA no within-pop var")
# ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_lowr_correlations.png", 
#        height = 6, width = 8, 
#        units = "in", device = "png")
# 
# 
# # a better null: maybe I should take the full covariance matrix for ind. ancestries, 
# # subtract off the diagonal for what I know the binomial sampling variance to be,
# # run a MVN, and then do binomial sampling at the end. 
# # Graham *might* say this is excessive since we don't know if some variance comes from the ancestry caller
# # but I think with the posteriors it's actually more conservative than a binomial
# # which is more likely to get to the bounds 0/1 than a posterior. 
# # (and I could use the ancestry call from ancestry_hmm if I wanted)
# # also if I have enough bees, this might just give me the same answer as the MVN
# 
# AbyPopr10 <- d %>%
#   gather("pop", "A", pops) %>%
#   group_by(., pop, r_bin10) %>%
#   summarise(meanA = mean(A))
# AbyPopr10 %>%
#   ggplot(aes(x = pop, y = meanA)) +
#   geom_point() +
#   facet_wrap(~r_bin10) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# AbyPopr5 <- d %>%
#   gather("pop", "A", pops) %>%
#   group_by(., pop, r_bin5) %>%
#   summarise(meanA = mean(A))
# AbyPopr5 %>%
#   ggplot(aes(x = pop, y = meanA, color = r_bin5)) +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))



# plot ancestry across all populations by mean latitude in that pop, separate by hybrid zone
# I used Riverside 2014 coordinates for Riverside 1999. 
# Note there are two Riverside collection sites, close in latitude, but different in temp etc. due to elevation
a <- cbind(A_AR_CA, A) %>%
  left_join(., dplyr::select(sites, c("chr", "pos", "chr_n", "snp_id")), by = "snp_id") %>%
  tidyr::gather(., "population", "A_ancestry", colnames(A)) %>%
  left_join(., meta.pop, by = "population")
a_mean <- a %>%
  group_by(population) %>%
  summarise(A_ancestry = mean(A_ancestry)) %>%
  left_join(meta.pop, by = "population") %>%
  mutate(abs_lat = abs(lat)) %>%
  mutate(abs_lat_c = abs_lat - mean(abs_lat))

# plot some ancestry clines across latitude for top SNPs:
# start with high shared outliers (there's only 15 regions):
high.shared.outliers
i = 1
a_shared_high_outlier_.05sig <- a %>%
  mutate(combined = (AR + CA)/2) %>%
  dplyr::select(combined, scaffold, pos, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_shared_high != "NA") %>%
  unique() %>%
  filter(scaffold == high.shared.outliers[i, "chr"] &
           pos >= high.shared.outliers[i, "start"] &
           pos <= high.shared.outliers[i, "end"]) %>%
  dplyr::arrange(combined) %>%
  .[1, "snp_id"]

# ID SNPs from very top outliers:
shared_high_top_snp <- a %>%
  mutate(combined = (AR + CA)/2) %>%
  dplyr::select(combined, scaffold, pos, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_shared_high != "NA") %>%
  unique() %>%
  dplyr::arrange(desc(combined)) %>%
  .[1, "snp_id"]
shared_low_top_snp <- a %>%
  mutate(combined = (AR + CA)/2) %>%
  dplyr::select(combined, scaffold, pos, FDR_shared_low, snp_id) %>%
  dplyr::filter(FDR_shared_low != "NA") %>%
  unique() %>%
  dplyr::arrange(combined) %>%
  .[1, "snp_id"]

# just AR
AR_high_only_top_snp <- a %>%
  dplyr::select(scaffold, pos, AR, CA, FDR_AR_high, FDR_CA_high, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_AR_high != "NA" & FDR_shared_high == "NA") %>%
  unique() %>%
  dplyr::arrange(desc(AR)) %>%
  .[1, "snp_id"]
AR_low_only_top_snp <- a %>%
  dplyr::select(scaffold, pos, AR, CA, FDR_AR_low, FDR_CA_low, FDR_shared_low, snp_id) %>%
  dplyr::filter(FDR_AR_low != "NA" & FDR_shared_low == "NA") %>%
  unique() %>%
  dplyr::arrange(AR) %>%
  .[1, "snp_id"]
# just CA
CA_high_only_top_snp <- a %>%
  dplyr::select(scaffold, pos, AR, CA, FDR_AR_high, FDR_CA_high, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_CA_high != "NA" & FDR_shared_high == "NA") %>%
  unique() %>%
  dplyr::arrange(desc(CA)) %>%
  .[1, "snp_id"]
# put together into data frame:
some_snp_outliers <- data.frame(
  snp_id = c(a_shared_high_outlier_.05sig,
             shared_high_top_snp,
             shared_low_top_snp,
             AR_high_only_top_snp,
             AR_low_only_top_snp,
             CA_high_only_top_snp),
  label = c("a 5% FDR outlier shared high A ancestry",
            "the top outlier shared high A ancestry",
            "the top outlier shared low A ancestry",
            "the top outlier Argentina-only high A ancestry",
            "the top outlier Argentina-only low A ancestry",
            "the top outlier California-only high A ancestry"),
  label_short = c("High A: shared (5% FDR outlier)",
                  "High A: shared",
                  "Low A: shared",
                  "High A: S. America",
                  "Low A: S. America",
                  "High A: N. America"),
  name = c("a_shared_high_outlier_.05sig",
           "shared_high_top_snp",
           "shared_low_top_snp",
           "AR_high_only_top_snp",
           "AR_low_only_top_snp",
           "CA_high_only_top_snp"))
  
# make plots:
for (i in 1:nrow(some_snp_outliers)){
my_plot <- a %>%
    filter(snp_id == some_snp_outliers[i, "snp_id"]) %>%
    ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
    geom_point() +
    geom_point(data = a_mean, shape = 2) +
    ggtitle(paste0("cline across latitude for ", some_snp_outliers[i, "label"], ": ", some_snp_outliers[i, "snp_id"])) +
    xlab("absolute latitude") +
    ylab("population mean African ancestry") +
    facet_grid(zone ~ .)
ggsave(paste0("plots/outlier_clines_", some_snp_outliers[i, "name"], ".png"),
       plot = my_plot,
       height = 5, 
       width = 10, 
       units = "in", 
       device = "png")
}
# plot all outliers in a grid
a %>%
  filter(snp_id %in% some_snp_outliers$snp_id) %>%
  left_join(., some_snp_outliers, by = "snp_id") %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
  geom_point() +
  geom_point(data = a_mean, shape = 2, color = "grey") +
  #ggtitle(paste0("cline across latitude for ", some_snp_outliers[i, "label"], ": ", some_snp_outliers[i, "snp_id"])) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  facet_grid(zone ~ label_short) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("plots/ancestry_clines_at_top_outlier_snps_all_pops_incl_mx.png",
       width = 8, height = 5, units = "in",
       device = "png")
# only plot the populations that are included (CA/AR 2014+2018),
# at the top outliers of each type:
a %>%
  filter(., population %in% c(AR_pops_included$population, CA_pops_included$population)) %>%
  filter(snp_id %in% some_snp_outliers$snp_id) %>%
  left_join(., some_snp_outliers, by = "snp_id") %>%
  filter(name != "a_shared_high_outlier_.05sig") %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
  geom_point() +
  geom_point(data = filter(a_mean, population %in% c(AR_pops_included$population, CA_pops_included$population)),
             color = "grey", 
             shape = 2) +
  #ggtitle(paste0("cline across latitude for ", some_snp_outliers[i, "label"], ": ", some_snp_outliers[i, "snp_id"])) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  facet_grid(zone ~ label_short) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("plots/ancestry_clines_at_top_outlier_snps.png",
       width = 8, height = 5, units = "in",
       device = "png")
ggsave("../../bee_manuscript/figures/ancestry_clines_at_top_outlier_snps.png",
       width = 8, height = 5, units = "in",
       device = "png")

# new approach: get 1 SNP per outlier. So find top SNP in each outlier region:
# outlier_sets load outlier sets from rank_genes.R
# and outlier_setnames too
# here I have all the outliers of each type:
top_SNP_outliers_all
top_SNP_outliers
outlier_set_names
top_SNP_outlier_labels = c("High A: shared",
                "Low A: shared",
                "High A: S. America",
                "Low A: S. America",
                "High A: N. America")
for (i in 1:length(top_SNP_outlier_labels)){
a %>%
  filter(., population %in% c(AR_pops_included$population, CA_pops_included$population)) %>%
  filter(snp_id %in% top_SNP_outliers[[i]]$top_snp) %>%
  left_join(., top_SNP_outliers[[i]], by=c("snp_id"="top_snp", "scaffold")) %>%
    mutate(snp = paste(substr(scaffold, start = 6, stop = 100), start+1, sep = ":")) %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = factor(FDR_region))) + #color = population, shape = factor(FDR_region))) +
  geom_point() +
  geom_point(data = filter(a_mean, population %in% c(AR_pops_included$population, CA_pops_included$population)),
             color = "grey", 
             shape = 2) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  ggtitle(top_SNP_outlier_labels[i]) +
  facet_grid(zone ~ snp) +
  theme_bw() +
  #scale_shape_manual(values = c(20, 1))+
  scale_color_manual(values = c("red", "orange")) +
  theme(legend.position = "bottom") +
  #guides(color = F) +
  #labs(shape = "FDR") +
  labs(color = "FDR")
  ggsave(paste0("plots/ancestry_clines_at_all_outliers_top_snp_", outlier_set_names[i],".png"),
         width = 16, height = 5, units = "in",
         device = "png")
  ggsave(paste0("../../bee_manuscript/figures/ancestry_clines_at_all_outliers_top_snp_", outlier_set_names[i],".png"),
         width = 16, height = 5, units = "in",
         device = "png")
  }
