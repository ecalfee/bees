library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)

# get metadata for individuals included in NGSadmix analysis
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")

# get ID's for pass1 bams
IDs <- read.table("../bee_samples_listed/pass1.list", stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") 

# scaffold lengths
genome <- read.table("../data/honeybee_genome/ordered_scaffolds.lengths",
                     stringsAsFactors = F, sep = "\t", header = F)
colnames(genome) <- c("scaffold", "length")
genome <- genome %>%
  separate(., "scaffold", into = c("chr", "n_scaffold"), 
           remove = F, sep = "[.]")
# view lengths of all chromosomes:
filter(genome, !is.na(n_scaffold)) %>% 
  group_by(., chr) %>% 
  summarize(., chr_length = sum(length))

# depth counts from ANGSD
chr_name = "Group3"
#chr_name = "Group2"
max_depth_bin = 20000
length_Group1.1 = genome[genome$scaffold == "Group1.1", "length"]# 1382403
length_chr = sum(genome[genome$chr == chr_name, "length"])
tot_length = length_chr # length of subset/chrom/scaffold
sampleDepth <- read.table(paste0("results/depth/", chr_name, ".depthSample"),
                          stringsAsFactors = F, header = F, sep = "\t")
globalDepth <- t(read.table(paste0("results/depth/", chr_name, ".depthGlobal"),
                            stringsAsFactors = F, header = F, sep = "\t")[1,-(max_depth_bin+2)])
globalBins <- 0:max_depth_bin
sampleBins <- 0:100

# set to path to mapping metrics file
raw_metrics_file <- "metrics/pass1.all.metrics.raw.Nreads" 
raw_metrics <- read.table(raw_metrics_file, header = T,
                      sep = "\t", stringsAsFactors = F)
Q30_read_count <- read.table("metrics/pass1.all.metrics.Q30.Nreads",
                             sep = "\t", header = F,
                             stringsAsFactors = F)
colnames(Q30_read_count) <- c("Bee_ID", "n_reads_pass")

# count number of spots with zero coverage:
globalDepth[1] <- tot_length - sum(globalDepth)

# plot global depth
plot(globalBins[-1], globalDepth[-1],
     main = paste0("non-zero total depth pass1: ", chr_name),
     cex = .1, xlab = "global depth", ylab = "# sites")

# get full distribution global depth
globalDepthVals <- unlist(mapply(rep, times = globalDepth, x = globalBins))
summary(globalDepthVals) # for some reason omits spots with depth 0
# what % of sites have depth between .5x and 2x where x is the global mean depth?
bounds <- summary(globalDepthVals)["Mean"]*c(.5, 2)
sum(globalDepthVals>bounds[1] & globalDepthVals<bounds[2])/length(globalDepthVals)
round(bounds)

# function to calculate mean depth from raw depth dataframe
meanDepth = function(dep, bins = sampleBins, length){
  dep[ , dim(dep)[2]] <- NULL
  colnames(dep)<- paste0("n", bins)
  dep[ , 1] = length - apply(dep, 1, sum) # count spots uncounted (zero depth)
  meanDepth <- apply(dep[ , paste0("n", bins)], 1,
                     function(r) r %*% bins/sum(r))
  return(meanDepth)
}
# calculate mean depth for each sample
bees$mean_depth = meanDepth(sampleDepth, length = tot_length)
hist(bees$mean_depth, breaks = 20)
ggplot(data = bees) +
  geom_histogram(aes(x = mean_depth))
ggplot(data = bees) +
  geom_histogram(aes(x = mean_depth,
                     fill = group)) + 
  facet_wrap(~source) +
  ggtitle(paste0("sample depth scaffolds ", chr_name))
ggsave(paste0("plots/sample_depth_by_group_", chr_name, ".png"), 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

summary(bees$mean_depth)
tapply(bees$mean_depth, bees$group, summary)
filter(bees, source == "Calfee") %>%
  select(mean_depth) %>%
  summary(.) # mean depth of new bees after quality filtering and de-duplication
sd(bees[bees$source == "Calfee", "mean_depth"])


# plot mapping statistics across bees: % mapped; % mapped > mapQ30:
# add in alignment metrics: mapping quality and sequencing coverage
metrics <- raw_metrics %>%
  mutate(n_reads_raw = 2*READ_PAIRS_EXAMINED + UNPAIRED_READS_EXAMINED) %>%
  mutate(prop_unmapped = UNMAPPED_READS/n_reads_raw) %>%
  mutate(n_reads_dedup = 2*(READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES) +
           UNPAIRED_READS_EXAMINED - UNPAIRED_READ_DUPLICATES) %>%
  mutate(Bee_ID = StudyID) %>%
  left_join(., bees, by = "Bee_ID") %>%
  left_join(., Q30_read_count, by = "Bee_ID") %>%
  # different sequencing sources have diff. lengths
  mutate(read_length = ifelse(source %in% c("Calfee", "Ramirez"),
                              150, NA)) %>%
  mutate(prop_filtTot = (n_reads_raw - n_reads_pass)/n_reads_raw) %>%
  mutate(n_reads_filtOutQ30 = n_reads_dedup-n_reads_pass) %>%
  mutate(prop_filtOutQ30 = n_reads_filtOutQ30/n_reads_dedup) %>%
  mutate(n_reads_dropped = n_reads_raw - n_reads_dedup + n_reads_filtOutQ30)


# plots for coverage and filtering effects
png(paste("plots/pass1_total_reads_dropped_by_group.png"),
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtTot, data=metrics, geom="density", fill=strain, alpha=I(.5),
      main="Prop. reads filtered total by group", xlab="Prop. reads dropped",
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
# same but in scatterplot across sequencing coverage
ggplot(data = metrics) + geom_point(aes(x = n_reads_raw,
                                            y = prop_filtTot,
                                            color = strain)) +
  theme(legend.title=element_blank())
png(paste("plots/pass1_reads_unmapped_by_group.png"),
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_unmapped, data=metrics, geom="density", fill=strain, alpha=I(.5),
      main="Prop. unmapped reads by group", xlab="Prop. reads unmapped",
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
ggplot(data = metrics) + geom_point(aes(x = n_reads_raw,
                                            y = prop_unmapped,
                                            color = group)) +
  theme(legend.title=element_blank()) +
  facet_wrap(~source) + 
  ggtitle("prop. unmapped reads by group")
ggsave("plots/pass1_prop_reads_unmapped_by_group.png", 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)

# percent read duplicates - picard
ggplot(data = metrics) + 
  geom_point(aes(x = n_reads_raw,
                 y = PERCENT_DUPLICATION,
                 color = group)) +
  theme(legend.title=element_blank()) +
  facet_wrap(~source) + 
  ggtitle("prop. duplicated reads by group")
ggsave("plots/pass1_prop_reads_duplicated_by_group.png", 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)
dev.off()
# raw depth from # reads sequenced (before filtering)
filter(metrics, source == "Calfee") %>%
  ggplot(.) +
  geom_histogram(aes(x = n_reads_raw*150/237000000,
                     # note: ref. panel sequenced w/ shorter reads
                     fill = group)) + 
  facet_wrap(~source) +
  ggtitle("est. depth from tot. raw reads")
ggsave("plots/pass1_depth_new_bees_est_tot_raw_reads.png", 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)
dev.off()
summary(metrics[metrics$source == "Calfee",
                      "n_reads_raw"]*150/237000000)

# depth for reads passing quality filters Q30:
# raw depth from # reads sequenced (before filtering)
filter(metrics, source == "Calfee") %>%
  ggplot(.) +
  geom_histogram(aes(x = n_reads_pass*150/237000000,
                     # note: ref. panel sequenced w/ shorter reads
                     fill = group)) + 
  ggtitle("est. depth from tot. passing reads")
ggsave("plots/pass1_depth_new_bees_est_tot_passing_reads.png", 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)
dev.off()
summary(metrics[metrics$source == "Calfee",
                    "n_reads_pass"]*150/237000000)

# are the high depth bees the ones with separate PCR?
# (because their libraries failed initially)
filter(metrics, source == "Calfee") %>%
  mutate(., PCR_group = ifelse(Bee_ID %in% 
                           c("AR1410", "AR2603", "AR0804", "AR0311", "AR1303"),
                            "group2", "group1")) %>%
  ggplot(.) +
  geom_histogram(aes(x = n_reads_pass*150/237000000,
                     # note: ref. panel sequenced w/ shorter reads
                     fill = PCR_group)) +
  ggtitle("est. depth from tot. passing reads")
ggsave("plots/pass1_depth_new_bees_by_PCR_Group.png", 
       device = png(), 
       width = 15, height = 8, units = "in",
       dpi = 200)
dev.off()

