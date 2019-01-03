library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)

# ID and population labels
pass1_refBees <- read.table("../bee_samples_listed/post_1994.meta", 
                   stringsAsFactors = F, header = T, sep = "\t")
pass1_newBees_ID <- read.table("../bee_samples_listed/pass1.list",
                    stringsAsFactors = F, header = F, sep = "\t")[99:161, ]
pass1_newBees = data.frame(ID = pass1_newBees_ID, geographic_location = substr(pass1_newBees_ID, 1, 2),
                           population = substr(pass1_newBees_ID, 1, 4), strain = "admixed", year = 2018)
meta <- bind_rows(pass1_refBees, pass1_newBees)

# depth counts from ANGSD
sampleDepth <- read.table("results/depth/Group1.1.depthSample",
                          stringsAsFactors = F, header = F, sep = "\t")
globalDepth <- t(read.table("results/depth/Group1.1.depthGlobal",
                          stringsAsFactors = F, header = F, sep = "\t")[1,-10002])
globalBins <- 0:10000
sampleBins <- 0:100
tot_length = 1382403 # length of Group1.1 scaffold

mapping_metrics_file <- NA # set to path to mapping metrics file

# count number of spots with zero coverage:
globalDepth[1] <- tot_length - sum(globalDepth)

# plot global depth
plot(globalBins, globalDepth, 
     main = "depth all pass1 ind's across scaffold Group1.1", 
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
meta$mean_depth = meanDepth(sampleDepth, length = tot_length)
hist(meta$mean_depth, breaks = 20)
ggplot(data = meta) +
  geom_histogram(aes(x = mean_depth))
ggplot(data = meta) +
  geom_histogram(aes(x = mean_depth, 
                     fill = strain))
#plot(pass1$est_coverage, mean_depthF, main = "Depth est. regions vs. # reads") # perfectly correlated
summary(meta$mean_depth)
filter(meta, strain == "admixed") %>%
  select(mean_depth) %>%
  summary(.) # mean depth of new bees after quality filtering and de-duplication
sd(meta[meta$strain == "admixed", "mean_depth"])


# plot mapping statistics across hilo: % mapped; % mapped > mapQ30:
# add in alignment metrics: mapping quality and sequencing coverage
raw_metrics_no_header <- read.table(mapping_metrics_file, 
                                    stringsAsFactors = F, header = F, skip = 1)
con = file(mapping_metrics_file, "r")
metrics_head <- strsplit(readLines(con, 1), split = "\t")[[1]]
close(con)
raw_metrics0 <- raw_metrics_no_header %>% 
  unite(., "LIB", V2, V3, sep= "_") # join 2 library_unknown columns
colnames(raw_metrics0) <- metrics_head # add header
raw_metrics <- raw_metrics0 %>%
  mutate(n_reads = 2*READ_PAIRS_EXAMINED + UNPAIRED_READS_EXAMINED) %>%
  mutate(prop_unmapped = UNMAPPED_READS/n_reads) %>%
  mutate(n_dedup = 2*(READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES) +
           UNPAIRED_READS_EXAMINED - UNPAIRED_READ_DUPLICATES) %>%
  mutate(n = as.integer(substr(StudyID, 5, 10))) %>%
  left_join(., pass1, by = "n") %>%
  mutate(n_filtOutQ30 = n_dedup-num_read_pass) %>%
  mutate(prop_filtOutQ30 = n_filtOutQ30/n_dedup) %>%
  mutate(prop_filtTot = (n_reads - num_read_pass)/n_reads) %>%
  mutate(n_dropped = n_reads-n_dedup+n_filtOutQ30)



# plots for coverage and filtering effects
png(paste("plots/pass1_total_reads_dropped_by_group.png"),
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtTot, data=raw_metrics, geom="density", fill=strain, alpha=I(.5), 
      main="Prop. reads filtered total by group", xlab="Prop. reads dropped", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
# same but in scatterplot across sequencing coverage
ggplot(data = raw_metrics) + geom_point(aes(x = n_reads, 
                                            y = prop_filtTot,
                                            color = strain)) +
  theme(legend.title=element_blank())
png(paste("plots/pass1_reads_unmapped_by_group.png"),
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_unmapped, data=raw_metrics, geom="density", fill=strain, alpha=I(.5), 
      main="Prop. unmapped reads by group", xlab="Prop. reads unmapped", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
ggplot(data = raw_metrics) + geom_point(aes(x = n_reads, 
                                            y = prop_unmapped,
                                            color = strain)) +
  theme(legend.title=element_blank())
# of reads that mapped and are not duplicates -- what proportion are mapQ <30?
png(paste("plots/pass1_prop_reads_mapped_but_belowQ30_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtOutQ30, data=raw_metrics, geom="density", fill=strain, alpha=I(.5), 
      main="Prop. reads filtered by mapQ<30 by group", xlab="Prop. reads < Q30", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()

png(paste("plots/pass1_prop_duplication_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(PERCENT_DUPLICATION, data=raw_metrics, geom="density", fill=strain, alpha=I(.5), 
      main="Prop. reads duplications by group", xlab="Prop. reads duplicates", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()

