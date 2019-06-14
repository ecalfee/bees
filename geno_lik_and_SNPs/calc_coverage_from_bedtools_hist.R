library(dplyr)
library(ggplot2)
file_in <- "coverage/backup_pass1_plus_kohn_and_wallberg_random_100_regions.coverage"
coverage <- read.table(file_in, stringsAsFactors = F)
colnames(coverage) <- c("scaffold", "start", "end", "x_coverage", "n_bp", "region_length", "perc_bp")
by_region <- coverage %>%
  group_by(paste(scaffold, start, end, sep = "_")) %>%
  summarize(mean_coverage = sum(x_coverage*perc_bp)) 
summary(by_region$mean_coverage)
hist(by_region$mean_coverage)
coverage %>%
  filter(x_coverage <= 900) %>%
  summarize(count_low = sum(n_bp))/sum(coverage$n_bp)
coverage %>%
  filter(x_coverage >= 4000)%>%
  summarize(count_high = sum(n_bp))/sum(coverage$n_bp)
summary(coverage$x_coverage)

# get total reads and then multiply by different read lengths for different studies:
reads <- read.table("../filtered_bams/metrics/Q30/CA_AR_MX_harpur_sheppard_kohn_wallberg.nreads",
                    stringsAsFactors = F, header = F) 
colnames(reads) <- "nreads_Q30"
IDs <- read.table("../bee_samples_listed/CA_AR_MX_harpur_sheppard_kohn_wallberg.list", stringsAsFactors = F)
colnames(IDs) <- c("Bee_ID")
# get metadata for individuals included in NGSadmix analysis
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                    header = T, sep = "\t")

# calculate coverage
ref_genome_size <- sum(read.table("../data/honeybee_genome/Amel_4.5_scaffolds.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including GroupUn scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
# slight overestimate of genome size due to including gaps, but ok
read_lengths <- data.frame(source = c("Calfee", "Ramirez", "Harpur", "Wallberg", "Kohn", "Sheppard"),
                           read_lengths = c(150, 100, 50, 75, 100, 100),
                           stringsAsFactors = F)

# join all together
bees <- bind_cols(IDs, reads) %>%
  dplyr::left_join(., meta, by = "Bee_ID") %>% # add in source for some of my bees missing from meta1 (to fix later)
  mutate(source = ifelse(is.na(source) & substr(Bee_ID, 1, 2) %in% c("CA", "AR", "MX"), "Calfee", source)) %>%
  left_join(., read_lengths, by = "source") %>%
  mutate(est_coverage = read_lengths*nreads_Q30/ref_genome_size)

# plot coverage across different genome sources
ggplot(bees, aes(x = est_coverage, fill = source)) + geom_histogram(bins = 100)

# mean total coverage
sum(bees$est_coverage)

bees %>% 
  filter(!(geographic_location == "Mexico" | source %in% c("Kohn", "Wallberg"))) %>%
           dplyr::select(est_coverage) %>%
  sum()

# write file with estimated coverage per bee
dplyr::select(bees, c("Bee_ID", "read_lengths", "nreads_Q30", "est_coverage")) %>%
  write.table(., "results/CA_AR_MX_harpur_sheppard_kohn_wallberg/coverage/est_coverage_by_reads_Q30.txt", quote = F, col.names = T, row.names = F, sep = "\t")


