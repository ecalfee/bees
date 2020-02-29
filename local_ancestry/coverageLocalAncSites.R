# what is the distribution of coverage across sites used for local ancestry inference?
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/AR0302.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/CA0502.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/AR1403.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/CA1120.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/AR2603.counts.txt", header = F, stringsAsFactors = F)
tot <- apply(a, 1, sum)
mean(tot)
hist(rpois(n=length(tot), lambda = mean(tot)), freq = F, col = "pink")
hist(rnbinom(n = length(tot), size = mean(tot)/(3 - 1), mu = mean(tot)), freq = F,
          col = "yellow", add = T)
hist(tot, freq = F, add = T, col = "blue")
hist(rnbinom(n = length(tot), size = mean(tot)/(2 - 1), mu = mean(tot)), freq = F,
     col = "green", add = T, bins = 20)
mean(tot)
var(tot) # true variance in coverage
mean(tot) + mean(tot)*(3-1) # expected variance with VMR=3
var(rnbinom(n = length(tot), size = mean(tot)/(3 - 1), mu = mean(tot))) # observed variance
mean(tot) + mean(tot)*(2-1) # expected variance with VMR=2
var(rnbinom(n = length(tot), size = mean(tot)/(2 - 1), mu = mean(tot))) # observed variance

# so it looks like at least for each of this small set of bees with varying coverage from ~1 to ~11x mean coverage,
# their distribution of coverage across included SNPs has mean variance and distribution similar to a negative binomial with VMR=2 
# (slightly less variable that what I modeled for neg binom but more variable than a poisson distribution of coverage)

# load metadata for all bees:
load("results/meta.RData")

mean_cov <- sapply(meta.ind$Bee_ID, function(b) # takes several minutes
  mean(apply(read.table(paste0("results/SNPs/combined_sept19/countsMajMin/", b, ".counts.txt"), 
                        header = F, stringsAsFactors = F), 1, sum)))
# save
data.frame(meta.ind$Bee_ID, hmm_coverage = mean_cov) %>%
  write.table("results/mean_depth_of_coverage_across_ancestry_hmm_SNPs.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
summary(mean_cov)
hist(mean_cov)
meta.ind %>%
  mutate(hmm_coverage = mean_cov) %>%
  ggplot(., aes(x = hmm_coverage, fill = lat > 0)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ year) +
  scale_x_continuous(breaks = seq(0, 26, by = 2)) +
  theme_light()

meta.ind %>% # shows median, upper and lower quartiles
  mutate(hmm_coverage = mean_cov) %>%
  ggplot(., aes(y = hmm_coverage, group = factor(year), x = factor(year))) +
  geom_boxplot() +
  theme_classic() +
  ylab("Mean depth individual bees") +
  xlab("Sample")
ggsave("plots/boxplot_mean_depth_ind_bees_hmm_SNPs.png",
       height = 4, width = 4, units = "in", dpi = 600, device = "png")



meta.ind %>%
  mutate(hmm_coverage = mean_cov) %>%
  group_by(year) %>%
  summarise(min = min(hmm_coverage),
            lower_quartile = quantile(hmm_coverage, 0.25),
            median = median(hmm_coverage),
            mean = mean(hmm_coverage),
            upper_quartile = quantile(hmm_coverage, 0.75),
            max = max(hmm_coverage)) %>%
  write.table("results/summary_mean_depth_of_coverage_across_ancestry_hmm_SNPs.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
