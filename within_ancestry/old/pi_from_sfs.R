# scratch analyses
# taken out of plot_pi_fst_from_allele_freqs.R
###---old SFS----####
# read in the SFS's that were successfully calculated (for comparison):
path_sfs_folded <- "results/folded_SFS_incorrect/non-outlier_regions/"
SFS <- lapply(pops, function(x) 
  read.table(paste0(path_sfs_folded, "/combined/", x, ".folded.sfs"),
             stringsAsFactors = F))
calc_pi_sfs_folded <- function(sfs){
  n <- 0:(length(sfs) - 1) # number of copies observed
  p <- n/((length(sfs) - 1)*2) # allele frequency
  f <- sfs/sum(sfs) # observation frequency of sfs bin
  het <- 2*p*(1-p)
  sum(het*f) # take weighted mean of pi across frequency bins
}

sfs1 <- unlist(read.table(paste0(path_sfs_folded, "/combined/AR01", ".folded.sfs"),
                          stringsAsFactors = F))
calc_pi_sfs_folded(sfs1)
lapply(SFS, calc_pi_sfs_folded)
thetas <- do.call(rbind,
                  lapply(ancestries_combined, function(a)
                    data.frame(population = pops,
                               ancestry = a,
                               theta = unlist(lapply(pops, function(x)
                                 calc_pi_sfs_folded(read.table(paste0(path_sfs_folded, "/", a, "/", x, ".folded.sfs"),
                                                               stringsAsFactors = F)))))))
# something's not quite right. populations with less positions with data have higher diversity estimates:
plot(lapply(SFS, sum), thetas[thetas$ancestry=="combined", "theta"], main = "neg. relationship ANGSD estimated pi and number of sites in SFS",
     xlab = "total sites in SFS", ylab = "theta estimate (combined)")

thetas_ACM <- data.frame(population = ACM,
                         ancestry = "combined",
                         theta = unlist(lapply(ACM, function(x)
                           calc_pi_sfs_folded(read.table(paste0(path_sfs_folded, "/", "combined", "/", x, ".folded.sfs"),
                                                         stringsAsFactors = F))))) %>%
  mutate(ancestry = paste0(population, population))

thetas_test <- data.frame(population = pops,
                          theta = unlist(lapply(SFS, calc_pi_sfs_folded)),
                          ancestry = "combined")

thetas %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = abs(lat), y = theta, color = ancestry, shape = factor(year))) +
  geom_point() +
  ggtitle("Allelic diversity at known SNPs within ancestries and combined across all ancestries") +
  xlab("Degrees latitude from the equator") +
  facet_grid(. ~ zone, scales = "free_x") +
  geom_abline(data = thetas_ACM, aes(intercept = theta, slope = 0, color = ancestry))
ggsave("plots/pi_by_latitude_from_folded_SFS.png", device = "png",
       width = 10, height = 5)
ggsave("../../bee_manuscript/figures/pi_by_latitude_from_folded_SFS.png", device = "png",
       width = 10, height = 5)

####--END old SFS---##########################333
