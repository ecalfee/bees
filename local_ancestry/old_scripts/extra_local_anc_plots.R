library(dplyr)
library(tidyr)
source("k_matrix.R")
# script to combined data and summarise ancestry_hmm results
# for later plotting

# load sample data
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
# meanA0 <- apply(A, 1, mean)
# meanC0 <- apply(C, 1, mean)
# meanM0 <- apply(M, 1, mean)


CA_pops_included <- meta.pop[meta.pop$zone == "N. America", ]
AR_pops_included <- meta.pop[meta.pop$zone == "S. America", ]
CA_A <- A[, CA_pops_included$population]
AR_A <- A[, AR_pops_included$population]

# take mean across individuals, not across populations:
meanA_CA <- apply(CA_A, 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
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

# K matrix
zAnc_bees = make_K_calcs(t(A[ , pops_by_lat]))

save(file = "results/zAnc.RData", list = "zAnc_bees")


# what are mean covariances/correlations?
# within CA
# within AR S
# within AR N
# CA -> AR S
# CA -> AR N
# AR S -> AR N

# mean ancestry covariances
mean_cov_k <- get_mean_from_K(zAnc_bees$K)
#mean_cov_k
mean_cov_k %>%
  write.table(., "results/mean_anc_cov_grouped.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")

# mean ancestry correlations
mean_corr_k <- get_mean_from_K(cov2cor(zAnc_bees$K))
#mean_corr_k
mean_corr_k %>%
  write.table(., "results/mean_anc_corr_grouped.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")

# Simulated ancestries
# simulate from MVN distribution:
n_sim <- 10^5
set.seed(101)
MVNsim <- MASS::mvrnorm(n = n_sim, 
                        mu = zAnc_bees$alpha, 
                        Sigma = zAnc_bees$K, 
                        tol = 1e-6, 
                        empirical = FALSE, 
                        EISPACK = FALSE)

# how many simulations exceed reasonable bounds [0,1]?
sum(MVNsim < 0)/sum(table(MVNsim<0)) # low outliers need to be set to bound
sum(MVNsim > 1)/sum(table(MVNsim>1))
# by population
perc_MVNsim_out_of_bounds = data.frame(population = colnames(MVNsim),
                                       Lower = apply(MVNsim, 2, function(x) sum(x < 0))/nrow(MVNsim),
                                       Upper = apply(MVNsim, 2, function(x) sum(x > 1))/nrow(MVNsim),
                                       stringsAsFactors = F)

# Now set bounds at 0 and 1
MVNsim_bounded <- data.frame(MVNsim, stringsAsFactors = F) 
MVNsim_bounded[MVNsim < 0] <- 0
MVNsim_bounded[MVNsim > 1] <- 1


MVNsim_AR_bounded <- MVNsim_bounded[ , AR_pops_included$population]
MVNsim_CA_bounded <- MVNsim_bounded[ , CA_pops_included$population]
# mean across populations
#meanA_MVNsim_AR_bounded0 <- apply(MVNsim_AR_bounded, 1, mean)
#meanA_MVNsim_CA_bounded0 <- apply(MVNsim_CA_bounded, 1, mean)
# mean across individuals
meanA_MVNsim_AR_bounded <- apply(MVNsim_AR_bounded[ , AR_pops_included$population], 1, 
                                 function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_bounded <- apply(MVNsim_CA_bounded[ , CA_pops_included$population], 1, 
                                 function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))

# combined mean across individuals
meanA_MVNsim_bounded <- apply(cbind(MVNsim_AR_bounded, MVNsim_CA_bounded)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                              function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))
#meanA_MVNsim <- apply(cbind(MVNsim_AR, MVNsim_CA)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
#                      function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))

# save MVN simulation as data objects
save(MVNsim_bounded, meanA_MVNsim_bounded, 
     meanA_MVNsim_AR_bounded, meanA_MVNsim_CA_bounded,
     file = "results/MVNsim_bounded.RData")

# save unbounded MVN results too for plots on effect of [0,1] truncation
MVNsim_AR <- MVNsim[ , AR_pops_included$population]
MVNsim_CA <- MVNsim[ , CA_pops_included$population]
meanA_MVNsim_AR_unbounded <- apply(MVNsim_AR[ , AR_pops_included$population], 1, 
                                   function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_unbounded <- apply(MVNsim_CA[ , CA_pops_included$population], 1, 
                                   function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))

save(meanA_MVNsim_AR_unbounded, meanA_MVNsim_CA_unbounded,
     perc_MVNsim_out_of_bounds,
     file = "results/MVNsim_unbounded.RData")



# Binomial sampling (no within or between population covariances)
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
# save
save(PoiBinsim_CA, PoiBinsim_AR, PoiBinsim_combined,
     file = "results/PoiBinsim.RData")


# MVN simulation bounded with no covariances
set.seed(101)
MVNsim_no_cov <- MASS::mvrnorm(n = n_sim, 
                               mu = zAnc_bees$alpha, 
                               Sigma = diag(diag(zAnc_bees$K), length(zAnc_bees$alpha), length(zAnc_bees$alpha)), 
                               tol = 1e-6, 
                               empirical = FALSE, 
                               EISPACK = FALSE)

# how many simulations exceed reasonable bounds [0,1]?
sum(MVNsim_no_cov < 0)/sum(table(MVNsim_no_cov<0)) # low outliers need to be set to bound
sum(MVNsim_no_cov > 1)/sum(table(MVNsim_no_cov>1))

# truncate at [0, 1]
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
# save
save(MVNsim_no_cov_bounded, meanA_MVNsim_AR_no_cov_bounded, meanA_MVNsim_CA_no_cov_bounded, meanA_MVNsim_no_cov_bounded,
     file = "results/MVNsim_no_cov_bounded.RData")

# Calculate false discovery rates based on MVN simulations
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

FDRs = data.frame(FDR_values = FDR_values,
                  CA_high = FDRs_CA_high2,
                  CA_low = FDRs_CA_low2,
                  AR_high = FDRs_AR_high2,
                  AR_low = FDRs_AR_low2,
                  stringsAsFactors = F)
write.table(FDRs, "results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", quote = F, col.names = T, row.names = F, sep = "\t")

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
