library(dplyr)
library(tidyr)
library(MASS) # for mvrnorm
source("k_matrix.R")
source("calc_FDRs.R") # functions to calculate false discovery rates

# script to combined data and summarise ancestry_hmm results for later plotting,
# simulate different null models,
# and calculate false discovery rates 

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


# write a bed file with outlier regions (to compare with honeybee genes):
# ancestry calls are extended to halfway between any two calls and end at the position of the last/first call for a scaffold
# note: this does not extend all the way to the ends of the chromosomes (so some genes may not have ancestry calls)
sites_bed <- read.table("results/SNPs/combined_sept19/chr.var.sites.bed", # created by ancestry_sites_to_tracts.R
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


