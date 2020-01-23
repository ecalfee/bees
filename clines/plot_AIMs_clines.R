library(dplyr)
library(tidyr)
library(ggplot2)
  
# plot clines in allele frequencies at outlier loci
# to see snps underlying ancestry skews and confirm ancestry calls are consistent wiht allelic clines
source("../colors.R") # get color palette
load("../local_ancestry/results/pops_by_lat.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 

# first load SNPs with ancestry calls:
load("../local_ancestry/results/sites_r.RData") # loads sites dataframe
# and find the outlier loci.

# load allele freq data for each population at all AIMs
sites_11 <- read.table("results/AIMs/M/Group11.ACM.freqs", header = T)
# then for chr11 find nearest M AIMs for outlier region & plot the allele freq clines NA and SA
mafs_11_list <- lapply(unique(meta.pop$population), function(pop) # for all sites, get pop freq
  left_join(sites_11[ , c("scaffold", "pos", "major", "minor")], 
                  read.table(paste0("results/AIMs/M/Group11/", pop, ".mafs.gz"), header = T),
                  by = c("scaffold"="chromo", "pos"="position", "major", "minor")) %>%
  dplyr::select(scaffold, pos, phat, nInd) %>%
  pivot_longer(data = ., cols = c("phat", "nInd"), names_to = "type", values_to = pop))
mafs_11 <- cbind(left_join(mafs_11_list[[1]][ , c("scaffold", "pos", "type")],
                 sites_11,
                 by = c("scaffold", "pos")),
                 # add in metadata for each site, then grab pop freqs:
                 do.call(cbind, 
                         lapply(mafs_11_list, function(m) m[ , 4]))) # take just pop freq
plot(meta.pop$lat, rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "latitude",
     ylab = "M allele frequency",
     main = "chr 11",
     ylim = c(0,1))

mafs_11 %>%
  filter(freq_M > .9 & type == "phat") %>%
  #filter(pos == 631741) %>%
  filter(pos > 1.45*10^7 & pos < 1.46*10^7) %>%
  apply(., 1, function(x)
    points(meta.pop$lat, x[meta.pop$population], col = "orange"))
mafs_11 %>%
  filter(freq_M > .9 & type == "phat") %>%
  filter(pos > 1.3*10^7 & pos < 1.31*10^7) %>%
  apply(., 1, function(x)
    points(meta.pop$lat, x[meta.pop$population], col = "blue"))
# split N and S America in plotting & plot only 1 AIM, nearest to peak:
par(mfrow = c(2, 1))
# N. America
plot(abs(meta.pop$lat), rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "Degrees latitude from the equator",
     ylab = "Allele frequency",
     main = "N. America",
     ylim = c(0,1))
mafs_11 %>%
  filter(freq_M > .9 & type == "phat") %>%
  filter(pos > 1.45*10^7 & pos < 1.46*10^7) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "N. America"]), 
           x[meta.pop$population[meta.pop$zone == "N. America"]], col = "orange"))
mafs_11 %>%
  filter(freq_M > .9 & type == "phat") %>%
  filter(pos > 1.32*10^7 & pos < 1.33*10^7) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "N. America"]), 
           x[meta.pop$population[meta.pop$zone == "N. America"]], col = "blue"))
# S. America
plot(abs(meta.pop$lat), rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "Degrees latitude from the equator",
     ylab = "Allele frequency",
     main = "S. America",
     ylim = c(0,1))
mafs_11 %>%
  filter(freq_M > .9 & type == "phat") %>%
  filter(pos > 1.45*10^7 & pos < 1.46*10^7) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "S. America"]), 
           x[meta.pop$population[meta.pop$zone == "S. America"]], col = "orange"))
mafs_11 %>%
  filter(freq_M > .9 & type == "phat") %>%
  filter(pos > 1.32*10^7 & pos < 1.33*10^7) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "S. America"]), 
           x[meta.pop$population[meta.pop$zone == "S. America"]], col = "blue"))
par(mfrow = c(1, 1))


# for chr1 find nearest A AIMs for outlier regions & plot the allele freq clines NA and SA
