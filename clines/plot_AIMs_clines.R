library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
  
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
                         lapply(mafs_11_list, function(m) m[ , 4]))) %>% # take just pop freq
  rename(chr = scaffold) %>%
  mutate(pos = as.numeric(pos), start = pos - 1, end = pos)
  
AR_11 <- mafs_11 %>%
  filter(type == "phat") %>%
  dplyr::select(meta.pop$population[meta.pop$zone == "S. America"]) %>%
  . %*% meta.pop$n_bees[meta.pop$zone == "S. America"]
  sapply(1:nrow(.), function(i)
  sum(.[i, meta.pop$population[meta.pop$zone == "S. America"]] * 
        meta.pop$n_bees[meta.pop$zone == "S. America"])/
    sum(meta.pop$n_bees[meta.pop$zone == "S. America"]))
overlap_aims_low_M_outliers_11 <- bedr(
  engine = "bedtools", 
  input = list(a = mafs_11[ , c("chr", "start", "end")] %>%
                 mutate(chr = as.character(chr),
                        start = as.integer(start),
                        end = as.integer(end)),
               b = read.table("../functional_analysis/results/outlier_regions/low_AR.bed", 
                              header = T) %>%
                 mutate(chr = as.character(chr),
                        min_FDR = as.numeric(min_FDR)) %>%
                 dplyr::select(chr, start, end, min_FDR)), 
  method = "map", #intersect
  params = "-sorted -g ../data/honeybee_genome/chr.lengths -c 4 -o min",
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR")) %>%
  filter(min_FDR != ".") %>%
  filter(min_FDR == "0.01") %>%
  mutate(pos = end,
         color = rainbow(nrow(.)),
         AIM = paste(chr, pos, sep = ":"))
dim(overlap_aims_low_M_outliers_11)


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

png("plots/AIMs_clines_chr11_outlier_ACM.png",
    height = 12, width = 12, units = "in", res = 600)
# plot AIMs:
par(mfrow = c(2, 1))
# N. America
plot(abs(meta.pop$lat), rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "Degrees latitude from the equator",
     ylab = "Allele frequency",
     main = "N. America",
     ylim = c(0,1))
mafs_11 %>%
  filter(pos %in% overlap_aims_low_M_outliers_11$pos) %>%
  filter(freq_M > .9 & type == "phat") %>%
  mutate(color = rainbow(nrow(.))) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "N. America"]), 
           x[meta.pop$population[meta.pop$zone == "N. America"]], 
           col = scales::alpha(x["color"], 0.5), pch = 20))
points(abs(meta.pop$lat[meta.pop$zone == "N. America"]), # add genomewide mean 
       apply(M, 2, mean)[meta.pop$population[meta.pop$zone == "N. America"]], col = "black", pch = 20)

# S. America
plot(abs(meta.pop$lat), rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "Degrees latitude from the equator",
     ylab = "Allele frequency",
     main = "S. America",
     ylim = c(0,1))
mafs_11 %>%
  filter(pos %in% overlap_aims_low_M_outliers_11$pos) %>%
  filter(freq_M > .9 & type == "phat") %>%
  mutate(color = rainbow(nrow(.))) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "S. America"]), 
           x[meta.pop$population[meta.pop$zone == "S. America"]], 
           col = scales::alpha(x["color"], 0.5), pch = 20))
points(abs(meta.pop$lat[meta.pop$zone == "S. America"]), 
       apply(M, 2, mean)[meta.pop$population[meta.pop$zone == "S. America"]], col = "black", pch = 20)

par(mfrow = c(1, 1))
dev.off()








# for chr1 find nearest A AIMs for outlier regions & plot the allele freq clines NA and SA
sites_1 <- read.table("results/AIMs/A/Group1.ACM.freqs", header = T)
# then for chr1 find nearest A AIMs for outlier region & plot the allele freq clines NA and SA
mafs_1_list <- lapply(unique(meta.pop$population), function(pop) # for all sites, get pop freq
  left_join(sites_1[ , c("scaffold", "pos", "major", "minor")], 
            read.table(paste0("results/AIMs/A/Group1/", pop, ".mafs.gz"), header = T),
            by = c("scaffold"="chromo", "pos"="position", "major", "minor")) %>%
    dplyr::select(scaffold, pos, phat, nInd) %>%
    pivot_longer(data = ., cols = c("phat", "nInd"), names_to = "type", values_to = pop))
mafs_1 <- cbind(left_join(mafs_1_list[[1]][ , c("scaffold", "pos", "type")],
                           sites_1,
                           by = c("scaffold", "pos")),
                 # add in metadata for each site, then grab pop freqs:
                 do.call(cbind, 
                         lapply(mafs_1_list, function(m) m[ , 4]))) %>% # take just pop freq
  filter(type == "phat") %>%
  mutate(chr = as.character(scaffold)) %>%
  mutate(start = pos - 1, end = pos)
# find overlap with AIMs and shared outlier regions chr1:
low.AR.intersect <- 
  #  bedr( 
  #  engine = "bedtools", 
  #  input = list(a = low.AR[ , c("chr", "start", "end")],
  #               b = low.shared.outliers[ , c("chr", "start", "end")]), 
  #  method = "intersect", 
  #  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths",
  #  check.chr = F
overlap_aims_shared_outliers <- bedr(
  engine = "bedtools", 
  input = list(a = mafs_1[ , c("chr", "start", "end")] %>%
                 mutate(start = as.integer(start)),
               b = read.table("../functional_analysis/results/outlier_regions/high_shared2.bed", 
                              header = T) %>%
                 mutate(chr = as.character(chr),
                        min_FDR = as.numeric(min_FDR)) %>%
                 dplyr::select(chr, start, end, min_FDR)), 
  method = "map", #intersect
  params = "-sorted -g ../data/honeybee_genome/chr.lengths -c 4 -o min",
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR")) %>%
  filter(min_FDR != ".") %>%
  mutate(pos = end,
         color = rainbow(nrow(.)),
         AIM = paste(overlap_aims_shared_outliers$chr, overlap_aims_shared_outliers$pos, sep = ":"))

png("plots/AIMs_clines_chr1_outlier_ACM.png",
    height = 12, width = 12, units = "in", res = 600)
# plot AIMs:
par(mfrow = c(2, 1))
# N. America
plot(abs(meta.pop$lat), rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "Degrees latitude from the equator",
     ylab = "Allele frequency",
     main = "N. America",
     ylim = c(0,1))
mafs_1 %>%
  filter(pos %in% overlap_aims_shared_outliers$pos) %>%
  filter(freq_A > .9 & type == "phat") %>%
  mutate(color = rainbow(nrow(.))) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "N. America"]), 
           x[meta.pop$population[meta.pop$zone == "N. America"]], 
           col = scales::alpha(x["color"], 0.5), pch = 20))
points(abs(meta.pop$lat[meta.pop$zone == "N. America"]), # add genomewide mean 
       apply(A, 2, mean)[meta.pop$population[meta.pop$zone == "N. America"]], col = "black", pch = 20)

# S. America
plot(abs(meta.pop$lat), rep(0.5, length(meta.pop$lat)), col = NULL,
     xlab = "Degrees latitude from the equator",
     ylab = "Allele frequency",
     main = "S. America",
     ylim = c(0,1))
mafs_1 %>%
  filter(pos %in% overlap_aims_shared_outliers$pos) %>%
  filter(freq_A > .9 & type == "phat") %>%
  mutate(color = rainbow(nrow(.))) %>%
  apply(., 1, function(x)
    points(abs(meta.pop$lat[meta.pop$zone == "S. America"]), 
           x[meta.pop$population[meta.pop$zone == "S. America"]], 
           col = scales::alpha(x["color"], 0.5), pch = 20))
points(abs(meta.pop$lat[meta.pop$zone == "S. America"]), 
       apply(A, 2, mean)[meta.pop$population[meta.pop$zone == "S. America"]], col = "black", pch = 20)

par(mfrow = c(1, 1))
dev.off()




head(mafs_1)
target_1 = which.min(abs(top_pos_shared_high - mafs_1$pos))
mafs_1[(target_1-6):target_1+6,] %>%
  filter(type == "phat")
mafs_1[target_1,]
mafs_1[target_1+2,]
mafs_1[target_1+4,]
top_pos_shared_high
mafs_1$pos[target_1]
mafs_1$pos[target_1+2]
top_pos_shared_high - mafs_1$pos[target_1]
top_pos_shared_high - mafs_1$pos[target_1+2]



