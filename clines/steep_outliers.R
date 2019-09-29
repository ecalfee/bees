library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
library(rethinking)
source("../colors.R") # for color palette
source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions

# script attempts to identify outlier loci with steep clines across latitude

# I need to load the local ancestry data from plotLocalAncestryTracts.R:
# objects: A, meta.ind, meta.pop, sites

A_loci <- cbind(sites, A) %>%
  mutate(snp_n = 1:nrow(sites)) # give each snp an integer identifier



A_loci_test <- A_loci %>%
  filter(snp_n <= 30) %>%
  tidyr::gather(., "population", "A", colnames(A)) %>%
  left_join(., meta.pop, by = "population") %>%
  mutate(abs_lat_c = abs(lat) - mean(abs(lat))) %>%
  filter(zone == "S. America") # only take Argentina for now, check steep outliers against CA data
dim(A_loci_test)

# with just a couple loci, this takes A LONG time to run! And >131 divergent transitions. Didn't mix well.
# besides, I can't load all my SNPs into memory to format the data
# I think I need to try a new strategy, namely fitting with MAP approximation
# and running 1 locus at a time
m_snp <- map2stan( # quadratic approximation of the posterior MAP
  alist(
    A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat[snp_n]*abs_lat_c + b_snp[snp_n],
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat[snp_n] ~ dnorm(b_lat_mean, b_lat_theta),
    b_lat_mean ~ dnorm(0, 5),
    b_lat_theta ~ dcauchy(0, 2),
    b_snp[snp_n] ~ dnorm(0, b_snp_theta),
    b_snp_theta ~ dcauchy(0, 2)
  ),
  data = A_loci_test,
  iter = 1000, warmup = 500, chains = 1)
precis(m_snp)
pairs(m_snp@stanfit@sim$samples[[1:5]]) # too big -- how do I view just some?
plot(m_snp)

# try again: fit each snp independently and save b_lat and b_snp (slope and intercept) ? and theta
# really I want the best fitting logistic, that minimizes theta.. I'll look up options in R

