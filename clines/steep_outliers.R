library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
library(rethinking)
library(betareg)
source("../colors.R") # for color palette
source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions

# script attempts to identify outlier loci with steep clines across latitude

# I need to load the local ancestry data from plotLocalAncestryTracts.R:
load("../local_ancestry/results/pops_by_lat.RData") # contains objects A sites pops_by_lat meta.pop and meta.AR.order.by.lat 

A_loci <- cbind(sites, A) %>%
  mutate(snp_n = 1:nrow(sites)) # give each snp an integer identifier


set.seed(101)
random_snps <- data.frame(snp_n = sample(1:max(A_loci$snp_n), size = 100, replace = F),
                          snp_index = 1:100)
A_loci_test <- A_loci %>% # take a subset of 100 random loci
  filter(snp_n %in% random_snps$snp_n) %>%
  tidyr::gather(., "population", "A", colnames(A)) %>%
  left_join(., meta.pop, by = "population") %>%
  filter(zone == "S. America") %>% # only take Argentina for now, check steep outliers against CA data
  mutate(abs_lat_c = abs(lat) - mean(abs(lat))) %>%
  left_join(., random_snps, by = "snp_n")
dim(A_loci_test)

# with just a couple loci, this takes A LONG time to run! And >131 divergent transitions. Didn't mix well.
# besides, I can't load all my SNPs into memory to format the data
# I think I need to try a new strategy, namely fitting with MAP approximation
# and running 1 locus at a time
m_snp <- map2stan( # mcmc sampling
  alist(
    A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b. 
    # 'theta' is a precision parameter (lower = less precise = higher variance)
    logit(p) <- mu + b_lat[snp_index]*abs_lat_c + b_snp[snp_index],
    mu ~ dnorm(0, 5),
    theta ~ dcauchy(0, 2),
    b_lat[snp_index] ~ dnorm(b_lat_mean, b_lat_theta),
    b_lat_mean ~ dnorm(0, 5),
    b_lat_theta ~ dcauchy(0, 2),
    b_snp[snp_index] ~ dnorm(0, b_snp_theta),
    b_snp_theta ~ dcauchy(0, 2)
  ),
  data = A_loci_test,
  iter = 2000, warmup = 1000, chains = 1)
precis(m_snp)
pairs(m_snp, subset = c("mu", "theta", "b_lat_mean", "b_lat_theta", "b_snp_theta", "b_lat[1]", "b_lat[2]",
                        "b_snp[1]", "b_snp[2]")) # too big -- how do I view just some? #m_snp@stanfit@sim$samples[[1:5]]
post_m_snp <- extract.samples(m_snp)
str(post_m_snp)
with(post_m_snp[[c("mu", "theta")]], plot(mu ~ theta)) # didn't work
with(post_m_snp[ , c("mu", "theta")], plot(mu ~ theta)) # didn't work
plot(m_snp)

# try again: fit each snp independently and save b_lat and b_snp (slope and intercept) ? and theta
# really I want the best fitting logistic, that minimizes theta.. I'll look up options in R
A_snp_all <- A[ , meta.AR.order.by.lat$population] # take only SA
system.time(A_snp_1 <- A_snp_all[1, ] %>% # take just 1 snp
  tidyr::gather(., "population", "A") %>%
  left_join(., meta.AR.order.by.lat, by = "population"))


m_snp1 <- map( # quadratic approximation of the posterior MAP
  alist(
    A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b. 
    # 'theta' is a precision parameter (lower = less precise = higher variance)
    logit(p) <- mu + b_lat*abs_lat_SA_c, # mu is now the mean for this particular snp
    mu ~ dnorm(0, 5),
    theta ~ dcauchy(0, 2),
    b_lat ~ dnorm(0, 5)
  ),
  data = A_snp_1,
  start = list(mu = 0, theta = 1, b_lat = 0))
precis(m_snp1)
with(A_snp_1, plot(abs_lat_SA_c, A, main = "cline at snp 1 in Argentina"))
curve(logistic(m_snp1@coef["mu"] + x*m_snp1@coef["b_lat"]), from = -4, to = 4, add = T, col = "blue")
# now fit with betareg R package. Very similar output (priors have weak effect) and 10x faster.
m_snp1_betareg <- betareg(A ~ abs_lat_SA_c, 
                          data = A_snp_1, 
                          link = "logit")
m_snp1_betareg$coefficients
curve(logistic(m_snp1_betareg$coefficients$mean["(Intercept)"] + 
                 x*m_snp1_betareg$coefficients$mean["abs_lat_SA_c"]), 
      from = -4, to = 4, add = T, col = "orange")