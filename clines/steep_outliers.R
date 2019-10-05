library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
library(rethinking)
library(betareg)
source("../colors.R") # for color palette
source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions
source("../local_ancestry/calc_FDRs.R") # scripts to calculate false discovery rates
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
A_snp_1 <- A_snp_all[1, ] %>% # take just 1 snp
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population")
A_snp_1 <- MVNsim_bounded[8, ] %>% # take a snp with some zeros
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population")

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

# read in results from betareg_clines_ind_loci.R:
MVN_clines <- read.table("results/ind_snp_clines/MVNsim_bounded.txt",
                         header = T, sep = "\t")
MVN_zero_clines <- read.table("results/ind_snp_clines/MVNsim_zero_bounded.txt",
                              header = T, sep = "\t")
# look at all individual clines across all loci
A_clines1 <- do.call(rbind,
                     lapply(1:5, function(i) read.table(paste0("results/ind_snp_clines/A", i, ".txt"),
                        header = T, sep = "\t")))
apply(MVN_clines, 2, hist)
with(MVN_clines, plot(mu, b_lat))
summary(MVN_clines$mu)
hist(MVN_clines$mu)
summary(MVN_zero_clines$mu)
hist(MVN_zero_clines$mu)
summary(MVN_clines$b_lat)
hist(MVN_clines$b_lat)
summary(MVN_zero_clines$b_lat)
hist(MVN_zero_clines$b_lat)
summary(A_clines1$b_lat)
hist(A_clines1$b_lat)

# plot a few clines from the simulations and from the data:
plot(1, 1, xlim = c(-4, 4), ylim = c(0,1), col = NULL, main = "some clines") # make empty plot
for (i in sample(1:nrow(A_clines1), size = 10)){
  curve(logistic(A_clines1$mu[i] + 
                   x*A_clines1$b_lat[i]), 
        from = -4, to = 4, add = T, col = "blue")
}
for (i in sample(1:nrow(MVN_clines), size = 10)){
  curve(logistic(MVN_clines$mu[i] + 
                   x*MVN_clines$b_lat[i]), 
        from = -4, to = 4, add = T, col = "red")
}
for (i in sample(1:nrow(MVN_zero_clines), size = 10)){
  curve(logistic(MVN_zero_clines$mu[i] + 
                   x*MVN_zero_clines$b_lat[i]), 
        from = -4, to = 4, add = T, col = "orange")
}
# plot the most extreme cline observed:
curve(logistic(A_clines1$mu[which.max(A_clines1$b_lat)] + 
                   x*A_clines1$b_lat[which.max(A_clines1$b_lat)]), 
        from = -4, to = 4, add = T, col = "blue", lwd = 3)
curve(logistic(A_clines1$mu[which.min(A_clines1$b_lat)] + 
                 x*A_clines1$b_lat[which.min(A_clines1$b_lat)]), 
      from = -4, to = 4, add = T, col = "blue", lwd = 3, lty = 2)
curve(logistic(A_clines1$mu[which.max(A_clines1$mu)] + 
                 x*A_clines1$b_lat[which.max(A_clines1$mu)]), 
      from = -4, to = 4, add = T, col = "blue", lwd = 3, lty = 3)
curve(logistic(A_clines1$mu[which.min(A_clines1$mu)] + 
                 x*A_clines1$b_lat[which.min(A_clines1$mu)]), 
      from = -4, to = 4, add = T, col = "blue", lwd = 3, lty = 4)
mean(A_clines1$mu)
var(A_clines1$mu)
mean(A_clines1$b_lat)
var(A_clines1$b_lat)
# plot extremes for simulations too:
curve(logistic(MVN_clines$mu[which.max(MVN_clines$b_lat)] + 
                 x*MVN_clines$b_lat[which.max(MVN_clines$b_lat)]), 
      from = -4, to = 4, add = T, col = "red", lwd = 3)
curve(logistic(MVN_clines$mu[which.min(MVN_clines$b_lat)] + 
                 x*MVN_clines$b_lat[which.min(MVN_clines$b_lat)]), 
      from = -4, to = 4, add = T, col = "red", lwd = 3, lty = 2)
curve(logistic(MVN_clines$mu[which.max(MVN_clines$mu)] + 
                 x*MVN_clines$b_lat[which.max(MVN_clines$mu)]), 
      from = -4, to = 4, add = T, col = "red", lwd = 3, lty = 3)
curve(logistic(MVN_clines$mu[which.min(MVN_clines$mu)] + 
                 x*MVN_clines$b_lat[which.min(MVN_clines$mu)]), 
      from = -4, to = 4, add = T, col = "red", lwd = 3, lty = 4)
mean(MVN_clines$mu)
var(MVN_clines$mu)
mean(MVN_clines$b_lat)
var(MVN_clines$b_lat)
# and simulation with all + or zero covariances:
curve(logistic(MVN_zero_clines$mu[which.max(MVN_zero_clines$b_lat)] + 
                 x*MVN_zero_clines$b_lat[which.max(MVN_zero_clines$b_lat)]), 
      from = -4, to = 4, add = T, col = "orange", lwd = 3)
curve(logistic(MVN_zero_clines$mu[which.min(MVN_zero_clines$b_lat)] + 
                 x*MVN_zero_clines$b_lat[which.min(MVN_zero_clines$b_lat)]), 
      from = -4, to = 4, add = T, col = "orange", lwd = 3, lty = 2)
curve(logistic(MVN_zero_clines$mu[which.max(MVN_zero_clines$mu)] + 
                 x*MVN_zero_clines$b_lat[which.max(MVN_zero_clines$mu)]), 
      from = -4, to = 4, add = T, col = "orange", lwd = 3, lty = 3)
curve(logistic(MVN_zero_clines$mu[which.min(MVN_zero_clines$mu)] + 
                 x*MVN_zero_clines$b_lat[which.min(MVN_zero_clines$mu)]), 
      from = -4, to = 4, add = T, col = "orange", lwd = 3, lty = 4)
mean(MVN_zero_clines$mu)
var(MVN_zero_clines$mu)
mean(MVN_zero_clines$b_lat)
var(MVN_zero_clines$b_lat)

# Q. do low recombination regions have steeper clines than high recombination regions?
# do loci with steep clines have low A in N. America? or average A?
# plot: maybe the mean trend and then a few example clines. possibly main text or possibly just supplement.
# how coupled would be expect these clines to be under neutrality?
# maybe ask if steep clines are more likely found in pollen vs. nectar qtls
shallow_cline <- which.max(A_clines1$b_lat)
steep_cline <- which.min(A_clines1$b_lat)
more_A <- which.max(A_clines1$mu)
less_A <- which.min(A_clines1$mu)
A[shallow_cline, ] %>%
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population") %>%
  with(., points(abs_lat_SA_c, A))
A[steep_cline, ] %>%
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population") %>%
  with(., points(abs_lat_SA_c, A, pch = 20))
A[less_A, ] %>%
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population") %>%
  with(., points(abs_lat_SA_c, A, pch = 1, col = "blue"))
A[more_A, ] %>%
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population") %>%
  with(., points(abs_lat_SA_c, A, pch = 20, col = "blue"))

# how surprised am I by my steepest clines?
# conclusion: We see an excess of both shallow and steep clines compared to MVN simulations, but mostly an excess of shallow clines
# these clines don't look dramatic to me though, especially the steep ones. And I'm not sure what to make of the shallow ones.
# tbd once I have all the data
blat_range = seq(from = range(A_clines1$b_lat)[1], 
                 to = range(A_clines1$b_lat)[2],
                 by = .1)
test_fdr_steep_clines <- sapply(blat_range, function(x) 
  fdr_1pop_low(a = x, pop = A_clines1$b_lat, 
                sims = MVN_zero_clines$b_lat))
FDRs_steep_clines_low <- sapply(FDR_values, function(p) 
  max(blat_range[test_fdr_steep_clines<p], na.rm = T))
fdr_1pop_low <- function(a, pop, sims){# takes in an ancestry, test for low ancestry in 1 zone
  obs <- sum(pop <= a)
  null <- sum(sims <= a)/length(sims)*length(pop)
  null/obs
}
obs <- sum(A_clines1$b_lat < -0.8)
null <- sum(MVN_zero_clines$b_lat <= -0.8)/length(MVN_zero_clines$b_lat)*length(A_clines1$b_lat)
head(A_clines1)
quantile(A_clines1$b_lat, .01)
summary(meanA_CA[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01)])
summary(meanA_CA[A_clines1$b_lat > quantile(A_clines1$b_lat, .01)])
A_clines1$center = -1*A_clines1$mu/A_clines1$b_lat + mean(meta.AR.order.by.lat$lat)
MVN_zero_clines$center = -1*MVN_zero_clines$mu/MVN_zero_clines$b_lat + mean(meta.AR.order.by.lat$lat)
summary(A_clines1$center)
quantile(A_clines1$center, probs = c(.01, .025, .05, .1, .5, .9, .95, .975, .99))
summary(MVN_zero_clines$center)
quantile(MVN_zero_clines$center, probs = c(.01, .025, .05, .1, .5, .9, .95, .975, .99))
table(cbind(sites_r, A_clines1)[A_clines1$center > quantile(A_clines1$center, .01), "r_bin5"])

# steep clines are possibly enriched in regions of the genome with low r:
table(sites_r$r_bin5[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01)])/table(sites_r$r_bin5)
sapply(1:16, function(i) sum(sites_r$chr_n[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01)] == i)/sum(sites_r$chr_n == i))
#table(sites_r$chr[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01)])/table(sites_r$chr)
# maybe chr6 has an excess. looks like low recombining region Group11 may be driving the effect
cbind(sites_r, A_clines1)[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01), ] %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = b_lat, color = r_bin5_factor)) +
  facet_wrap(~chr)
cbind(sites_r, A_clines1, AR = meanA_AR)[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01), ] %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = AR, color = r_bin5_factor)) +
  facet_wrap(~chr)
# Group11 possibly has some large inversion or something, not the largest peak, but the widest.
# also not colocalized with the other region of high M
cbind(sites_r, A_clines1, AR = meanA_AR)[A_clines1$b_lat <= quantile(A_clines1$b_lat, .01), ] %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = AR, color = r_bin5_factor)) +
  facet_wrap(~chr)

qqplot(MVN_clines$b_lat, A_clines1$b_lat)
abline(0, 1, col = "blue")
qqplot(MVN_clines$mu, A_clines1$mu)
abline(0, 1, col = "blue")

# taking out near zero (but neg.) covariances has barely any effect
qqplot(MVN_clines$b_lat, MVN_zero_clines$b_lat)
abline(0, 1, col = "blue")
qqplot(MVN_clines$mu, MVN_zero_clines$mu)
abline(0, 1, col = "blue")