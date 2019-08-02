# script to plot clines in ancestry, e.g. with latitude and environment
library(sp)
library(raster)
library(geosphere)
library(rethinking)
library(dplyr)
library(ggplot2)
library(betareg) # alternative ML fitting

retrieve_data_new <- F

if (retrieve_data_new){
  # Bioclimatic variables: https://www.worldclim.org/bioclim . 
  # monthly data also available here: http://worldclim.org/version2
  # get climate data - bioclim2
  climate <- getData("worldclim", var = "bio", res = 2.5) # res 2.5 is about 4.5 km^2 at the equator
  
  # get bee lat/long
  meta.ind <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
    left_join(bees, ., by = c("Bee_ID", "population"))
  meta.ind.included <- meta.ind[meta.ind$population %in% c("Avalon_2014", "Placerita_2014", "Riverside_2014",
                                                           "Stanislaus_2014", "Stebbins_2014", "MX10") |
                                  meta.ind$group %in% c("CA_2018", "AR_2018"), ]
  
  bees.clim <- data.frame(extract(climate, SpatialPoints(data.frame(meta.ind.included[ , c("long", "lat")], stringsAsFactors = F)))/10, 
                          stringsAsFactors = F) %>%
    data.table::setnames(., c("AnnualMeanTemp", "MeanDiurnal", "Isothermality", "TempSeasonality",
                              "MaxTempWarmestMonth", "MinTempColdestMonth", "TempAnnualRange",
                              "MeanTempWettestQuarter", "MeanTempDriestQuarter",
                              "MeanTempWarmestQuarter", "MeanTempColdestQuarter",
                              "AnnualPrecip", "PrecipWettestMonth", "PrecipDriestMonth",
                              "PrecipSeasonality", "PrecipWettestQuarter", "PrecipDriestQuarter",
                              "PrecipWarmestQuarter", "PrecipColdestQuarter")) %>%
    bind_cols(meta.ind.included, .)
  
  # add distance to sao paulo, Brazil
  # Rio Claro, Sao Paulo Brazil: 22.4149° S, 47.5651° W (Google maps)
  sao_paulo <- data.frame(long = -47.5651, lat = -22.4149)
  bees.clim$km_from_sao_paulo <- apply(bees.clim[ , c("long", "lat")], 1,
                                       function(x) distm(x, sao_paulo, 
                                                         fun = distGeo))/1000
  distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
  ggplot(bees.clim, aes(x = km_from_sao_paulo, y = lat)) +
    geom_point()
  
  # write data
  write.table(bees.clim, "results/BioclimVar_bees_2.5min.txt", sep = "\t",
              quote = F, col.names = T, row.names = F)
} else{ # read in data already processed
  # read in data
  bees.clim <- read.table("results/BioclimVar_bees_2.5min.txt", sep = "\t",
                          header = T, stringsAsFactors = F)
}

# add admixture data
prefix <- "CA_AR_MX_harpur_sheppard_kohn_wallberg"
# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table(paste0("../bee_samples_listed/", prefix, ".list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
K = 3 # 3 admixing populations
n = 250 # snps thinned to 1 every nth
prefix1 = paste0("ordered_scaffolds_", prefix, "_prunedBy", n)
name = paste0("K", K, "_", prefix1)
file = paste0("../global_ancestry/results/NGSAdmix/", name, ".qopt")
admix <- read.table(file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3)
# get meta data for all individuals included in NGSadmix analysis (plus extras)
ids.pops <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  dplyr::select(Bee_ID, population)
admix2 <- bind_cols(IDs, admix) %>%
  left_join(ids.pops, by = "Bee_ID")
# label ancestries
anc_labels <- data.frame(ancestry = colnames(admix),
                  ancestry_label = sapply(colnames(admix), 
                                          function(x) names(which.max(tapply(admix2[ , x], admix2$population, sum)))),
                  stringsAsFactors = F)

# 'tidy' formatted data
d <- bees.clim %>%
  left_join(., admix2, by = c("Bee_ID", "population")) %>% 
  tidyr::gather(., "ancestry", "alpha", colnames(admix)) %>%
  left_join(., anc_labels, by = "ancestry") %>%
  mutate(continent = ifelse(geographic_location == "Argentina", "S. America", "N. America"))


# make some new variables, standardized for better model fitting
d_A <- d %>%
  filter(geographic_location != "Mexico") %>% # filter out mexico
  # limit to A ancestry
  filter(ancestry_label == "A") %>%
  # absolute latitude
  mutate(abs_lat = abs(lat)) %>%
  # mean absolute latitude, centered
  mutate(abs_lat_c = abs(lat) - mean(abs(lat))) %>%
  # distance from sao paulo, centered, and in units of 1000km
  mutate(dist_c = (km_from_sao_paulo - mean(km_from_sao_paulo))/1000) %>%
  mutate(S_America = ifelse(continent == "S. America", 1, 0)) %>%
  mutate(S_America_c = S_America - mean(S_America)) %>% # center for better fit and reduce correlation with mu estimates
  mutate(temp_c = AnnualMeanTemp - mean(AnnualMeanTemp)) %>%
  mutate(cold_c = MeanTempColdestQuarter - mean(MeanTempColdestQuarter)) %>%
  mutate(precip_c = AnnualPrecip - mean(AnnualPrecip))


# plots
d_A %>%
  ggplot(., aes(y = MeanTempColdestQuarter, x = AnnualMeanTemp, color = geographic_location_short)) +
  geom_point()
d_A %>%
  ggplot(., aes(y = MeanTempColdestQuarter, x = abs(lat), color = geographic_location_short)) +
  geom_point()
d_A %>%
  ggplot(., aes(y = AnnualMeanTemp, x = abs(lat), color = geographic_location_short)) +
  geom_point()
d_A %>%
  ggplot(., aes(x = abs(lat), y = AnnualPrecip, color = geographic_location_short)) +
  geom_point() # basically S. America is wet and N. America is dry

# plot some climate variables:
d_A %>%
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "AnnualPrecip")) %>%
  ggplot(., aes(x = abs(lat), y = value, color = geographic_location_short)) +
  geom_point(size = .5) + 
  facet_wrap(~climate_var, scales = "free_y", ncol = 1) +
  xlab("Degrees latitude from the equator") + 
  theme(legend.position="bottom") + 
  labs(color = "")
ggsave("plots/climate_variables_across_latitude.png", 
       height = 5, width = 5, units = "in")
# save in figures for manuscript:
ggsave("../../bee_manuscript/figures/climate_variables_across_latitude.png", 
       height = 5, width = 5, units = "in")

# A ancestry vs. latitude:
d_A %>%
  ggplot(., aes(x = abs(lat), y = alpha, color = continent)) +
  geom_point()

# do linear models using latitude, temp, and distance from sao paulo brazil as variables.
m <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2))
precis(m)
pairs(m)

m_cold <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_cold*cold_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0))
precis(m_cold)
pairs(m_cold)

m_temp <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0))
precis(m_temp)
pairs(m_temp)

m_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_precip = 0))
precis(m_precip)
pairs(m_precip)

# just latitude
m_lat <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_lat = 0))
precis(m_lat)
pairs(m_lat)
# alternative ML model fit:
m_lat_betareg <- betareg(alpha ~ abs_lat_c,
                     link = "logit",
                     data = d_A) # basically gives the same result


# latitude and continent:
m_lat_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_lat = 0, b_SA = 0))
precis(m_lat_SA)
pairs(m_lat_SA)


m_lat_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_precip*precip_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_precip ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_precip = 0, b_lat = 0))
precis(m_lat_precip)
pairs(m_lat_precip)

m_lat_cold <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_cold*cold_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_lat = 0))
precis(m_lat_cold)
pairs(m_lat_cold)

# also add in different variance for S. America vs. N. America (doesn't work):
m_lat_varSA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = (1 - S_America)*theta_NA + S_America*theta_SA), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu +  b_lat*abs_lat_c + mu_SA*S_America,
    mu ~ dnorm(0, 5),
    mu_SA ~ dnorm(0, 5),
    theta_NA ~ dunif(0, 50),
    theta_SA ~ dunif(0, 50),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta_NA = 7, theta_SA = 7, b_lat = 0))
precis(m_lat_varSA)
pairs(m_lat_varSA)


m_lat_temp <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0))
precis(m_lat_temp)
pairs(m_lat_temp)

m_lat_temp_cold <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_cold*cold_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_cold = 0))
precis(m_lat_temp_cold)
pairs(m_lat_temp_cold)

m_lat_temp_cold_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_cold*cold_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_cold = 0, b_precip = 0))
precis(m_lat_temp_cold_precip)
pairs(m_lat_temp_cold_precip)

m_lat_temp_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_precip = 0))
precis(m_lat_temp_precip)
pairs(m_lat_temp_precip)

m_lat_cold_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu +  b_lat*abs_lat_c + b_cold*cold_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_lat = 0, b_cold = 0, b_precip = 0))
precis(m_lat_cold_precip)
pairs(m_lat_cold_precip)

m_lat_temp <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0))
precis(m_lat_temp)
pairs(m_lat_temp)



# add in latitude as a predictor:
m1 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0))
precis(m1)
pairs(m1) # predictions for distance and latitude are negatively correlated




# alpha ~ logit(a) # to constrain from 0 to 1 
# a ~ mu + b_dist*km_from_sao_paulo + b_lat*lat # linear model
# should I center lat or dist?
m0 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0))
precis(m0)
pairs(m0)

m0_cold <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0))
precis(m0_cold)
pairs(m0_cold)

m0_temp <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_temp*temp_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_temp = 0))
precis(m0_temp)
pairs(m0_temp)

m0_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_precip = 0))
precis(m0_precip)
pairs(m0_precip)

# add in latitude as a predictor:
m1 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0))
precis(m1)
pairs(m1) # predictions for distance and latitude are negatively correlated

# plot model predictions:
m1.post <- sim(m1, d_A, n = 100)
m1.post.mu <- apply(m1.post, 2, mean)
m1.post.HPDI <- apply(m1.post, 2, HPDI, prob=0.95)
d_A %>%
  mutate(m1 = m1.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1 = m1.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1, color = continent)) +
  geom_point()


# compare fit using latitude with fit using other environmental predictors (e.g. temp)
# add in latitude as a predictor:
m1_cold <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0))
precis(m1_cold)
pairs(m1_cold)
m1_cold.post <- sim(m1_cold, d_A, n = 100)
m1_cold.post.mu <- apply(m1_cold.post, 2, mean)
d_A %>%
  mutate(m1_cold = m1_cold.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1_cold")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1_cold = m1_cold.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1_cold, color = continent)) +
  geom_point()

# mean temperature:
m1_temp <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_temp*temp_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_temp = 0))
precis(m1_temp)
pairs(m1_temp)
m1_temp.post <- sim(m1_temp, d_A, n = 100)
m1_temp.post.mu <- apply(m1_temp.post, 2, mean)
d_A %>%
  mutate(m1_temp = m1_temp.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1_temp")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1_temp = m1_temp.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1_temp, color = continent)) +
  geom_point()

# mean precipitation
m1_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_precip = 0))
precis(m1_precip)
pairs(m1_precip)
m1_precip.post <- sim(m1_precip, d_A, n = 100)
m1_precip.post.mu <- apply(m1_precip.post, 2, mean)

m1_temp_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_temp = 0, b_precip = 0))
precis(m1_temp_precip)
pairs(m1_temp_precip)

m1_cold_precip <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_col ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_precip = 0))
precis(m1_cold_precip)
pairs(m1_cold_precip)

m1_cold_temp <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_temp*temp_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_col ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0))
precis(m1_cold_temp)
pairs(m1_cold_temp)


d_A %>%
  mutate(m1_precip = m1_precip.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1_precip")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1_precip = m1_precip.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1_precip, color = continent)) +
  geom_point()

# plot instead with shading
# plot raw data
# fading out points to make line and interval more visible
plot(alpha ~ abs(lat), data = d_A, col = col.alpha(rangi2, 0.5) )
# plot the MAP line, aka the mean mu for each weight
lines(abs(d_A$lat)[order(abs(d_A$lat))], m1.postmu[order(abs(d_A$lat))])
# plot a shaded region for 89% HPDI
shade(m1.post.HPDI[ , order(abs(d_A$lat))], abs(d_A$lat)[order(abs(d_A$lat))])


m2 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0, b_precip = 0))
precis(m2)
coef(m2)

# plot some predictions of the model against my data

# compare some models:
rethinking::compare(m0, m1, m1_cold, m1_temp, m1_precip, m2)
# temp and cold are correlated w/ lat and so having both doesn't add much to prediction
rethinking::compare(m0, m1, m1_cold, m1_temp, m1_precip) 
# latitude is the best individual predictor when starting w/ distance; 
# though precipitation is collinear w/ distance
rethinking::compare(m0, m1, m0_cold, m0_temp, m0_precip) # precipitation is correlated w/ distance
# what if I don't include distance from sao paulo?
rethinking::compare(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_precip, m1, m0_temp, m0_cold, m0_precip, m2)
# the bets model is the most complicated (m2) with all predictors,
# but from the simpler models, latitude is the single best predictor,
# and latitude + cold (mean temp in coldest quarter) is the best double-predictor (despite correlations)
rethinking::compare(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_precip, m_lat_cold, m_lat_temp, m1, m0_temp, m0_cold, m0_precip)
rethinking::compare(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_SA, m_lat_precip, m_lat_cold, m_lat_temp_cold, m_lat_temp_cold_precip, m_lat_cold_precip, m_lat_temp, m1, m0_temp, m0_cold, m0_precip, m2, m1_cold, m1_temp, m1_precip, m1_cold_precip, m1_temp_precip, m1_cold_temp, m_lat_varSA)
# I'm not sure it really makes sense to give the two continents different variances,
# but what it's capturing is more divergence from the logit in SA, not variance in alpha for a given latitude
# which is what I (naively) thought it might capture

postcheck(m2)
rethinking::compare(m, m_lat, m_lat_SA)

# alternatively, I could cross-validate predictors if I use one hybrid zone 
# to predict the other (?) using either latitude or temp, dist etc.



# e.g. beta distribution for mean of all predictors:
curve(rethinking::dbeta2(x, prob = logistic(coef(m1)["mu"]), 
                         theta = coef(m1)["theta"]), from = 0, to = 1,
       ylab = "Beta Density", xlab = "A ancestry proportion", ylim=c(0, 3), lwd=2)

# plot raw data and a line for model prediction based on latitude. Color by cold. Shape = S. America or N. America.
# m_lat best single-predictor model:
full_range_data = list(abs_lat_c = seq(from = range(d_A$abs_lat_c)[1], 
                                       to = range(d_A$abs_lat_c)[2], 
                                       length.out = 20),
                       abs_lat = seq(from = range(d_A$abs_lat_c)[1], 
                                     to = range(d_A$abs_lat_c)[2], 
                                     length.out = 20) + mean(abs(d_A$lat)))
m_lat.post <- sim(m_lat, full_range_data, n = 10000)
m_lat.post.mu <- apply(m_lat.post, 2, mean)
m_lat.post.mode <- apply(m_lat.post, 2, mode)
m_lat.post.median <- apply(m_lat.post, 2, median)
m_lat.post.HDPI <- apply(m_lat.post, 2, HPDI, prob=0.95)

# plot model prediction for m_lat:
png("plots/m_lat_model_prediction.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs(lat), data = d_A, col = ifelse(d_A$continent == "S. America", rangi2, "skyblue"),
     ylim = c(0, 1), main = "African ancestry predicted by latitude", xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
# plot the mean mu for each observation as a line
#lines(abs(d_A$lat)[order(abs(d_A$lat))], m_lat.post.mu[order(abs(d_A$lat))])

# plot MAP line from model coefficients:
curve(logistic(coef(m_lat)["mu"] + (x - mean(abs(d_A$lat)))*coef(m_lat)["b_lat"]), range(abs(d_A$lat)), n = 1000, add = T)
# plot MAP line from ML fit -- fits the same 
# (so it's not the prior, but an asymmetric posterior that makes the 
# mean posterior lower than the MAP model fit line)
#curve(logistic(coef(m_lat_betareg)["(Intercept)"] + 
#                 (x - mean(abs(d_A$lat)))*coef(m_lat_betareg)["abs_lat_c"]), 
#      range(abs(d_A$lat)), n = 1000, col = "red", add = T)

# plot a shaded region for 95% HPDI
shade(m_lat.post.HDPI, full_range_data$abs_lat)
legend("topright", c("S. America", "N. America", "model prediction (MAP)", "95% confidence (HPDI)"),
        pch = c(1, 1, NA, 15), lty = c(NA, NA, 1, NA), col = c(rangi2, "skyblue", "black", "grey"))
dev.off()
# make figure again for manuscript folder:
png("../../bee_manuscript/figures/m_lat_model_prediction.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs(lat), data = d_A, col = ifelse(d_A$continent == "S. America", rangi2, "skyblue"),
     ylim = c(0, 1), main = "African ancestry predicted by latitude", xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
# plot MAP line from model coefficients:
curve(logistic(coef(m_lat)["mu"] + (x - mean(abs(d_A$lat)))*coef(m_lat)["b_lat"]), range(abs(d_A$lat)), n = 1000, add = T)
# plot a shaded region for 95% HPDI
shade(m_lat.post.HDPI, full_range_data$abs_lat)
legend("topright", c("S. America", "N. America", "model prediction (MAP)", "95% confidence (HPDI)"),
       pch = c(1, 1, NA, 15), lty = c(NA, NA, 1, NA), col = c(rangi2, "skyblue", "black", "grey"))
dev.off()




d_A %>%
  ggplot(., aes(x = abs(lat), y = alpha, color = MeanTempColdestQuarter)) +
  geom_point() +
  facet_wrap(~continent) + # I'm not sure how to plot this curve
  stat_function(aes(x = abs(lat), fun = logistic(coef(m_lat)["mu"] + coef(m_lat)["b_lat"]*abs(lat))))
# it would be nicer to just plot minimum temperature in coldest quarter on a map


d_A %>%
  mutate(m1 = m1.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1 = m1.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1, color = continent)) +
  geom_point()




