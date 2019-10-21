# script to plot clines in ancestry, e.g. with latitude and environment
library(sp)
library(raster)
library(geosphere)
library(rethinking)
library(dplyr)
library(ggplot2)
library(betareg) # alternative ML fitting
source("../colors.R") # get color palette

#retrieve_data_new <- T
retrieve_data_new <- F
#res_bioclim <- 2.5
res_bioclim <- 0.5 # res 2.5 is about 4.5 km^2 at the equator; 0.5 is < 1 km^2

# get admixture data
#prefix <- "CA_AR_MX_harpur_sheppard_kohn_wallberg"

# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table(paste0("../bee_samples_listed/combined_sept19.list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
K = 3 # 3 admixing populations
n = 250 # snps thinned to 1 every nth
#prefix1 = paste0("ordered_scaffolds_", prefix, "_prunedBy", n)
#name = paste0("K", K, "_", prefix1)
name = "K3_combined_sept19_chr_prunedBy250"
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


# get climate data - bioclim1.4 with climatic variable averages 1960-1990
# Bioclimatic variables: https://www.worldclim.org/bioclim . 
# monthly data also available here: http://worldclim.org/version2
bioclim19_names <- c("AnnualMeanTemp", "MeanDiurnal", "Isothermality", "TempSeasonality",
                     "MaxTempWarmestMonth", "MinTempColdestMonth", "TempAnnualRange",
                     "MeanTempWettestQuarter", "MeanTempDriestQuarter",
                     "MeanTempWarmestQuarter", "MeanTempColdestQuarter",
                     "AnnualPrecip", "PrecipWettestMonth", "PrecipDriestMonth",
                     "PrecipSeasonality", "PrecipWettestQuarter", "PrecipDriestQuarter",
                     "PrecipWarmestQuarter", "PrecipColdestQuarter")
if (retrieve_data_new){
  # first get bee lat/long
  meta.ind <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
    left_join(ids.pops, ., by = c("Bee_ID", "population"))
  meta.ind.included <- meta.ind[meta.ind$population %in% c("Avalon_2014", "Placerita_2014", "Riverside_2014",
                                                           "Stanislaus_2014", "Stebbins_2014", "MX10") |
                                  meta.ind$group %in% c("CA_2018", "AR_2018"), ]
  
  bees.clim <- do.call(rbind,
          lapply(1:nrow(meta.ind.included), function(i) raster::extract(getData('worldclim', var = "bio", res = res_bioclim, 
                                                                                lon = meta.ind.included$long[i], 
                                                                                lat = meta.ind.included$lat[i]),
                                                         SpatialPoints(meta.ind.included[i, c("long", "lat")]))/10)) %>%
    data.frame(., stringsAsFactors = F) %>%
    #data.table::setnames(., paste0("bio", 1:19)) %>%
    data.table::setnames(., bioclim19_names) %>%
    bind_cols(meta.ind.included, .)
  
  # add distance to sao paulo, Brazil
  # Rio Claro, Sao Paulo Brazil: 22.4149° S, 47.5651° W (Google maps)
  sao_paulo <- data.frame(long = -47.5651, lat = -22.4149)
  bees.clim$km_from_sao_paulo <- apply(bees.clim[ , c("long", "lat")], 1,
                                       function(x) distm(x, sao_paulo, 
                                                         fun = distGeo))/1000
  
  # write data
  write.table(bees.clim, 
              paste0("results/BioclimVar_bees_", res_bioclim, "min.txt"), sep = "\t",
              quote = F, col.names = T, row.names = F)
} else{ # read in data already processed
  # read in admixture data
  #load("../global_ancestry/results/NGSAdmix/ACM_K3_combined_sept19_chr_prunedBy250.rData")
  # read in climate data
  bees.clim <- read.table(paste0("results/BioclimVar_bees_", res_bioclim, "min.txt"), 
                          sep = "\t",
                          header = T, stringsAsFactors = F)
}

# 'tidy' formatted data
d <- bees.clim %>%
  left_join(admix2, ., by = c("Bee_ID", "population")) %>% 
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
  mutate(precip_c = AnnualPrecip - mean(AnnualPrecip)) %>%
  mutate(from2014 = year - 2014) %>% # any evidence the hybrid zone has moved 2014 to 2018?
  mutate(from2014_c = from2014 - mean(from2014)) %>%
  mutate(year = factor(year))



# plots
d_A %>%
  ggplot(., aes(y = MeanTempColdestQuarter, x = AnnualMeanTemp, color = continent, shape = year)) +
  geom_point()
d_A %>%
  ggplot(., aes(y = MeanTempColdestQuarter, x = abs(lat), color = continent, shape = year)) +
  geom_point()
d_A %>%
  ggplot(., aes(y = AnnualMeanTemp, x = abs(lat), color = continent, shape = year)) +
  geom_point()
d_A %>%
  ggplot(., aes(x = abs(lat), y = AnnualPrecip, color = continent, shape = year)) +
  geom_point() # basically S. America is wet and N. America is dry

# plot some climate variables:
d_A %>%
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "AnnualPrecip")) %>%
  ggplot(., aes(x = abs(lat), y = value, color = continent, shape = year)) +  
                #color = geographic_location_short)) + 
  geom_point(size = 1, alpha = .75) + 
  facet_wrap(~climate_var, scales = "free_y", ncol = 1) +
  xlab("Degrees latitude from the equator") + 
  theme(legend.position="bottom") + 
  labs(color = "") +
  scale_shape_manual(values = c(17, 19)) +
  scale_color_manual(values = col_NA_SA_both) +
  theme_classic()
ggsave("plots/climate_variables_across_latitude.png", 
       height = 5, width = 5, units = "in")
# save in figures for manuscript:
ggsave("../../bee_manuscript/figures/climate_variables_across_latitude.png", 
       height = 5, width = 5, units = "in")

# A ancestry vs. latitude:
d_A %>%
  ggplot(., aes(x = abs(lat), y = alpha, color = continent)) +
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)

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
summary(m_lat_betareg)

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

# all individual variables with SA
m_cold_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_cold*cold_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_SA = 0))

m_temp_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_temp*temp_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_SA = 0))

m_precip_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_precip = 0, b_SA = 0))

m_dist_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_SA = 0))

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

m_lat_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_precip*precip_c + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_precip ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_precip = 0, b_lat = 0, b_SA = 0))

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

m_lat_cold_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_cold*cold_c + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_lat = 0, b_SA = 0))

m_lat_dist_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_SA = 0))

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

m_lat_temp_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_SA = 0))

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

m_lat_temp_cold_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_cold*cold_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_cold = 0, b_SA = 0))
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

m_lat_temp_cold_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_cold*cold_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_cold = 0, b_precip = 0, b_SA = 0))


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

m_lat_temp_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_precip = 0, b_SA = 0))


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

m_lat_cold_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu +  b_lat*abs_lat_c + b_cold*cold_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_lat = 0, b_cold = 0, b_precip = 0, b_SA = 0))

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

m_lat_temp_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_temp*temp_c + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_lat = 0, b_SA = 0))


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

# add in latitude as a predictor:
m0_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_SA = 0))

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

m0_cold_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_SA = 0))

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
m0_temp_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_temp*temp_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_temp = 0, b_SA = 0))

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
m0_precip_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_precip = 0, b_SA = 0))

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

m1_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_SA = 0))




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
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)


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

m1_cold_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_SA = 0))



d_A %>%
  mutate(m1_cold = m1_cold.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1_cold")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1_cold = m1_cold.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1_cold, color = continent)) +
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)

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

m1_temp_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_temp*temp_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_temp = 0, b_SA = 0))



d_A %>%
  mutate(m1_temp = m1_temp.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1_temp")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1_temp = m1_temp.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1_temp, color = continent)) +
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)

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

m1_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_precip = 0, b_SA = 0))


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

m1_temp_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0,5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_temp = 0, b_precip = 0, b_SA = 0))



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
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_precip = 0))
precis(m1_cold_precip)
pairs(m1_cold_precip)

m1_cold_precip_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0,5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_precip = 0, b_SA = 0))

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
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0))
precis(m1_cold_temp)
pairs(m1_cold_temp)

m1_cold_temp_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_temp*temp_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_SA ~ dnorm(0,5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0, b_SA = 0))


d_A %>%
  mutate(m1_precip = m1_precip.post.mu) %>%
  tidyr::gather(., "data", "A_ancestry", c("alpha", "m1_precip")) %>%
  ggplot(., aes(x = abs(lat), y = A_ancestry, color = data)) +
  geom_point() +
  facet_wrap(~continent)
d_A %>%
  mutate(m1_precip = m1_precip.post.mu) %>%
  ggplot(., aes(x = alpha, y = m1_precip, color = continent)) +
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)

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
# with separate SA variable
m2_SA <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0, b_precip = 0, b_SA = 0))

# filling in missing combinations:
#m_cold_temp_precip_SA
m_cold_temp_precip_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_temp = 0, b_precip = 0, b_SA = 0))

# m_lat_temp_precip_SA
m_lat_temp_precip_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_lat*abs_lat_c + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_lat = 0, b_temp = 0, b_precip = 0, b_SA = 0))

# m_temp_precip_SA 
m_temp_precip_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_precip = 0, b_SA = 0))

#m_cold_precip_SA
m_cold_precip_SA <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_cold*cold_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_precip = 0, b_SA = 0))

#m0_temp_cold_precip_SA 
m0_temp_cold_precip_SA <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_temp = 0, b_precip = 0, b_SA = 0))

#m0_temp_cold_SA
m0_temp_cold_SA <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_temp*temp_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_temp = 0, b_SA = 0))
#m0_temp_precip_SA 
m0_temp_precip_SA <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_temp*temp_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_temp = 0, b_precip = 0, b_SA = 0))
#m0_cold_precip_SA
m0_cold_precip_SA <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_precip*precip_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_precip = 0, b_SA = 0))

#--without SA--:
#m_cold_temp_precip
m_cold_temp_precip <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_temp = 0, b_precip = 0))

# m_lat_temp_precip
m_lat_temp_precip <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_lat*abs_lat_c + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_lat = 0, b_temp = 0, b_precip = 0))

# m_temp_precip
m_temp_precip <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_temp = 0, b_precip = 0))

#m_cold_precip
m_cold_precip <- map(
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), 
    logit(p) <- mu + b_cold*cold_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_precip = 0))

#m0_temp_cold_precip
m0_temp_cold_precip <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_temp = 0, b_precip = 0))

#m0_temp_cold
m0_temp_cold <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_temp*temp_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_temp = 0))
#m0_temp_precip 
m0_temp_precip <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_temp*temp_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_temp = 0, b_precip = 0))
#m0_cold_precip
m0_cold_precip <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_dist*dist_c + b_cold*cold_c + b_precip*precip_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_cold = 0, b_precip = 0))

m_cold_temp <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_cold*cold_c + b_temp*temp_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_temp = 0))

m_cold_temp_SA <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_cold*cold_c + b_temp*temp_c + b_SA*S_America_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_SA ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_cold = 0, b_temp = 0, b_SA = 0))

# is there a change 2014 to 2018?

m_lat_2014 <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_time*from2014_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_time ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_time = 0, b_lat = 0))
precis(m_lat_2014)
pairs(m_lat_2014)

m_lat_2014_CA_only <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_time*from2014_c + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_time ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5)
  ),
  data = filter(d_A, !S_America), # CA only
  start = list(mu = -1, theta = 2, b_time = 0, b_lat = 0))
precis(m_lat_2014_CA_only, prob = .95)
pairs(m_lat_2014_CA_only)

m_lat_CA_only <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5)
  ),
  data = filter(d_A, !S_America), # CA only
  start = list(mu = -1, theta = 2, b_lat = 0))
precis(m_lat_CA_only, prob = .95)
pairs(m_lat_CA_only)

m_CA_only <- map( 
  alist(
    alpha ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = filter(d_A, !S_America), # CA only
  start = list(mu = -1, theta = 2))
precis(m_CA_only, prob = .95)
pairs(m_CA_only)
rethinking::compare(m_CA_only, m_lat_CA_only, m_lat_2014_CA_only)

m2_CA_only <- map( # quadratic approximation of the posterior MAP
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
  data = filter(d_A, !S_America), # CA only
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0, b_precip = 0))
precis(m2_CA_only)
coef(m2_CA_only)

m2_2014_CA_only <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c + b_time*from2014_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_time ~ dnorm(0,5)
  ),
  data = filter(d_A, !S_America), # CA only,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0, b_precip = 0, b_time = 0))
precis(m2_2014_CA_only, prob = .95)
coef(m2_2014_CA_only)

m2_2014 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_cold*cold_c + b_temp*temp_c + b_precip*precip_c + b_time*from2014_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_cold ~ dnorm(0, 5),
    b_temp ~ dnorm(0, 5),
    b_precip ~ dnorm(0, 5),
    b_time ~ dnorm(0,5)
  ),
  data = d_A,
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_cold = 0, b_temp = 0, b_precip = 0, b_time = 0))
precis(m2_2014, prob = .95)
coef(m2_2014)


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
# so excluding varSA:
all_models_with_SA <- list(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_precip, m_lat_cold, m_lat_temp_cold, m_lat_temp_cold_precip, m_lat_cold_precip, m_lat_temp, m1, m0_temp, m0_cold, m0_precip, m2, m1_cold, m1_temp, m1_precip, m1_cold_precip, m1_temp_precip, m1_cold_temp, m_temp_precip, m_cold_precip, m_cold_temp_precip, m_lat_temp_precip, m0_temp_cold_precip, m0_temp_cold, m0_temp_precip, m0_cold_precip, m_cold_temp,
                                           m_SA, m0_SA, m_cold_SA, m_temp_SA, m_precip_SA, m_lat_SA, m_lat_precip_SA, m_lat_cold_SA, m_lat_temp_cold_SA, m_lat_temp_cold_precip_SA, m_lat_cold_precip_SA, m_lat_temp_SA, m1_SA, m0_temp_SA, m0_cold_SA, m0_precip_SA, m2_SA, m1_cold_SA, m1_temp_SA, m1_precip_SA, m1_cold_precip_SA, m1_temp_precip_SA, m1_cold_temp_SA, m_temp_precip_SA, m_cold_precip_SA, m_cold_temp_precip_SA, m_lat_temp_precip_SA, m0_temp_cold_precip_SA, m0_temp_cold_SA, m0_temp_precip_SA, m0_cold_precip_SA, m_cold_temp_SA)
all_model_names_with_SA <- c("m", "m0", "m_cold", "m_temp", "m_precip", "m_lat", "m_lat_precip", "m_lat_cold", "m_lat_temp_cold", "m_lat_temp_cold_precip", "m_lat_cold_precip", "m_lat_temp", "m1", "m0_temp", "m0_cold", "m0_precip", "m2", "m1_cold", "m1_temp", "m1_precip", "m1_cold_precip", "m1_temp_precip", "m1_cold_temp", "m_temp_precip", "m_cold_precip", "m_cold_temp_precip", "m_lat_temp_precip", "m0_temp_cold_precip", "m0_temp_cold", "m0_temp_precip", "m0_cold_precip", "m_cold_temp",
                     "m_SA", "m0_SA", "m_cold_SA", "m_temp_SA", "m_precip_SA", "m_lat_SA", "m_lat_precip_SA", "m_lat_cold_SA", "m_lat_temp_cold_SA", "m_lat_temp_cold_precip_SA", "m_lat_cold_precip_SA", "m_lat_temp_SA", "m1_SA", "m0_temp_SA", "m0_cold_SA", "m0_precip_SA", "m2_SA", "m1_cold_SA", "m1_temp_SA", "m1_precip_SA", "m1_cold_precip_SA", "m1_temp_precip_SA", "m1_cold_temp_SA", "m_temp_precip_SA", "m_cold_precip_SA", "m_cold_temp_precip_SA", "m_lat_temp_precip_SA", "m0_temp_cold_precip_SA", "m0_temp_cold_SA", "m0_temp_precip_SA", "m0_cold_precip_SA", "m_cold_temp_SA")
compare_all_with_SA <- rethinking::compare(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_precip, m_lat_cold, m_lat_temp_cold, m_lat_temp_cold_precip, m_lat_cold_precip, m_lat_temp, m1, m0_temp, m0_cold, m0_precip, m2, m1_cold, m1_temp, m1_precip, m1_cold_precip, m1_temp_precip, m1_cold_temp, m_temp_precip, m_cold_precip, m_cold_temp_precip, m_lat_temp_precip, m0_temp_cold_precip, m0_temp_cold, m0_temp_precip, m0_cold_precip, m_cold_temp,
                                   m_SA, m0_SA, m_cold_SA, m_temp_SA, m_precip_SA, m_lat_SA, m_lat_precip_SA, m_lat_cold_SA, m_lat_temp_cold_SA, m_lat_temp_cold_precip_SA, m_lat_cold_precip_SA, m_lat_temp_SA, m1_SA, m0_temp_SA, m0_cold_SA, m0_precip_SA, m2_SA, m1_cold_SA, m1_temp_SA, m1_precip_SA, m1_cold_precip_SA, m1_temp_precip_SA, m1_cold_temp_SA, m_temp_precip_SA, m_cold_precip_SA, m_cold_temp_precip_SA, m_lat_temp_precip_SA, m0_temp_cold_precip_SA, m0_temp_cold_SA, m0_temp_precip_SA, m0_cold_precip_SA, m_cold_temp_SA)
# doesn't really make sense to include continent in the mix of possible mechanistic covariantes like environment and distance from sao paoulo
all_models <- list(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_precip, m_lat_cold, m_lat_temp_cold, m_lat_temp_cold_precip, m_lat_cold_precip, m_lat_temp, m1, m0_temp, m0_cold, m0_precip, m2, m1_cold, m1_temp, m1_precip, m1_cold_precip, m1_temp_precip, m1_cold_temp, m_temp_precip, m_cold_precip, m_cold_temp_precip, m_lat_temp_precip, m0_temp_cold_precip, m0_temp_cold, m0_temp_precip, m0_cold_precip, m_cold_temp)
all_model_names <- c("m", "m0", "m_cold", "m_temp", "m_precip", "m_lat", "m_lat_precip", "m_lat_cold", "m_lat_temp_cold", "m_lat_temp_cold_precip", "m_lat_cold_precip", "m_lat_temp", "m1", "m0_temp", "m0_cold", "m0_precip", "m2", "m1_cold", "m1_temp", "m1_precip", "m1_cold_precip", "m1_temp_precip", "m1_cold_temp", "m_temp_precip", "m_cold_precip", "m_cold_temp_precip", "m_lat_temp_precip", "m0_temp_cold_precip", "m0_temp_cold", "m0_temp_precip", "m0_cold_precip", "m_cold_temp")
compare_all <- rethinking::compare(m, m0, m_cold, m_temp, m_precip, m_lat, m_lat_precip, m_lat_cold, m_lat_temp_cold, m_lat_temp_cold_precip, m_lat_cold_precip, m_lat_temp, m1, m0_temp, m0_cold, m0_precip, m2, m1_cold, m1_temp, m1_precip, m1_cold_precip, m1_temp_precip, m1_cold_temp, m_temp_precip, m_cold_precip, m_cold_temp_precip, m_lat_temp_precip, m0_temp_cold_precip, m0_temp_cold, m0_temp_precip, m0_cold_precip, m_cold_temp)


# to add (always SA version too). Then check for last missing 1 or 2 models from set of 64 (!).
# m_temp_precip, m_cold_precip, m_cold_temp_precip, m_lat_temp_precip, m0_temp_cold_precip, m0_temp_cold, m0_temp_precip, m0_cold_precip
print(compare_all)
compare_all_df <- compare_all@output %>%
  mutate(model = rownames(.)) %>%
  dplyr::select(c("model", "WAIC", "pWAIC", "dWAIC", "weight", "SE", "dSE"))

compare_all_vals <- do.call(bind_rows, lapply(all_models, function(x) x@optim$par)) %>%
                                 mutate(model = all_model_names) %>%
  left_join(compare_all_df, ., by = "model")
View(compare_all_vals)
compare_all_vals %>%
  arrange(., is.na(b_dist), is.na(b_lat), is.na(b_temp), is.na(b_cold), is.na(b_precip), is.na(b_SA)) %>%
  dplyr::select(c("model", "b_dist", "b_lat", "b_temp", "b_cold", "b_precip", "b_SA")) %>%
  View()

write.table(compare_all_vals, 
            "results/model_comparison_table.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)

png("plots/m2_pairs_plot.png", height = )
pairs(m2, main = "Posterior estimates and their correlations for parameters of the full model")

rethinking::compare(m, m_lat, m_lat_SA)
postcheck(m2)


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
plot(alpha ~ abs(lat), data = d_A, 
     col = NULL,
     #col = ifelse(d_A$continent == "S. America", rangi2, "skyblue"),
     ylim = c(0, 1), main = "African ancestry predicted by latitude", xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
# plot the mean mu for each observation as a line
#lines(abs(d_A$lat)[order(abs(d_A$lat))], m_lat.post.mu[order(abs(d_A$lat))])
shade(m_lat.post.HDPI, full_range_data$abs_lat, col = alpha(col_blind[1], .25)) # plot a shaded region for 95% HPDI

points(alpha ~ abs(lat), data = d_A, # plot points again on top
     col = ifelse(d_A$continent == "S. America", col_NA_SA_both["S. America"], col_NA_SA_both["N. America"]), add = T)
# plot MAP line from model coefficients:
curve(logistic(coef(m_lat)["mu"] + (x - mean(abs(d_A$lat)))*coef(m_lat)["b_lat"]), range(abs(d_A$lat)), n = 1000, add = T)
# plot MAP line from ML fit -- fits the same 
# (so it's not the prior, but an asymmetric posterior that makes the 
# mean posterior lower than the MAP model fit line)
#curve(logistic(coef(m_lat_betareg)["(Intercept)"] + 
#                 (x - mean(abs(d_A$lat)))*coef(m_lat_betareg)["abs_lat_c"]), 
#      range(abs(d_A$lat)), n = 1000, col = "red", add = T)

legend("topright", c("S. America", "N. America", "model prediction (MAP)", "95% confidence (HPDI)"),
        pch = c(1, 1, NA, 15), lty = c(NA, NA, 1, NA), col = c(col_NA_SA_both[1:2], "black", col_blind[1]))
dev.off()
# make figure again for manuscript folder:
png("../../bee_manuscript/figures/m_lat_model_prediction.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs(lat), data = d_A, 
     col = NULL,
     ylim = c(0, 1), main = "African ancestry predicted by latitude", xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
shade(m_lat.post.HDPI, full_range_data$abs_lat, col = alpha(col_blind[1], .25)) # plot a shaded region for 95% HPDI

points(alpha ~ abs(lat), data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", col_NA_SA_both["S. America"], col_NA_SA_both["N. America"]), add - T)
# plot MAP line from model coefficients:
curve(logistic(coef(m_lat)["mu"] + (x - mean(abs(d_A$lat)))*coef(m_lat)["b_lat"]), range(abs(d_A$lat)), n = 1000, add = T)
legend("topright", c("S. America", "N. America", "model prediction (MAP)", "95% confidence (HPDI)"),
       pch = c(1, 1, NA, 15), lty = c(NA, NA, 1, NA), col = c(col_NA_SA_both[1:2], "black", col_blind[1]))
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

# do a PCA of the climate variables
bioclim19 <- as.matrix(d_A[ , bioclim19_names])
# perform PCA on the correlation matrix (b/c bioclim variables are on diff scales)
pca_bioclim <- prcomp(bioclim19, scale = T)
# sdev is sqrt of eigenvalues; to get % variance explained by each PC we need it's eigenvalue/total eigenvalues
pca_bioclim_perc_variance <- pca_bioclim$sdev^2/sum(pca_bioclim$sdev^2) # after first 3 to 5, there's not a lot of variance left to explain
png("plots/bioclim19_pca12.png", height = 16, width = 16, units = "in", res = 300)
biplot(pca_bioclim,
       col = c("grey", "black"),
       xlabs = round(d_A$lat, 1),
       xlab = paste("PC1:", round(pca_bioclim_perc_variance[1], 2)),
       ylab = paste("PC2:", round(pca_bioclim_perc_variance[2], 2)),
       main = "PCA of Bioclim19 variables")
dev.off()
png("plots/bioclim19_pca34.png", height = 16, width = 16, units = "in", res = 300)
biplot(pca_bioclim,
       col = c("grey", "black"),
       xlabs = round(d_A$lat, 1),
       choices = c(3,4), # pc3 and pc4
       xlab = paste("PC3:", round(pca_bioclim_perc_variance[3], 2)),
       ylab = paste("PC4:", round(pca_bioclim_perc_variance[4], 2)),
       main = "PCA of Bioclim19 variables")
dev.off()
png("plots/bioclim19_pca56.png", height = 16, width = 16, units = "in", res = 300)
biplot(pca_bioclim,
       col = c("grey", "black"),
       xlabs = round(d_A$lat, 1),
       choices = c(5,6),
       xlab = paste("PC5:", round(pca_bioclim_perc_variance[5], 2)),
       ylab = paste("PC6:", round(pca_bioclim_perc_variance[6], 2)),
       main = "PCA of Bioclim19 variables")
dev.off()

plot(pca_bioclim$x[ , "PC1"], pca_bioclim$x[ , "PC2"])
plotArrows(pca_bioclim$rotation[ , "PC1"], pca_bioclim$rotation[ , "PC2"], col = "blue")
points(pca_bioclim$rotation, plot(PC1, PC2), col = "blue")

#eigenvectors are the columns of pca_bioclim$rotation
pca_bioclim$rotation
d_A %>% 
  bind_cols(., data.frame(pca_bioclim$x, stringsAsFactors = F)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = abs(lat), shape = continent)) +
  geom_point()
d_A %>% 
  bind_cols(., data.frame(pca_bioclim$x, stringsAsFactors = F)) %>%
  ggplot(., aes(x = PC3, y = PC4, color = abs(lat), shape = continent)) +
  geom_point()
d_A %>% 
  bind_cols(., data.frame(pca_bioclim$x, stringsAsFactors = F)) %>%
  ggplot(., aes(x = PC5, y = PC6, color = abs(lat), shape = continent)) +
  geom_point()

# linear model with first 6 PCs:
m_pc6_dist_lat <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_dist*dist_c + b_lat*abs_lat_c + b_pc1*PC1 + b_pc2*PC2 + b_pc3*PC3 + b_pc4*PC4 + b_pc5*PC5 + b_pc6*PC6,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_dist ~ dnorm(0, 5),
    b_lat ~ dnorm(0, 5),
    b_pc1 ~ dnorm(0, 5),
    b_pc2 ~ dnorm(0, 5),
    b_pc3 ~ dnorm(0, 5),
    b_pc4 ~ dnorm(0, 5),
    b_pc5 ~ dnorm(0, 5),
    b_pc6 ~ dnorm(0, 5)   
  ),
  data = bind_cols(d_A, data.frame(pca_bioclim$x, stringsAsFactors = F)),
  start = list(mu = -1, theta = 2, b_dist = 0, b_lat = 0, b_pc1 = 0, b_pc2 = 0, b_pc3 = 0, b_pc4 = 0, b_pc5 = 0, b_pc6 = 0))
precis(m_pc6_dist_lat)
coef(m_pc6_dist_lat)
pairs(m_pc6_dist_lat)

m_pc6_lat <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat*abs_lat_c + b_pc1*PC1 + b_pc2*PC2 + b_pc3*PC3 + b_pc4*PC4 + b_pc5*PC5 + b_pc6*PC6,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_pc1 ~ dnorm(0, 5),
    b_pc2 ~ dnorm(0, 5),
    b_pc3 ~ dnorm(0, 5),
    b_pc4 ~ dnorm(0, 5),
    b_pc5 ~ dnorm(0, 5),
    b_pc6 ~ dnorm(0, 5)   
  ),
  data = bind_cols(d_A, data.frame(pca_bioclim$x, stringsAsFactors = F)),
  start = list(mu = -1, theta = 2, b_lat = 0, b_pc1 = 0, b_pc2 = 0, b_pc3 = 0, b_pc4 = 0, b_pc5 = 0, b_pc6 = 0))
precis(m_pc6_lat)
coef(m_pc6_lat)
pairs(m_pc6_lat)

m_pc6 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_pc1*PC1 + b_pc2*PC2 + b_pc3*PC3 + b_pc4*PC4 + b_pc5*PC5 + b_pc6*PC6,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_pc1 ~ dnorm(0, 5),
    b_pc2 ~ dnorm(0, 5),
    b_pc3 ~ dnorm(0, 5),
    b_pc4 ~ dnorm(0, 5),
    b_pc5 ~ dnorm(0, 5),
    b_pc6 ~ dnorm(0, 5)   
  ),
  data = bind_cols(d_A, data.frame(pca_bioclim$x, stringsAsFactors = F)),
  start = list(mu = -1, theta = 2, b_pc1 = 0, b_pc2 = 0, b_pc3 = 0, b_pc4 = 0, b_pc5 = 0, b_pc6 = 0))
precis(m_pc6)
coef(m_pc6)
pairs(m_pc6)

m_pc4_lat <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat*abs_lat_c + b_pc1*PC1 + b_pc2*PC2 + b_pc3*PC3 + b_pc4*PC4,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_pc1 ~ dnorm(0, 5),
    b_pc2 ~ dnorm(0, 5),
    b_pc3 ~ dnorm(0, 5),
    b_pc4 ~ dnorm(0, 5)
  ),
  data = bind_cols(d_A, data.frame(pca_bioclim$x, stringsAsFactors = F)),
  start = list(mu = -1, theta = 2, b_lat = 0, b_pc1 = 0, b_pc2 = 0, b_pc3 = 0, b_pc4 = 0))
precis(m_pc4_lat)
pairs(m_pc4_lat)

m_pc4 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_pc1*PC1 + b_pc2*PC2 + b_pc3*PC3 + b_pc4*PC4,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_pc1 ~ dnorm(0, 5),
    b_pc2 ~ dnorm(0, 5),
    b_pc3 ~ dnorm(0, 5),
    b_pc4 ~ dnorm(0, 5)
  ),
  data = bind_cols(d_A, data.frame(pca_bioclim$x, stringsAsFactors = F)),
  start = list(mu = -1, theta = 2, b_pc1 = 0, b_pc2 = 0, b_pc3 = 0, b_pc4 = 0))
precis(m_pc4)
pairs(m_pc4)

m_pc1 <- map( # quadratic approximation of the posterior MAP
  alist(
    alpha ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_pc1*PC1,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_pc1 ~ dnorm(0, 5)
  ),
  data = bind_cols(d_A, data.frame(pca_bioclim$x, stringsAsFactors = F)),
  start = list(mu = -1, theta = 2, b_pc1 = 0))
precis(m_pc1)
pairs(m_pc1)

# using pc's instead of selected environmental variables doesn't improve the model fit.
rethinking::compare(m, m_lat, m1, m2, m_pc6_lat, m_pc6_dist_lat, m_pc4_lat, m_pc6, m_pc4, m_pc1)

############-----logistic clines w/ nls()------#########
# define logistic curve
logistic3 <- function(x, mu, b){
  1/(1 + exp(-b*(x - mu)))
}

# because bees in brazil have ~ 84% A ancestry, the true
# asymptote is 84% A not 100% A for the cline, so
# I can optionally use this to rescale the inferred logistic curves
rescale = .84
rescale = 1

# fit_d just has distance brazil
# fit_lat just has latitude
start_lat <- getInitial(alpha/rescale ~ SSlogis(abs_lat, Asym, 
                                   xmid, scal), 
                       d = d_A)
fit_lat <- nls(alpha ~ rescale*logistic3(x = abs_lat, b = b, mu = mu),
              start = list(b = 1/unname(start_lat["scal"]),
                           mu = unname(start_lat["xmid"])),
              data = d_A,
              trace = F)
sum_lat <- summary(fit_lat)
info_lat <- c(converged = sum_lat$convInfo$isConv,
              mu = sum_lat$coefficients["mu", "Estimate"],
              b = sum_lat$coefficients["b", "Estimate"],
              residual_error = sum_lat$sigma,
              AIC = AIC(fit_lat),
              w = sum_lat$coefficients["b", "Estimate"])
# how well does the model fit?
cor(d_A$alpha, predict(fit_lat)) # very high correlation!
# high correlation for S. America
cor(d_A$alpha[d_A$S_America == 1], 
    sapply(d_A$abs_lat[d_A$S_America == 1], function(x)
           logistic3(x, b = info_lat["b"], 
                     mu = info_lat["mu"])))
# pure lat. model has lower correlation for N. America
cor(d_A$alpha[d_A$S_America == 0], 
    sapply(d_A$abs_lat[d_A$S_America == 0], function(x)
      logistic3(x, b = info_lat["b"], 
                mu = info_lat["mu"]))) 
#plot
with(d_A, plot(abs_lat, alpha))
curve(logistic3(x, b = info_lat["b"], mu = info_lat["mu"]), 
      from = min(d_A$abs_lat), to = max(d_A$abs_lat), add = T,
      col = "blue", lwd = 3)
lines(d_A$abs_lat, predict(fit_lat), lty=2, col="red", lwd=3)
# what width would be expect for neutral diffusion?
km_per_degree_lat = 110.567 # estimate from the equator
curve((info_lat["w"]*km_per_degree_lat)^2/(2*pi*x), from = 20, to = 90,
      xlab = "generations since admixture",
      ylab = "neutral s.d. parent-offspring dispersal")
#1957 to 2018 is approx 30 to 60 generations depending on whether you believe 1 gen every year or every other year
(2018-1970)/2 # ~24 generations if you think of time since reaching current cline center
# how do I add other predictors to an nls() model?
# start by just fitting separate clines for NA and SA based on lat.
start_lat_NA <- getInitial(alpha/rescale ~ SSlogis(abs_lat, Asym, 
                                        xmid, scal), 
                        d = d_A[d_A$S_America == 0, ])
fit_lat_NA <- nls(alpha ~ rescale*logistic3(x = abs_lat, b = b, mu = mu),
               start = list(b = 1/unname(start_lat_NA["scal"]),
                            mu = unname(start_lat_NA["xmid"])),
               data = d_A[d_A$S_America == 0, ],
               trace = F)
start_lat_SA <- getInitial(alpha/rescale ~ SSlogis(abs_lat, Asym, 
                                           xmid, scal), 
                           d = d_A[d_A$S_America == 1, ])
fit_lat_SA <- nls(alpha ~ rescale*logistic3(x = abs_lat, b = b, mu = mu),
                  start = list(b = 1/unname(start_lat_SA["scal"]),
                               mu = unname(start_lat_SA["xmid"])),
                  data = d_A[d_A$S_America == 1, ],
                  trace = F)
summary(fit_lat)
summary(fit_lat_NA)
summary(fit_lat_SA)

# plot the two clines:
#plot
with(d_A, plot(abs_lat, alpha))
curve(logistic3(x, b = info_lat["b"], mu = info_lat["mu"]), 
      from = min(d_A$abs_lat), to = max(d_A$abs_lat), add = T,
      col = "blue", lwd = 3)

#png("../../bee_manuscript/figures/m_lat_model_prediction.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by latitude", 
     xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
points(alpha ~ abs_lat, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(rescale*logistic3(x, b = summary(fit_lat_NA)$coefficients["b", "Estimate"], 
                mu = summary(fit_lat_NA)$coefficients["mu", "Estimate"]), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, add = T)
curve(rescale*logistic3(x, b = summary(fit_lat_SA)$coefficients["b", "Estimate"], 
                mu = summary(fit_lat_SA)$coefficients["mu", "Estimate"]), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, add = T)
curve(rescale*logistic3(x, b = info_lat["b"], mu = info_lat["mu"]), 
      range(d_A$abs_lat), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)

legend("topright", c("S. America", "N. America", "model prediction (MAP)", "95% confidence (HPDI)"),
       pch = c(1, 1, NA, 15), lty = c(NA, NA, 1, NA), col = c(col_NA_SA_both[1:2], "black", col_blind[1]))
#dev.off()




# distance from brazil
start_lat_dist <- getInitial(alpha ~ SSlogis(km_from_sao_paulo, Asym, 
                                        xmid, scal), 
                        d = d_A)
fit_lat <- nls(alpha ~ logistic3(x = abs_lat, b = b, mu = mu),
               start = list(b = 1/unname(start_lat["scal"]),
                            mu = unname(start_lat["xmid"])),
               data = d_A,
               trace = F)
sum_lat <- summary(fit_lat)




