# script to plot clines in ancestry, e.g. with latitude and environment
library(sp)
library(raster)
library(geosphere)
library(rethinking)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(purrr)
library(nlstools)
library(betareg) # alternative ML fitting
source("../colors.R") # get color palette
#install.packages("nls.multstart")
library(nls.multstart)
source("../local_ancestry/calc_FDRs.R")
load("../local_ancestry/results/pops_by_lat.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 
load("../local_ancestry/results/A.RData") # ancestry data for each pop

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
  # load and merge rasters (merging is slow!)
  # rasters for argentina:
  ar1 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -70, 
                lat = -50)
  ar2 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -50, 
                lat = -20)
  ar3 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -50, 
                lat = -50)
  ar4 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -70, 
                lat = -20)
  # rasters for ca:
  ca1 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -120, 
                lat = 30)
  ca2 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -120, 
                lat = 50)
  ca3 = getData('worldclim', var = "bio", res = res_bioclim, 
                lon = -150, 
                lat = 50)
  #sapply(list(ca1, ca2, ca3), function(x) names(x)[1])
  #plot(merge(ca1[[1]], ca2[[1]], ca3[[1]]))
  #plot(merge(ar1[[1]], ar2[[1]], ar3[[1]], ar4[[1]]))
  keep_bioclim <- c(1,6,11,12) # variables of interest
  raster_SA = merge(ar1[[keep_bioclim]], ar2[[keep_bioclim]], ar3[[keep_bioclim]], ar4[[keep_bioclim]])
  raster_NA = merge(ca1[[keep_bioclim]], ca2[[keep_bioclim]], ca3[[keep_bioclim]])
  bees.clim <- do.call(rbind,
          lapply(1:nrow(meta.ind.included), function(i) {if (meta.ind.included$group[i] == "AR_2018") {
            # avg. within 5km buffer radius. divide output by 10 for more intuitive units
            raster::extract(raster_SA, SpatialPoints(meta.ind.included[i, c("long", "lat")]), buffer = 5000, fun = mean)/10}
            else{ # north america raster
              raster::extract(raster_NA, SpatialPoints(meta.ind.included[i, c("long", "lat")]), buffer = 5000, fun = mean)/10
            }})) %>%
    data.frame(., stringsAsFactors = F) %>%
    data.table::setnames(., bioclim19_names[keep_bioclim]) %>%
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

iHot <- which.max(d_A$AnnualMeanTemp)
iCold <- which.min(d_A$AnnualMeanTemp)
bees.clim[c(iHot, iCold), ]
distm(d_A[iHot, c("long", "lat")], d_A[iCold, c("long", "lat")], 
      fun = distGeo)/1000 # distance in km between hottest and coldest site (CA mountains and desert)

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
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "MinTempColdestMonth", "AnnualPrecip")) %>%
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
ggsave("../../bee_manuscript/figures/climate_variables_across_latitude.pdf", 
       height = 5, width = 5, units = "in")

# A ancestry vs. latitude:
d_A %>%
  ggplot(., aes(x = abs(lat), y = alpha, color = continent)) +
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)

# do linear models using latitude, temp, and distance from sao paulo brazil as variables.
m <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_cold <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_temp <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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
m_lat <- rethinking::map( # quadratic approximation of the posterior MAP
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
m_lat_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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
m_cold_SA <- rethinking::map(
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

m_temp_SA <- rethinking::map(
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

m_precip_SA <- rethinking::map(
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

m_dist_SA <- rethinking::map(
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

m_lat_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_cold <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_cold_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_dist_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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
m_lat_varSA <- rethinking::map( # quadratic approximation of the posterior MAP
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


m_lat_temp <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_cold <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_cold_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_cold_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_cold_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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


m_lat_temp_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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


m_lat_cold_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_cold_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_lat_temp_SA <- rethinking::map(
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
m0 <- rethinking::map( # quadratic approximation of the posterior MAP
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
m0_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m0_cold <- rethinking::map( # quadratic approximation of the posterior MAP
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

m0_cold_SA <- rethinking::map(
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

m0_temp <- rethinking::map( # quadratic approximation of the posterior MAP
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
m0_temp_SA <- rethinking::map(
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

m0_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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
m0_precip_SA <- rethinking::map(
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
m1 <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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
m1_cold <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_cold_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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
m1_temp <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_temp_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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
m1_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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


m1_temp_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_temp_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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



m1_cold_precip <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_cold_precip_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_cold_temp <- rethinking::map( # quadratic approximation of the posterior MAP
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

m1_cold_temp_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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


m2 <- rethinking::map( # quadratic approximation of the posterior MAP
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
m2_SA <- rethinking::map( # quadratic approximation of the posterior MAP
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
m_cold_temp_precip_SA <- rethinking::map(
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
m_lat_temp_precip_SA <- rethinking::map(
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
m_temp_precip_SA <- rethinking::map(
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
m_cold_precip_SA <- rethinking::map(
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
m0_temp_cold_precip_SA <- rethinking::map( 
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
m0_temp_cold_SA <- rethinking::map( 
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
m0_temp_precip_SA <- rethinking::map( 
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
m0_cold_precip_SA <- rethinking::map( 
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
m_cold_temp_precip <- rethinking::map(
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
m_lat_temp_precip <- rethinking::map(
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
m_temp_precip <- rethinking::map(
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
m_cold_precip <- rethinking::map(
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
m0_temp_cold_precip <- rethinking::map( 
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
m0_temp_cold <- rethinking::map( 
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
m0_temp_precip <- rethinking::map( 
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
m0_cold_precip <- rethinking::map( 
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

m_cold_temp <- rethinking::map( 
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

m_cold_temp_SA <- rethinking::map( 
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

m_lat_2014 <- rethinking::map( 
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

m_lat_2014_CA_only <- rethinking::map( 
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

m_lat_CA_only <- rethinking::map( 
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

m_CA_only <- rethinking::map( 
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

m2_CA_only <- rethinking::map( # quadratic approximation of the posterior MAP
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

m2_2014_CA_only <- rethinking::map( # quadratic approximation of the posterior MAP
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

m2_2014 <- rethinking::map( # quadratic approximation of the posterior MAP
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
m_pc6_dist_lat <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_pc6_lat <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_pc6 <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_pc4_lat <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_pc4 <- rethinking::map( # quadratic approximation of the posterior MAP
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

m_pc1 <- rethinking::map( # quadratic approximation of the posterior MAP
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
#rescale1 = .84
rescale1 = 1
d_A <- mutate(d_A, rescale = rescale1)

# fit_d just has distance brazil
# fit_lat just has latitude
start_lat <- getInitial(alpha/rescale ~ SSlogis(abs_lat, Asym, 
                                   xmid, scal), 
                       d = d_A)
fit_lat0 <- nls(alpha ~ rescale*logistic3(x = abs_lat, b = b, mu = mu),
              start = list(b = 1/unname(start_lat["scal"]),
                           mu = unname(start_lat["xmid"])),
              data = d_A,
              trace = F)
sum_lat0 <- summary(fit_lat0)
info_lat0 <- c(converged = sum_lat0$convInfo$isConv,
               mu = sum_lat0$coefficients["mu", "Estimate"],
               b = sum_lat0$coefficients["b", "Estimate"],
               residual_error = sum_lat0$sigma,
               AIC = AIC(fit_lat0),
               w = sum_lat0$coefficients["b", "Estimate"]*4)



fit_lat <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                                              b = b, mu = mu),
                            start_lower = list(b = -1, mu = 25),
                            start_upper = list(b = 1, mu = 40),
                            supp_errors = 'Y',
                            iter = 250,
                            convergence_count = 100,
                            data = d_A)
glance(fit_lat)
summary(fit_lat)
coef_lat = tidy(fit_lat) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)

# how well does the model fit?
cor(d_A$alpha, predict(fit_lat)) # very high correlation!
# high correlation for S. America
cor(d_A$alpha[d_A$S_America == 1], predict(fit_lat, newdata = d_A[d_A$S_America == 1, ]))

# pure lat. model has lower correlation for N. America
cor(d_A$alpha[d_A$S_America == 0], predict(fit_lat, newdata = d_A[d_A$S_America == 0, ]))

#plot
with(d_A, plot(abs_lat, alpha))
curve(rescale1*logistic3(x, b = coef_lat$b, mu = coef_lat$mu), 
      from = min(d_A$abs_lat), to = max(d_A$abs_lat), add = T,
      col = "blue", lwd = 3)
# can do the same w/ 'predict()' but need to sort data
lines(d_A$abs_lat[order(d_A$abs_lat)], predict(fit_lat)[order(d_A$abs_lat)], lty=2, col="red", lwd=3)
# what width would be expect for neutral diffusion?
km_per_degree_lat = 110.567 # estimate from the equator
curve((coef_lat$b*4*km_per_degree_lat)^2/(2*pi*x), from = 20, to = 90,
      xlab = "generations since admixture",
      ylab = "neutral s.d. parent-offspring dispersal")
years_tot = 2018-1957
dispersal_kernel <- (coef_lat$b*4*km_per_degree_lat)^2/(2*pi*years_tot)
dispersal_kernel # variance in distance between parent and offspring
660/(1.68*sqrt(years_tot))
# a dispersal kernel variance of ~130km/year seems like a lot,  
# but may still be low for the rate of spread we saw in the Africanized honey bee invasion:
# in a gaussian dispersal kernel, this would only gain bees ~5 km a year
x = matrix(0, years_tot+1, 1000)
for (i in 2:nrow(x)){
  x[i, ] <- sapply(x[i - 1, ], function(m) rnorm(n = 1, mean = m, sd = sqrt(dispersal_kernel)))
}
max_pos = apply(x, 1, max)
year = 1:nrow(x)
plot(year, max_pos)
abline(lm(max_pos ~ year))
max(max_pos)/nrow(x) # mean distance gained per year

#1957 to 2018 is approx 30 to 60 generations depending on whether you believe 1 gen every year or every other year
(2018-1970)/2 # ~24 generations if you think of time since reaching current cline center
# how do I add other predictors to an nls() model?
# start by just fitting separate clines for NA and SA based on lat.
start_lat_NA <- getInitial(alpha/rescale ~ SSlogis(abs_lat, Asym, 
                                        xmid, scal), 
                        d = d_A[d_A$S_America == 0, ])
fit_lat_NA0 <- nls(alpha ~ rescale*logistic3(x = abs_lat, b = b, mu = mu),
               start = list(b = 1/unname(start_lat_NA["scal"]),
                            mu = unname(start_lat_NA["xmid"])),
              trace = F,
              data = d_A[d_A$S_America == 0, ])
# can alternatively fit w/ multstart
fit_lat_NA <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                              b = b, mu = mu),
                start_lower = list(b = -1, mu = 25),
                start_upper = list(b = 1, mu = 40),
                supp_errors = 'Y',
                iter = 250,
                convergence_count = 100,
                data = d_A[d_A$S_America == 0, ])
glance(fit_lat_NA)
summary(fit_lat_NA)
coef(fit_lat_NA0)

start_lat_SA <- getInitial(alpha/rescale ~ SSlogis(abs_lat, Asym, 
                                           xmid, scal), 
                           d = d_A[d_A$S_America == 1, ])
fit_lat_SA <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                                                            b = b, mu = mu),
                                          start_lower = list(b = -1, mu = 25),
                                          start_upper = list(b = 1, mu = 40),
                                          supp_errors = 'Y',
                                          iter = 250,
                                          convergence_count = 100,
                                          data = d_A[d_A$S_America == 1, ])

glance(fit_lat_SA)
summary(fit_lat_SA)
summary(fit_lat)
summary(fit_lat_NA)
summary(fit_lat_SA)
tidy(fit_lat)
info(fit_lat)

coef_lat_NA = tidy(fit_lat_NA) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)
coef_lat_SA = tidy(fit_lat_SA) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)

# ok what about fitting 1 model with a SA variable? joint fit.

fit_lat_zone_mu_and_b <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                                                                 b = b + b_SA*S_America, 
                                                                 mu = mu + mu_SA*S_America),
                                       start_lower = list(b = -1, mu = 25, b_SA = -1, mu_SA = -5),
                                       start_upper = list(b = 1, mu = 40, b_SA = 1, mu_SA = 5),
                                       supp_errors = 'Y',
                                       iter = 250,
                                       convergence_count = 100,
                                       data = d_A)
summary(fit_lat_zone_mu_and_b) # b_SA is not significant, drop:
coef_lat_zone_mu_b = tidy(fit_lat_zone_mu_and_b) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)
# only fit different centers
fit_lat_zone <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                                                        b = b, 
                                                        mu = mu + mu_SA*S_America),
                              start_lower = list(b = -1, mu = 25, mu_SA = -5),
                              start_upper = list(b = 1, mu = 40, mu_SA = 5),
                              supp_errors = 'Y',
                              iter = 250,
                              convergence_count = 100,
                              data = d_A)
summary(fit_lat_zone)
glance(fit_lat_zone_mu_and_b)
glance(fit_lat_zone)
glance(fit_lat)
coef_lat_zone = tidy(fit_lat_zone) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)


# plot the two clines:
#plot
# need to make pdf version eventually.
png("../../bee_manuscript/figures/m_lat_model_prediction.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     #main = "African ancestry predicted by latitude", 
     xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
points(alpha ~ abs_lat, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(rescale1*logistic3(x, b = coef_lat_zone_mu_b$b + coef_lat_zone_mu_b$b_SA, 
                         mu = coef_lat_zone_mu_b$mu + coef_lat_zone_mu_b$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 1, add = T)
curve(rescale1*logistic3(x, b = coef_lat_zone_mu_b$b, mu = coef_lat_zone_mu_b$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 1, add = T)
#curve(rescale1*logistic3(x, b = coef_lat$b, mu = coef_lat$mu), 
#      range(d_A$abs_lat), n = 1000,
#      col = "black", lwd = 2, add = T, lty = 2)

legend("topright", c("N. America", "S. America"), pch = c(1, 1), col = col_NA_SA_both[1:2])
dev.off()




# distance from brazil
start_dist <- getInitial(alpha ~ SSlogis(km_from_sao_paulo, 
                                             Asym, 
                                        xmid, scal), 
                        d = d_A[d_A$S_America == 1, ])
fit_dist0 <- nls(alpha/rescale ~ logistic3(x = km_from_sao_paulo, b = b, mu = mu),
               start = list(b = 1/unname(start_dist["scal"]),
                            mu = unname(start_dist["xmid"])),
               data = d_A,
               trace = F)
fit_dist <- nls_multstart(alpha ~ rescale*logistic3(x = km_from_sao_paulo, 
                                        b = b, mu = mu),
              start_lower = list(b = -5, mu = 1000),
              start_upper = list(b = 5, mu = 10000),
              supp_errors = 'Y',
              iter = 250,
              convergence_count = 100,
              data = d_A)
glance(fit_dist)
sum_dist <- summary(fit_dist)
plot(alpha ~ km_from_sao_paulo, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by dispersal distance", 
     xlab = "Distance from Sao Paulo (km)",
     ylab = "A ancestry proportion")
points(alpha ~ km_from_sao_paulo, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(rescale1*logistic3(x, b = summary(fit_dist)$coefficients["b", "Estimate"], 
                        mu = summary(fit_dist)$coefficients["mu", "Estimate"]), 
      range(d_A$km_from_sao_paulo), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)
# clearly this is a silly exercise


# also add in 2014 (vs. 2018) as something to fit
# I removed Avalon_2014 because it's such an outlier for ancestry
# and there's no 2018 sampling. I should explore if other outliers are driving this effect
# also to compare to other models I need to use same data (also excluding Avalon).
# This suggests hybrid zone advanced a little in the past 4 years.
fit_lat_zone_2014 <- nls(alpha ~ rescale*logistic3(x = abs_lat, b = b, 
                                              mu = mu + mu_SA*S_America + mu_2014*from2014),
                    start = list(b = 1/unname(start_lat["scal"]),
                                 mu = unname(start_lat["xmid"]),
                                 mu_SA = 0,
                                 mu_2014 = 0),
                    data = d_A[d_A$population != "Avalon_2014", ],
                    trace = F)
summary(fit_lat_zone_2014)

fit_winter <- nls_multstart(alpha ~ rescale*logistic3(x = MeanTempColdestQuarter, 
                                                              b = b, mu = mu),
                                    start_lower = list(b = -5, mu = min(d_A$MeanTempColdestQuarter)),
                                    start_upper = list(b = 5, mu = max(d_A$MeanTempColdestQuarter)),
                                    supp_errors = 'Y',
                                    iter = 250,
                                    convergence_count = 100,
                                    data = d_A)

# what if I fit on mean temp. instead of latitude?
start_temp <- getInitial(alpha/rescale ~ SSlogis(AnnualMeanTemp, Asym, 
                                                     xmid, scal), 
                             d = d_A)
fit_temp_no_2014 <- nls_multstart(alpha ~ rescale*logistic3(x = AnnualMeanTemp, 
                                                    b = b, mu = mu),
                          start_lower = list(b = -5, mu = min(d_A$AnnualMeanTemp)),
                          start_upper = list(b = 5, mu = max(d_A$AnnualMeanTemp)),
                          supp_errors = 'Y',
                          iter = 250,
                          convergence_count = 100,
                          data = d_A[d_A$year != 2014, ])
fit_lat_no_2014 <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                                                    b = b, mu = mu),
                          start_lower = list(b = -5, mu = min(d_A$abs_lat)),
                          start_upper = list(b = 5, mu = max(d_A$abs_lat)),
                          supp_errors = 'Y',
                          iter = 250,
                          convergence_count = 100,
                          data = d_A[d_A$year != 2014, ])
fit_winter_no_2014 <- nls_multstart(alpha ~ rescale*logistic3(x = MeanTempColdestQuarter, 
                                                           b = b, mu = mu),
                                 start_lower = list(b = -5, mu = min(d_A$MeanTempColdestQuarter)),
                                 start_upper = list(b = 5, mu = max(d_A$MeanTempColdestQuarter)),
                                 supp_errors = 'Y',
                                 iter = 250,
                                 convergence_count = 100,
                                 data = d_A[d_A$year != 2014, ])
fit_lat_zone_no_2014 <- nls_multstart(alpha ~ rescale*logistic3(x = abs_lat, 
                                                        b = b, 
                                                        mu = mu + mu_SA*S_America),
                              start_lower = list(b = -1, mu = 25, mu_SA = -5),
                              start_upper = list(b = 1, mu = 40, mu_SA = 5),
                              supp_errors = 'Y',
                              iter = 250,
                              convergence_count = 100,
                              data = d_A[d_A$year != 2014, ])

tidy(fit_temp)
glance(fit_temp)
glance(fit_winter)
glance(fit_lat)
glance(fit_lat_zone)

glance(fit_temp_no_2014)
glance(fit_winter_no_2014)
glance(fit_lat_no_2014)
glance(fit_lat_zone_no_2014)

nls(alpha ~ rescale*logistic3(x = AnnualMeanTemp, b = b, 
                                              mu = mu),
                    start = list(b = 1/unname(start_temp["scal"]),
                                 mu = unname(start_temp["xmid"])),
                    data = d_A,
                    trace = F)
summary(fit_temp)

# split up temp by zone (I'm not sure if this makes sense..)
# only fit different centers
fit_temp_zone <- nls(alpha ~ rescale*logistic3(x = AnnualMeanTemp, b = b, 
                                              mu = mu + mu_SA*S_America),
                    start = list(b = 1/unname(start_temp["scal"]),
                                 mu = unname(start_temp["xmid"]),
                                 mu_SA = 0),
                    data = d_A,
                    trace = F)
summary(fit_temp_zone)


start_winter <- getInitial(alpha/rescale ~ SSlogis(MeanTempColdestQuarter, Asym, 
                                                 xmid, scal), 
                         d = d_A)
fit_winter <- nls(alpha ~ rescale*logistic3(x = MeanTempColdestQuarter, b = b, 
                                          mu = mu),
                start = list(b = 1/unname(start_winter["scal"]),
                             mu = unname(start_winter["xmid"])),
                data = d_A,
                trace = F)
summary(fit_winter)

start_coldest <- getInitial(alpha/rescale ~ SSlogis(MinTempColdestMonth, Asym, 
                                                   xmid, scal), 
                           d = d_A)
fit_coldest <- nls(alpha ~ rescale*logistic3(x = MinTempColdestMonth, b = b, 
                                            mu = mu),
                  start = list(b = 1/unname(start_coldest["scal"]),
                               mu = unname(start_coldest["xmid"])),
                  data = d_A,
                  trace = F)
summary(fit_coldest)

# compare LL: 
all_models <- list(fit_lat, fit_lat_zone, fit_lat_zone_mu_and_b, fit_lat_zone_2014,
                   fit_temp, fit_temp_zone, fit_winter, fit_coldest)
sapply(all_models, AIC)
str(fit_lat)

# TO DO: draw predictions for 2014, CA vs. AR zone, all together.
# how diff. are these?
# also plot predictions for mean temp and min annual temp

# look up how to use tails vs. center of hybrid zone to estimate dispersal/strength of sel.
# try to fit all snps to clines -- id which ones don't fit..


# load in individual snp clines
a0 <- a
a <- read.table("results/ind_snp_nls_clines/A1.txt", # data subset
           header = T,
           stringsAsFactors = F)
str(a)
table(is.na(a$mu))
no_fit <- which(is.na(a$mu))
plot(a0$mu + mean(d_A$abs_lat[d_A$group == "AR_2018"]) ~ a$mu)
# plot a 'good' cline estimate
load(paste0("../local_ancestry/results/A.RData"))
load(paste0("../local_ancestry/results/pops_by_lat.RData"))
i = sample(no_fit, 1)
i = sample(1:nrow(a)[-no_fit], 1)
plot_a_cline(sample(1:nrow(a)[-no_fit], 1))
plot_a_cline(i = sample(no_fit, 1))
plot_a_cline <- function(i){
  plot(meta.AR.order.by.lat$lat, 
     A[i, meta.AR.order.by.lat$population],
     ylim = c(0, 1),
     xlim = range(d_A$lat[d_A$group == "AR_2018"]),
     ylab = "A ancestry - logistic",
     xlab = "latitude",
     main = paste("snp", i))
curve(logistic3(mu = a$mu[i], b = a$b[i], x), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      ylab = "A ancestry - logistic",
      xlab = "latitude",
      col = "black", 
      ylim = c(0, 1),
      add = T,
      lwd = 2, 
      lty = 2)
curve(logistic(a0$mu[i] + a0$b[i]*
                 (abs(x) - mean(d_A$abs_lat[d_A$group == "AR_2018"]))), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      ylab = "A ancestry - logistic",
      xlab = "latitude",
      col = "blue", 
      ylim = c(0, 1),
      add = T,
      lwd = 2, 
      lty = 2)
}

i0 <- 31104
d0 <- data.frame(A = unname(t(A[i0, meta.AR.order.by.lat$population])),
                 lat = meta.AR.order.by.lat$lat)
start0 <- getInitial(A ~ SSlogis(lat, Asym, 
                                 xmid, scal),
                     control = nls.control(#tol
                       minFactor = .00001),
                     d = d0)
fit0 <- nls(A ~ logistic3(x = lat, b = b, mu = mu),
            start = list(b = 1/unname(start0["scal"]),
                         mu = unname(start0["xmid"])),
            data = snp,
            trace = F)
fit0cc <- nls(A ~ logistic3(x = lat, b = b, mu = mu),
            start = list(b = 1/unname(cc["b"]),
                         mu = unname(cc["a"])),
            data = d0,
            trace = F)
fit0cc
sum0 <- summary(fit0)


x <- d0$lat
z <- d0$A
z <- z/(1.05 * max(z))  ## scale to max=1/(1.05)
zlogit <- log(z/(1-z))
cc <- setNames(coef(lm(x~zlogit)),c("a","b"))
#This is the linear function that SSlogis() fits (!!)
qplot(zlogit, x) + 
  geom_smooth(method="lm", se=FALSE, colour="red")
plot(zlogit, x)
predfun <- function(x) {
  with(as.list(cc),
       plogis((x - a)/b)*1.05*max(z))
}
ggplot(d0, aes(x = lat, y = A)) + 
  geom_point() +
  stat_function(fun = predfun, colour="red")
ggplot(d0, aes(x = lat, y = A)) + 
  geom_point() +
  stat_function(fun = logistic3(x = d0$lat, 
                                b = coef(fit0cc)["b"],
                                mu = coef(fit0cc)["mu"]), colour="red")
plot(meta.AR.order.by.lat$lat, 
     A[i0, meta.AR.order.by.lat$population],
     ylim = c(0, 1),
     xlim = range(d_A$lat[d_A$group == "AR_2018"]),
     ylab = "A ancestry - logistic",
     xlab = "latitude",
     main = paste("snp", i0))
curve(logistic3(mu = coef(fit0cc)["mu"], 
                b = coef(fit0cc)["b"], 
                x), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      ylab = "A ancestry - logistic",
      xlab = "latitude",
      col = "orange", 
      ylim = c(0, 1),
      add = T,
      lwd = 2, 
      lty = 2)
# maybe for tricky loci I should start with many initial values and then pick the one with lowest error


### --------------------- INDIVIDUAL SNP CLINES --------------------------------- #####
mvn_zero <- list(info = read.table("results/ind_snp_nls_multistart_clines/MVNsim_zero_bounded.info", # data simulation
                             header = T, stringsAsFactors = F),
            params = read.table("results/ind_snp_nls_multistart_clines/MVNsim_zero_bounded.params",
                header = T, stringsAsFactors = F))
mvn <- list(info = read.table("results/ind_snp_nls_multistart_clines/MVNsim_bounded.info", # data simulation
                                   header = T, stringsAsFactors = F),
                 params = read.table("results/ind_snp_nls_multistart_clines/MVNsim_bounded.params",
                                     header = T, stringsAsFactors = F))
sets <- paste0("A", 1:5)
clines <- list(info = do.call(rbind, lapply(sets, function(s)
                                            read.table(paste0("results/ind_snp_nls_multistart_clines/", s, ".info"), # data subset
                                 header = T, stringsAsFactors = F))),
               params = do.call(rbind, lapply(sets, function(s)
                 read.table(paste0("results/ind_snp_nls_multistart_clines/", s, ".params"), # data subset
                            header = T, stringsAsFactors = F))))
# check all models converged -- great!
table(mvn$info$isConv)
table(mvn_zero$info$isConv)
table(clines$info$isConv)

# spread of simulated clines
rbind(mutate(mvn_zero$params, data = "MVN_sim_zero_bounded"), 
      mutate(mvn$params, data = "MVN_sim_bounded"), 
      mutate(clines$params, data = "observed")) %>%
  ggplot(., aes(x = estimate, fill = data)) +
  geom_histogram() +
  facet_wrap(~term, scales = "free_x")
qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "mu"], clines$params$estimate[clines$params$term == "mu"])
#qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "mu"], mvn$params$estimate[mvn$params$term == "mu"])
abline(a = 0, b = 1, col = "blue")
qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "b"], clines$params$estimate[clines$params$term == "b"])
qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "b"], mvn$params$estimate[mvn$params$term == "b"])
abline(a = 0, b = 1, col = "blue")

# what is the FDR?
#FDR_values = c(.1, .05, .01)
# make seq of values to calculate FDR at
b_range = seq(min(clines$params$estimate[clines$params$term == "b"]), 
              max(clines$params$estimate[clines$params$term == "b"]), 
              length.out = 1000)
mu_range = seq(min(clines$params$estimate[clines$params$term == "mu"]), 
             max(clines$params$estimate[clines$params$term == "mu"]), 
             length.out = 1000)

# false discovery rate calculations
fdr_high_b <- sapply(b_range, function(b) fdr_1pop_high(a = b, pop = clines$params$estimate[clines$params$term == "b"], 
                                                        sims = mvn_zero$params$estimate[mvn_zero$params$term == "b"]))
fdr_low_b <- sapply(b_range, function(b) fdr_1pop_low(a = b, pop = clines$params$estimate[clines$params$term == "b"], 
                                                        sims = mvn_zero$params$estimate[mvn_zero$params$term == "b"]))
fdr_high_mu <- sapply(mu_range, function(m) fdr_1pop_high(a = m, pop = clines$params$estimate[clines$params$term == "mu"], 
                                                        sims = mvn_zero$params$estimate[mvn_zero$params$term == "mu"]))
fdr_low_mu <- sapply(mu_range, function(m) fdr_1pop_low(a = m, pop = clines$params$estimate[clines$params$term == "mu"], 
                                                      sims = mvn_zero$params$estimate[mvn_zero$params$term == "mu"]))
# combine to set FDR threshold
FDR_clines = data.frame(FDR_values = FDR_values,
                  high_b = sapply(FDR_values, function(p) ifelse(min(fdr_high_b) <= p, # are there any sig snps?
                                                                     max(b_range[fdr_high_b > p], na.rm = T),
                                                                 NA)),
                  low_b = sapply(FDR_values, function(p) ifelse(min(fdr_low_b) <= p, # are there any sig snps?
                                                                min(b_range[fdr_low_b > p], na.rm = T),
                                                                NA)),
                  high_mu = sapply(FDR_values, function(p) ifelse(min(fdr_high_mu) <= p, # are there any sig snps?
                                                                  max(mu_range[fdr_high_mu > p], na.rm = T),
                                                                  NA)),
                  low_mu = sapply(FDR_values, function(p) ifelse(min(fdr_low_mu) <= p, # are there any sig snps?
                                                                 min(mu_range[fdr_low_mu > p], na.rm = T),
                                                                 NA)),
                  stringsAsFactors = F)
# visualize tails:
#plot(mu_range[850:1000], fdr_high_mu[850:1000])
#points(mvn_zero$params$estimate[mvn_zero$params$term == "mu"], rep(.01, 100000), col = "blue")

write.table(FDR_clines, "results/FDRs_ind_snp_nls_multistart_clines.txt", quote = F, col.names = T, row.names = F, sep = "\t")



##########3 ---------------------- MAKE SOME PLOTS OF CLINES ---------------- ##############


# plot a set of clines from simulations:
mean_cline_mvn <- group_by(mvn$params, term) %>%
  summarise(mean = mean(estimate))
plot(meta.AR.order.by.lat$lat, 
     rep(0, length(meta.AR.order.by.lat$lat)),
     col = NULL,
     ylim = c(0, 1),
     xlim = range(d_A$lat[d_A$group == "AR_2018"]),
     ylab = "A ancestry",
     xlab = "latitude")
curve(logistic3(mu = mean_cline_mvn$mean[mean_cline_mvn$term == "mu"], 
                b = mean_cline_mvn$mean[mean_cline_mvn$term == "b"], 
                x), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      add = T,
      lwd = 2, 
      lty = 2)
# sample 10 random clines
set.seed(100)
sample_10 <- sample(1:(nrow(mvn$params)/2), 100, replace = F)
sapply(sample_10, function(i)
       curve(logistic3(mu = mvn$params$estimate[mvn$params$term == "mu" & mvn$params$snp_index == i], 
                       b = mvn$params$estimate[mvn$params$term == "b" & mvn$params$snp_index == i], 
                       x), 
             from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
             to = range(d_A$lat[d_A$group == "AR_2018"])[2],
             n = 100,
             col = "grey",
             add = T,
             lwd = 1, 
             lty = 1))

# plot 1 steep cline from my data vs. 1 steep cline in the simulated data and median snp cline
outlier_clines <- c(which.max(clines$params$estimate[clines$params$term == "b"]),
                    which.min(clines$params$estimate[clines$params$term == "b"]),
                    which.max(clines$params$estimate[clines$params$term == "mu"]),
                    which.min(clines$params$estimate[clines$params$term == "mu"]))
sapply(1:length(outlier_clines), function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = rainbow(4)[i],
        add = T,
        lwd = 1, 
        lty = 1))
outlier_mvn_clines <- c(which.max(mvn$params$estimate[mvn$params$term == "b"]),
                    which.min(mvn$params$estimate[mvn$params$term == "b"]),
                    which.max(mvn$params$estimate[mvn$params$term == "mu"]),
                    which.min(mvn$params$estimate[mvn$params$term == "mu"]))
sapply(1:length(outlier_mvn_clines), function(i)
  curve(logistic3(mu = mvn$params$estimate[mvn$params$term == "mu" & mvn$params$snp_index == outlier_mvn_clines[i]], 
                  b = mvn$params$estimate[mvn$params$term == "b" & mvn$params$snp_index == outlier_mvn_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = rainbow(4)[i],
        add = T,
        lwd = 1, 
        lty = 2))
plot_random_clines <- function(d, n, seed = 200, color = "darkgrey"){
  # set frame
  plot(meta.AR.order.by.lat$lat, 
       rep(0, length(meta.AR.order.by.lat$lat)),
       col = NULL,
       ylim = c(0, 1),
       xlim = range(d_A$lat[d_A$group == "AR_2018"]),
       ylab = "A ancestry",
       xlab = "latitude")
  # set seed
  set.seed(seed)
  # sample data
  sample_data <- sample(1:(nrow(d$params)/2), 100, replace = F)
  # plot
  sapply(sample_data, function(i)
    curve(logistic3(mu = d$params$estimate[d$params$term == "mu" & d$params$snp_index == i], 
                    b = d$params$estimate[d$params$term == "b" & d$params$snp_index == i], 
                    x), 
          from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
          to = range(d_A$lat[d_A$group == "AR_2018"])[2],
          n = 100,
          col = "darkgrey",
          add = T,
          lwd = 1, 
          lty = 1))
}
# just the background random clines
png("plots/random_clines_100.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
dev.off()
viridis::magma(n = 3)

png("plots/random_clines_100_plus_steep.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
sapply(1, function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1))
dev.off()

png("plots/random_clines_100_plus_introgress.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
sapply(2:3, function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1))
dev.off()

# all 3 outliers
png("plots/random_clines_100_plus_3outliers.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
}
  )
dev.off()

# ADD Points to all of these plots!
png("plots/random_clines_100_plus_steep_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
sapply(1, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(meta.AR.order.by.lat$lat, 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
  })
dev.off()

png("plots/random_clines_100_plus_introgress_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
sapply(2:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(meta.AR.order.by.lat$lat, 
       A[outlier_clines[i], meta.AR.order.by.lat$population],
       pch = 20,
       col = viridis::viridis(n = 3)[i])
})
dev.off()

# all 3 outliers
png("plots/random_clines_100_plus_3outliers_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(meta.AR.order.by.lat$lat, 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
}
)
dev.off()


# what does the N. American cline look like for the 'steep' SNP? Also looks steep.
# should exclude Avalon -- such a consistent outlier point can really affect cline shape
points(-meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[1], meta.pop$zone == "N. America"],
       pch = 20,
       col = "deeppink")

fit_pink = nls_multstart(A ~ logistic3(x = lat, 
                                       b = b, mu = mu),
                                       start_lower = list(b = -1, mu = 25),
                                       start_upper = list(b = 1, mu = 40),
                                       supp_errors = 'Y',
                                       iter = 250,
                                       convergence_count = 100,
                                       data = data.frame(A = unname(t(A[outlier_clines[1], 
                                                                        meta.CA.order.by.lat$population[meta.CA.order.by.lat$population != "Avalon_2014"]])), 
                                                         lat = meta.CA.order.by.lat$lat[meta.CA.order.by.lat$population != "Avalon_2014"], 
                                                         stringsAsFactors = F))
curve(logistic3(mu = unlist(-tidy(fit_pink)[tidy(fit_pink)$term == "mu", "estimate"]), 
                b = unlist(-tidy(fit_pink)[tidy(fit_pink)$term == "b", "estimate"]), 
                x), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      col = "deeppink",
      add = T,
      lwd = 3, 
      lty = 1)





# I can alternatively plot absolute latitude:
curve(logistic3(mu = -mean_cline_mvn$mean[mean_cline_mvn$term == "mu"], 
                b = -mean_cline_mvn$mean[mean_cline_mvn$term == "b"], 
                x), 
      from = abs(range(d_A$lat[d_A$group == "AR_2018"])[2]), 
      to = abs(range(d_A$lat[d_A$group == "AR_2018"])[1]),
      n = 100,
      col = "blue",
      add = F,
      lwd = 2, 
      lty = 2)

# the top clines with data from N. America too:
meta.pop.A <- meta.pop %>%
  mutate(meanA = apply(A[ , meta.pop$population], 2, mean)) %>%
  mutate(abs_lat = abs(lat))
png("plots/high_A_SA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$lat[meta.pop.A$zone == "S. America"], 
       meta.pop.A$meanA[meta.pop.A$zone == "S. America"],
       ylim = c(0, 1),
       pch = 20,
     cex = 1,
       col = col_blind[1],
        ylab = "A ancestry",
        xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "S. America"], 
       A[outlier_clines[2], meta.pop$zone == "S. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["S. America"])
dev.off()

png("plots/high_A_NA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$abs_lat[meta.pop.A$zone == "N. America"], 
       meta.pop.A$meanA[meta.pop.A$zone == "N. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[2], meta.pop$zone == "N. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["N. America"])
dev.off()

# low SA outlier:
png("plots/low_A_SA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$lat[meta.pop.A$zone == "S. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "S. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "S. America"], 
       A[outlier_clines[3], meta.pop$zone == "S. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["S. America"])
dev.off()

png("plots/low_A_NA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$abs_lat[meta.pop.A$zone == "N. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "N. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[3], meta.pop$zone == "N. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["N. America"])
dev.off()






##### -------------------- ENRICHMENT OF STEEP CLINES IN REGIONS OF THE GENOME WITH LOW RECOMB? ---------------######


# re-do calculation of enrichment for steep clines in regions of the genome with low recombination rates
load("../local_ancestry/results/sites_r.RData") # load recombination rate data by sites
steepest_clines <- which(clines$params$estimate[clines$params$term == "b"] >= 
                           quantile(clines$params$estimate[clines$params$term == "b"], .99))
table(sites_r$r_bin5[steepest_clines])

cbind(sites_r, clines$params[clines$params$term == "b", ]) %>%
  filter(estimate >= quantile(estimate, .99)) %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = estimate, color = r_bin5_factor)) +
  facet_wrap(~chr)

# (!) Here high b is steep, but for slopes the other way (on + lat), high b will be large & negative
clines_plus_sites <- sites_r %>%
  mutate(b = clines$params[clines$params$term == "b", "estimate"]) %>%
  mutate(mu = clines$params[clines$params$term == "mu", "estimate"]) %>%
  mutate(top_b = b >= quantile(b, .99)) %>%
  mutate(top_mu = mu >= quantile(mu, .99)) %>%
  mutate(bottom_mu = mu <= quantile(mu, .01))
mean_b_by_r <- clines_plus_sites %>%
  group_by(r_bin5_factor) %>%
  summarise(mean_b = mean(b))
clines_plus_sites %>%
  group_by(r_bin5_factor) %>%
  summarise(mean_b = mean(mu))
clines_plus_sites %>%
  group_by(r_bin5_factor) %>%
  summarise(perc_top_b = mean(top_b),
            perc_top_mu = mean(top_mu),
            perc_bottom_mu = mean(bottom_mu))
b_by_r <- with(clines_plus_sites, lm(b ~ cM_Mb))
summary(b_by_r)
mu_by_r <- with(clines_plus_sites, lm(mu ~ cM_Mb))
summary(mu_by_r)
# subsample to plot faster
dplyr::sample_frac(clines_plus_sites, size = .1) %>%
  ggplot(., aes(x = cM_Mb, y = b)) +
  geom_point() +
  geom_smooth()
# ok so the slopes go down with higher recomb, but is this a big enough effect ot be interesting? not sure. let's plot:
#mean_b_by_r
png("plots/lower_r_regions_have_very_slightly_steeper_mean_slopes.png", height = 6, width = 6, units = "in", res = 300)
plot_random_clines(d = clines, n = 100, seed = 200, color = "darkgrey")
for (i in 1:5){
  curve(logistic3(mu = mean(clines_plus_sites$mu), 
                  b = mean_b_by_r$mean_b[i], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::magma(5)[i],
        add = T,
        lwd = 2, 
        lty = 1)
}
dev.off()

# Group11 possibly has some large inversion or something, not the largest peak, but the widest.
# also not colocalized with the region of high M
cbind(sites_r, clines$params[clines$params$term == "b", ]) %>%
  mutate(AR = apply(A[ , meta.AR.order.by.lat$population], 1, mean)) %>%
  filter(estimate >= quantile(estimate, .99)) %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = AR, color = r_bin5_factor)) +
  facet_wrap(~chr)
