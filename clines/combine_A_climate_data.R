library(sp)
library(raster)
library(geosphere)
library(tidyr)
# simple script retrieves climate data and global ancestry data. Creates d_A data object.
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

# load wing length data
wings <- read.table(file = "../wing_analysis/results/bee_wing_lengths_cm_1.28.20.txt",
                    header = T, sep = "\t", stringsAsFactors = F)

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
  meta.ind.included <- meta.ind[meta.ind$population %in% c("Avalon_2014", "Davis_2014", "Placerita_2014", "Riverside_2014",
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
  mutate(year = factor(year)) %>%
  left_join(., wings[ , c("Bee_ID", "wing_cm")], by = "Bee_ID") # add wing length data

# save d_A
save(file = "results/d_A.RData", list = c("d_A"))


# distance between 2 riverside points:
skyvalley <- data.frame(long = -116.27, lat = 33.84)
idyllwild <- data.frame(long = -116.72, lat = 33.74)
distm(skyvalley, idyllwild)/1000 # km, geodesic distance, equivalent to distGeo()
distm(c(max(d_A$long[d_A$population == "Riverside_2014"]), max(d_A$lat[d_A$population == "Riverside_2014"])),
      c(min(d_A$long[d_A$population == "Riverside_2014"]), min(d_A$lat[d_A$population == "Riverside_2014"])))/1000
