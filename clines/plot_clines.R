# script to plot clines in ancestry, e.g. with latitude and environment
library(sp)
library(raster)
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
write.table(bees.clim, "results/BioclimVar_bees_2.5min.txt", sep = "\t",
            quote = F, col.names = T, row.names = F)

bees.clim %>%
  ggplot(., aes(x = MeanTempColdestQuarter, y = AnnualMeanTemp, color = geographic_location_short)) +
  geom_point()
bees.clim %>%
  ggplot(., aes(x = MeanTempColdestQuarter, y = abs(lat), color = geographic_location_short)) +
  geom_point()
bees.clim %>%
  ggplot(., aes(x = AnnualMeanTemp, y = abs(lat), color = geographic_location_short)) +
  geom_point()

# add admixture.
# do linear models using latitude, temp, and distance from sao paulo brazil as variables.

