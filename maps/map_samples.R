# this script makes a map of GPS locations of individually sampled bees 
# from the Africanized honey bee project
library(dplyr)
library(tidyr)
library(ggplot2)
# R packages to plot GPS
library(maps)
library(mapproj)
library(geosphere)
# another maps package
library(ggmap)

# import bees from field notes
ca_bees <- read.csv("~/Dropbox/grad school/bee project/CA_Bees/CA_bee_fieldnotes.csv", 
                      header = T,
                      stringsAsFactors = F)
ar_bees <- read.csv("~/Dropbox/grad school/bee project/Argentina/Field Notes Argentina bees/Argentina_fieldnotes_bee_gps_formatted.csv",
                      header = T,
                      stringsAsFactors = F)
# combine data sets
bees <- full_join(ca_bees, ar_bees) %>%
  separate(., col = Bee_ID, sep = c(2,4), 
           into = c("state", "popN", "indN"),
           remove = F) %>%
  # convert degree/min/sec to decimal-based lat & long while keeping those already converted
  # which have min/sec = NA or blank in data
  mutate(., lat = ifelse(is.na(Min_Lat), Deg_Lat, Deg_Lat + Min_Lat/60 + Sec_Lat/360) * 
           ifelse(bees$state=="AR", -1 , 1)) %>%
  mutate(., long = ifelse(is.na(Min_Long), Deg_Long,
                          Deg_Long + Min_Long/60 + Sec_Long/360) * -1)

##get a map
par(mar = c(1, 1, 1, 1) + 0.1)
map(database = "world", fill = T, projection = "fisheye", par = "0",
    col = "darkgreen", bg = "darkblue",
    xlim = c(-124, -57), # long
    ylim = c(-38, 38)) # lat


##get points on the map
coord = mapproject(x = bees$long, y = bees$lat) # convert to flat coordinates
points(x = coord$x, y = coord$y, col = "yellow", pch = 1, cex = .01) # add to map

# map of CA bees only
par(mar = c(1, 1, 1, 1) + 0.1)
map(database = "county", fill = T, projection = "fisheye", par = "0",
    col = "darkgreen", bg = "darkblue",
    xlim = c(-124, -118), # long
    ylim = c(28, 38)) # lat

coord = mapproject(x = bees$long, y = bees$lat) # convert to flat coordinates
points(x = coord$x, y = coord$y, col = "yellow", pch = 1, cex = .01) # add to map

# new map with different background
map <- get_map(location = 'Santa Fe, Argentina', zoom =  6)
mapPoints <- ggmap(map) +
  geom_point(aes(x = long, 
                     y = lat,
                 col = popN), 
                     data = bees[bees$state == "AR",],
             cex = .5,
                 alpha = .5)
# plot all of Argentina's samples
mapPoints

# diff background; all Argentina samples
ggmap(get_map(location = 'Santa Fe, Argentina', 
              maptype = "toner-lite",
zoom =  6)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
             data = bees[bees$state == "AR",],
             cex = .5,
             alpha = .5)

# Northern Argentina samples
ggmap(get_map(location = 'Reconquista, Argentina', zoom =  8)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
             data = bees[bees$state == "AR",],
             cex = .5,
             alpha = .5)

# CA samples
ggmap(get_map(location = 'Ojai, California', zoom =  7)) +
  geom_point(aes(x = long, 
                 y = lat), 
             data = bees[bees$state == "CA",],
             cex = .5,
             col = "black",
             alpha = .5)

ggmap(get_map(location = 'Ojai, California', 
              #maptype = "watercolor", #nice!
              maptype = "toner-lite",
              zoom =  7)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
             data = bees[bees$state == "CA",],
             cex = .5,
             #col = "black",
             alpha = .5)
 
