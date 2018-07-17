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
           remove = F) 
bees <- bees %>%
  # convert degree/min/sec to decimal-based lat & long while keeping those already converted
  # which have min/sec = NA or blank in data
  mutate(., lat = ifelse(is.na(Min_Lat), Deg_Lat, Deg_Lat + Min_Lat/60 + Sec_Lat/3600) * 
           ifelse(bees$state=="AR", -1 , 1)) %>%
  mutate(., long = ifelse(is.na(Min_Long), Deg_Long,
                          Deg_Long + Min_Long/60 + Sec_Long/3600) * -1)

# add in which ones have had DNA extracted
extract1 <- c("CA0906", "CA0914", "CA0102", "CA0108")
extract2 <- c("CA0201", "CA0303", "CA0401", "CA0502", "CA0205",
              "CA0805", "CA1207", "CA1106", "CA1007", "CA1309",
              "CA1402", "CA0308")
#planned <- 
bees$toSequence <- FALSE
bees$toSequence[bees$Bee_ID %in% c(extract1, extract2)] <- TRUE


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
# color by pop
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
             alpha = .5) +
  geom_point(aes(x = long,
                 y = lat),
             data = bees[bees$toSequence == T,],
             cex = 1,
             pch = 5,
             col = "black")

extractToDo_CA <- c("CA1416", "CA1315", "CA1211", "CA1120", "CA1013", "CA0816", "CA0509", "CA0412")

# color by individual - add individuals to be sequenced
ggmap(get_map(location = 'Ojai, California', 
              maptype = "toner-lite",
              zoom =  7)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = indN), 
             data = bees[bees$state == "CA",],
             cex = 1,
             #col = "black",
             alpha = .5) +
  geom_point(aes(x = long,
                 y = lat),
             data = bees[bees$toSequence == T,],
             cex = 1,
             pch = 5,
             col = "black") +
  geom_point(aes(x = long,
                 y = lat),
             data = bees[bees$Bee_ID %in% extractToDo,],
             cex = 1,
             pch = 5,
             col = "red")

# bees sequence Argentina -- too many to check by hand if they are too close together
#distance matrix of all bees in a population (in km, ideally > 4km or at least 3.21 (~2miles)):
getPairFar <- function(data, popN, state, minDist, maxDist){ # takes in all data, a pop number, state
  # and returns a pair of individuals meeting minDist in km away from each other
  X = distm(x = bees[bees$popN == popN & bees$state == state, c("long", "lat")], fun = distGeo)/1000
  #colnames(X) <- bees[bees$popN == popN & bees$state == state, "Bee_ID"]
  #rownames(X) <- bees[bees$popN == popN& bees$state == state, "Bee_ID"]
  Xfar = which(X >= minDist & X <= maxDist, arr.ind=T)
  bees[bees$popN == popN & bees$state == state, "Bee_ID"][Xfar[sample(1:nrow(Xfar), 1),]]
}
pairs = lapply(unique(bees[bees$state=="AR", "popN"]), function(N) 
  getPairFar(bees, popN = N, state = "AR", minDist = 4, maxDist = 20))
sapply(pairs, function(i) distGeo(bees[bees$Bee_ID == i[1], c("long", "lat")],
                                  bees[bees$Bee_ID == i[2], c("long", "lat")])/1000)
ids_in_pairs = unlist(pairs)
extractToDo = c(extractToDo_CA, ids_in_pairs)
bees[bees$Bee_ID %in% extractToDo, "toSequence"] <- T
# skip pop1 for now (furthest South in Argentina; along with furthest North in CA)
bees[bees$popN == "01" & bees$state == "AR" & bees$toSequence == T, "toSequence"] <- F
# combine pops 27 & 29 in Argentina (furthest North); omit AR2906 and AR2716 for now
bees[bees$Bee_ID %in% c("AR2906", "AR2716"), "toSequence"] <- F

# distances between CA bees
#distm(x = bees[bees$state == "CA" & bees$Bee_ID %in% c(extract1, extract2, extractToDo), 
#               c("long", "lat")], fun = distGeo)/1000
# distances between CA bee pairs to be sequenced
sapply(unique(bees[bees$state == "CA" & bees$toSequence, "popN"]), 
       function(N) distm(bees[bees$toSequence & bees$state == "CA" & 
                                bees$popN == N, c("long", "lat")], 
                         fun = distGeo)/1000)[2,]
# Argentina color by individual -- add individuals to be sequenced
ggmap(get_map(location = 'Santa Fe, Argentina', 
              maptype = "toner-lite",
              zoom =  6)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
             data = bees[bees$state == "AR",],
             cex = 1,
             #col = "black",
             alpha = .5) +
  geom_point(aes(x = long,
                 y = lat),
             data = bees[bees$toSequence == T,],
             cex = 1,
             pch = 5,
             col = "red")
# look at most Northern pops to decide how to combine -- taking middle 2 of randomly selected ind's from each pop
ggmap(get_map(location = 'Villa Ocampo, Argentina', 
              maptype = "toner-lite",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
             data = bees[bees$state == "AR",],
             cex = 1,
             #col = "black",
             alpha = .5) +
  geom_point(aes(x = long,
                 y = lat),
             data = bees[bees$toSequence == T,],
             cex = 1,
             pch = 5,
             col = "red")
# which bees to extract next:
extract3 = sample(bees[bees$toSequence & ! bees$Bee_ID %in% c(extract1, extract2), "Bee_ID"], 24)
write.table(extract3, "ids2extract_3", quote = F)
extract4 = bees[bees$toSequence & ! bees$Bee_ID %in% c(extract1, extract2, extract3), "Bee_ID"]
write.table(extract4, "ids2extract_4", quote = F)
