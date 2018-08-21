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
# for scale bar
library(ggsn)

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
extractToDo_CA <- c("CA1410", "CA1315", "CA1211", "CA1120", "CA1013", 
                    "CA0816", "CA0509", "CA0412")

# 8 bees total from feral swarms (enjambres) -- marked for extraction
extractEnjambre <- c("AR1002", "AR1403", "AR1511", "AR2002", 
                     "AR2202", "AR2512", "AR2603", "AR2913")
extractToDo_AR <- c("AR0302", "AR0311", "AR0501", "AR0511", "AR0602", "AR0617",
                    "AR0804", "AR0814", "AR0902", "AR0910", "AR1016", "AR1103", 
                    "AR1115", "AR1201", "AR1216", "AR1303", "AR1310", "AR1410", 
                    "AR1506", "AR1603", "AR1612", "AR1910", "AR1915", "AR2009", 
                    "AR2212", "AR2305", "AR2316", "AR2502", "AR2614", "AR2711", "AR2704")
# decided to do two bees from pop 14 and none from pop15 in this first round of sequencing;
# similarly in Argentina in the Southern most samples, 2 bees from pop 3 and none from pop 1
# decided to seq. one bee from the Southern end of Pop 27, one from the Northern end, and just 1 from pop 29 (enjambre)
# dropped "AR2701" and "AR2909" .. so the 3 bees are spaced between the 2 pops for maximum information across latitude
bees$toSequence <- FALSE                 
bees$toSequence[bees$Bee_ID %in% c(extract1, extract2, extractEnjambre,
                                   extractToDo_CA, extractToDo_AR)] <- TRUE


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
map(database = "state", fill = T, projection = "fisheye", par = "0",
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
              source = "stamen",
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
              source = "stamen",
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
             col = "black") +
  ggsave(filename = "samples_CA.png", device = png(),
         plot = last_plot(), 
         width = 8, height = 8, units = "in",
         dpi = 300)

# color by individual - add individuals to be sequenced
ggmap(get_map(location = 'Ojai, California', 
              maptype = "toner-lite",
              source = "stamen",
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
             col = "black") #+
#  geom_point(aes(x = long,
#                 y = lat),
#             data = bees[bees$Bee_ID %in% extractToDo_CA,],
#             cex = 1,
#             pch = 5,
#             col = "red")

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
# finding pairs
#pairs = lapply(unique(bees[bees$state=="AR", "popN"]), function(N) 
#  getPairFar(bees, popN = N, state = "AR", minDist = 4, maxDist = 20))
#sapply(pairs, function(i) distGeo(bees[bees$Bee_ID == i[1], c("long", "lat")],
#                                  bees[bees$Bee_ID == i[2], c("long", "lat")])/1000)
#ids_in_pairs = unlist(pairs)
#extractToDo = c(extractToDo_CA, ids_in_pairs)
#bees[bees$Bee_ID %in% extractToDo, "toSequence"] <- T
# skip pop1 for now (furthest South in Argentina; along with furthest North in CA)
#bees[bees$popN == "01" & bees$state == "AR" & bees$toSequence == T, "toSequence"] <- F
# combine pops 27 & 29 in Argentina (furthest North); omit AR2906 and AR2716 for now
#bees[bees$Bee_ID %in% c("AR2906", "AR2716"), "toSequence"] <- F
# the bees from pairs are in the list above extractToDo_AR; 
# I replaced by hand some bees randomly chosen with those from enjambres (n=8)

# distances between CA bees
#distm(x = bees[bees$state == "CA" & bees$Bee_ID %in% c(extract1, extract2, extractToDo), 
#               c("long", "lat")], fun = distGeo)/1000
# distances between CA bee pairs to be sequenced
sapply(sort(unique(bees[bees$state == "CA" & bees$toSequence, "popN"])), 
       function(N) distm(bees[bees$toSequence & bees$state == "CA" & 
                                bees$popN == N, c("long", "lat")], 
                         fun = distGeo)/1000)
# distances all AR bees to be sequenced
lapply(sort(unique(bees[bees$state == "AR" & bees$toSequence, "popN"])), 
       function(N) distm(bees[bees$toSequence & bees$state == "AR" & 
                                bees$popN == N, c("long", "lat")], 
                         fun = distGeo)/1000)
distm(x = bees[bees$state == "AR" & bees$toSequence & bees$popN %in% c("01", "02", "27", "29"), 
               c("long", "lat")], fun = distGeo)/1000

# Argentina color by pop -- add individuals to be sequenced
ggmap(get_map(location = 'Santa Fe, Argentina', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  6)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
             data = bees[bees$state == "AR",],
             cex = 1,
             alpha = .5) +
  geom_point(aes(x = long,
                 y = lat),
             data = bees[bees$toSequence == T,],
             cex = 1,
             pch = 5,
             col = "black") +
  geom_point(aes(x = long, # color red bees sampled directly from feral hives
                 y = lat),
             data = bees[bees$Bee_ID %in% extractEnjambre,],
             cex = 1,
             pch = 5,
             col = "red") +
  ggsave(filename = "samples_argentina_with_enjambres.png", device = png(),
         plot = last_plot(), 
         width = 8, height = 8, units = "in",
         dpi = 300)
# look at most Northern pops to decide how to combine -- taking middle 2 of randomly selected ind's from each pop
ggmap(get_map(location = 'Villa Ocampo, Argentina', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = indN), 
             data = bees[bees$state == "AR" & bees$popN %in% c("26", "27", "29"),],
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
write.table(extract1, "../labwork/ids2extract_1", quote = F)
write.table(extract2, "../labwork/ids2extract_2", quote = F)
extract3 = sample(bees[bees$toSequence & ! bees$Bee_ID %in% c(extract1, extract2), "Bee_ID"], 24)
write.table(extract3, "../labwork/ids2extract_3", quote = F)
extract4 = bees[bees$toSequence & ! bees$Bee_ID %in% c(extract1, extract2, extract3), "Bee_ID"]
write.table(extract4, "../labwork/ids2extract_4", quote = F)



# map of CA Africanize Honey bees from Lin 2018 PlosOne paper
plos <- read.csv("Lin2018_plosOne_AHB_CA_supplement_GPS.csv", stringsAsFactors = F)
plos1 = plos[ , 1:5]
colnames(plos1) = c("County", "lat", "long", "AHBs", "Tot")
plos1$p = plos1$AHBs/plos1$Tot
minsec = plos1[grepl("'", plos1$lat, fixed =  T), ]
minsec_approx = mutate(minsec, lat = as.numeric(substr(lat, 1, 5))) %>% 
  mutate(., long = -1*as.numeric(substr(long, 1, 5)))
grepl("\"", plos1$lat, fixed =  T) # has a degree/min/sec notation
plos2 = plos1[!grepl("\"", plos1$lat, fixed =  T) & ! plos1$lat %in% c(" ", "N/A", "n/a"), ] %>% # for those with decimal lat/long, get as numeric
  mutate(lat = abs(as.numeric(lat))) %>%
  mutate(long = -1*abs(as.numeric(long)))
# visualize approximate coordinates of the minsec entries:
ggmap(get_map(location = 'Hollister, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = 1-p,
                 cex = Tot), 
             data = minsec_approx,
             alpha = .5) +
  scale_color_distiller(palette = "RdPu")

ggmap(get_map(location = 'Hollister, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = 1-p,
                 cex = Tot), 
             data = plos2,
             alpha = .5) +
  scale_color_distiller(palette = "RdPu")


# visualize bees Africanization (=p but easier to see 1-p scale)
ggmap(get_map(location = 'Monterey, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  7)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = 1-p,
                 cex = Tot), 
             data = plos2,
             alpha = .5) +
  ggtitle("Africanized bees CA (Lin 2018); p = proportion European mtDNA") + 
  scale_color_distiller(palette = "RdPu") +
  ggsave(filename = "CA_mtDNA_Lin_2018.png", device = png(),
         plot = last_plot(), 
         width = 8, height = 8, units = "in",
         dpi = 300)
  
#scale_color_gradientn(colours = terrain.colors(15))


# zoom in Salinas
ggmap(get_map(location = 'Salinas, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = 1-p,
                 cex = Tot), 
             data = plos2,
             alpha = .5) +
  scale_color_distiller(palette = "RdPu")

# zoom in San Luis Obispo
ggmap(get_map(location = 'San Luis Obispo, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = 1-p,
                 cex = Tot), 
             data = plos2,
             alpha = .5) +
  scale_color_distiller(palette = "RdPu")

# zoom in Lompoc
ggmap(get_map(location = 'Lompoc, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  8)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = 1-p,
                 cex = Tot), 
             data = plos2,
             alpha = .5) +
  scale_color_distiller(palette = "RdPu")

# map of world with all samples
ggmap(get_map(location = 'Costa Rica', 
              #maptype = "satellite", #nice!
              maptype = "watercolor",
              zoom =  3)) +
  geom_point(aes(x = long,
                 y = lat), #,
                 #col = popN), 
             data = bees,
             cex = .25,
             col = "black",
             alpha = .5) +
  ggsave(filename = "samples_world_watercolor.png", device = png(),
         plot = last_plot(), 
         width = 10, height = 10, units = "in",
         dpi = 300)


# zoom in samples already collected socal
ggmap(get_map(location = 'El Cajon, California', 
             maptype = "roadmap",
             # source = "stamen",
              zoom =  9)) +
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
  scalebar(data = bees, dist = 5, dd2km = TRUE, model = "WGS84")


