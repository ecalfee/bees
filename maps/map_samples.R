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

my_api <- read.table("google_maps_api_EC2018.txt",
                     header = F, stringsAsFactors = F)$V1
ggmap::register_google(key = my_api)
# import bees from field notes
ca_bees <- read.csv("~/Dropbox/grad school/Research/bee project/CA_Bees/CA_bee_fieldnotes.csv", 
                      header = T,
                      stringsAsFactors = F)
ar_bees <- read.csv("~/Dropbox/grad school/Research/bee project/Argentina/Field Notes Argentina bees/Argentina_fieldnotes_bee_gps_formatted.csv",
                      header = T,
                      stringsAsFactors = F)
mx_bees <- read.csv("~/Dropbox/grad school/Research/bee project/Mexico/Baja_bee_fieldnotes.csv",
                    header = T,
                    stringsAsFactors = F) %>%
  mutate(Time = as.character(Time)) %>%
  mutate(Plant = as.character(Plant))
# combine data sets
bees_raw <- full_join(ca_bees, ar_bees) %>%
  full_join(., mx_bees) %>%
  separate(., col = Bee_ID, sep = c(2,4), 
           into = c("state", "popN", "indN"),
           remove = F) %>%
  filter(., Bee_ID!="") # exclude empty records
bees <- bees_raw %>%
  # convert degree/min/sec to decimal-based lat & long while keeping those already converted
  # which have min/sec = NA or blank in data
  mutate(., lat = ifelse(is.na(Min_Lat), Deg_Lat, Deg_Lat + Min_Lat/60 + Sec_Lat/3600) * 
           ifelse(state=="AR", -1 , 1)) %>%
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
# DNA extraction groups
# which bees to extract next:
#write.table(extract1, "../labwork/ids2extract_1", quote = F)
#write.table(extract2, "../labwork/ids2extract_2", quote = F)
#extract3 = sample(bees[bees$toSequence & ! bees$Bee_ID %in% c(extract1, extract2), "Bee_ID"], 24)
#write.table(extract3, "../labwork/ids2extract_3", quote = F)
extract3 = read.table("../labwork/ids2extract_3", stringsAsFactors = F)$x
#extract4 = bees[bees$toSequence & ! bees$Bee_ID %in% c(extract1, extract2, extract3), "Bee_ID"]
#write.table(extract4, "../labwork/ids2extract_4", quote = F)
extract4 = read.table("../labwork/ids2extract_4", stringsAsFactors = F)$x

# I decided to do two bees from pop 14 and none from pop15 in this first round of sequencing for CA;
# similarly in Argentina in the Southern most samples, 2 bees from pop 3 and none from pop 1
# decided to seq. one bee from the Southern end of Pop 27, one from the Northern end, and just 1 from pop 29 (enjambre)
# dropped "AR2701" and "AR2909" .. so the 3 bees are spaced between the 2 pops for maximum information across latitude
bees$toSequence <- FALSE                 
bees$toSequence[bees$Bee_ID %in% c(extract1, extract2, extractEnjambre,
                                   extractToDo_CA, extractToDo_AR)] <- TRUE
bees$enjambre <- bees$Bee_ID %in% extractEnjambre # bee from a feral colony "enjambre"
# a few bees I dried to see how that went getting 
# a dry weight -- will not use these bees for DNA
bees_dried <- c("CA0101", "CA0104", "CA0107", "CA0110",
                "CA0406", "CA0415", "CA0909", "CA0911")
bees$dried <- bees$Bee_ID %in% bees_dried

# assign bees to sequence (which lane) & which not to sequence as well as to the morphology study
# load set of random bees with wings already measured
wing_done <- unlist(read.table("../labwork/ids2morph_10N_CA", stringsAsFactors = F,
                               sep = " ", header = F))
bees$wingPass1 <- bees$Bee_ID %in% wing_done
set.seed(341)
bees_random_order_CA <- filter(bees, state == "CA") %>%
  .[sample(1:nrow(.)), ] %>% # randomize order so bees w/ duplicate gps coordinates are selected out at random
  arrange(., !wingPass1) %>%
  arrange(., !toSequence) %>% # always include bees already specified to sequence
  mutate(., duplicate = duplicated(.[ , c("lat", "long")]))
no_dup_CA <- filter(bees_random_order_CA, !duplicate)
n_wing = 12 # select twelve bees per population for wing analysis
wing_to_do_CA <- unlist(lapply(unique(no_dup_CA$popN), function(pop)
  no_dup_CA[no_dup_CA$popN == pop, 
            "Bee_ID"][1:n_wing]))
no_dup_CA$toWing <- no_dup_CA$Bee_ID %in% wing_to_do_CA
# from the subset of bees for wing analysis,
# select 8-9 bees per population for sequencing
n_seq_CA = data.frame(popN = sort(unique(no_dup_CA$popN)),
                      n_seq = c(8,9,8,8,9,8,8,9,8,8,9,8,0)) # skip population #15 for sequencing
# random subset argentina bees:
set.seed(342)
bees_random_order_AR <- filter(bees, state == "AR") %>%
  .[sample(1:nrow(.)), ] %>% # randomize order so bees w/ duplicate gps coordinates are selected out at random
  arrange(., !wingPass1) %>%
  arrange(., !toSequence) %>% # always include bees already specified to sequence
  mutate(., duplicate = duplicated(.[ , c("lat", "long")]))
no_dup_AR <- filter(bees_random_order_AR, !duplicate, Bee_ID != "AR1515") %>%# better to exclude this bee without full sampling information
  mutate(., toWing = toSequence) # mark bees sequenced as ones with wings to do

# select 8-9 bees per population for sequencing
n_seq_AR = data.frame(popN = sort(unique(no_dup_AR$popN)),
                      n_seq = c(9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9))

# e.g. pop 13 -- how can I choose 8 bees to maximize distances between them?
# actually what I want to do is choose the bees that minimize
# the number of pairwise distances < 1.6km (1 mile)
# CA bees I can likely limit to all pairs > 3.2km apart (2 miles)
# I could additionally set a maximum distance w/in pop to 30km
# After applying these restrictions, I want to choose bees randomly
# For CA, see if I can meet these restrictions w/in the 12 I selected for wing analysis
# For both CA & AR, require that the already sequenced bees be included

# function takes in popN, data (e.g. no_dup_AR), 
# number of bees to select for sequencing from that pop,
# and all bee id's (if any) for bees already sequenced
# returns a random set that minimizes pairwise distances under minDist
# and has no pairs with distances over maxDist
# note: set seed ahead of time
find_set_toSeq <- function(d, N, n, 
                           already_seq = c(extract1, extract2, extract3, extract4), 
                           minDist, maxDist,
                           seed = 100){
  # set seed
  set.seed(seed)
  # subset data to population of interest
  beePop <- filter(d, popN == N)
  # get pairwise distances between all bees in that population
  beeDist <- distm(x = beePop[ , c("long", "lat")], 
                   fun = distGeo)/1000
  diag(beeDist) <- NA
  # which bees (if any) have already been sequenced?
  must_incl <- which(beePop$Bee_ID %in% already_seq)
  #print(paste("already sequenced: ", beePop$Bee_ID[must_incl]))
  # get all random combinations that sum to n total and include 'must_incl' bees
  combo_list <- apply(combn(x = (1:nrow(beeDist))[!(beePop$Bee_ID %in% already_seq)], 
                            m = n - sum(beePop$Bee_ID %in% already_seq)),
                      2, function(i) c(must_incl, i))
  n_too_close <- apply(combo_list, 2,
                       function(i) sum(beeDist[i,i] < minDist, na.rm = T)/2)
  n_too_far <- apply(combo_list, 2,
                     function(i) sum(beeDist[i,i] > maxDist, na.rm = T)/2)
  good_dist <- combo_list[ , which(n_too_far == 0 & n_too_close == min(n_too_close))]
  
  # if there are some bees below min dist, pick from the sets that maximize the close distances
  if (is.null(dim(good_dist)) | min(n_too_close) == 0){
    good_dist2 <- good_dist
  }else{
    # what's the mean of the distances less than minDist?
    mean_small_dist <- apply(good_dist, 2,
                           function(i) mean(beeDist[i,i][beeDist[i,i] < minDist], na.rm = T)/2)
    # find the set with the maximum mean small distances (ie reduce # bee pairs super close when possible)
    good_dist2 <- good_dist[ , which(mean_small_dist == max(mean_small_dist))]
}
  # of the set at a good distance, how many have wing analysis done?
  if (is.null(dim(good_dist2))){
    pick <- good_dist2 # there is only one possible choice
  }else{
    # sample 1 set that's at a good set of distances and has large # already wing analyzed
    n_toWing <- apply(good_dist2, 2,
                      function(i) sum(beePop[i,"toWing"]))
    pick <- good_dist2[ , sample(which(n_toWing == max(n_toWing)), 1)]
  }
  
  return(data.frame(Bee_ID = beePop$Bee_ID[pick],
                    toWing = beePop$toWing[pick],
                    n_too_close = min(n_too_close), 
                    minDist = min(beeDist[pick, pick], na.rm = T),
                    maxDist = max(beeDist[pick, pick], na.rm = T),
                    meanDist = mean(beeDist[pick,pick],na.rm = T),
                    state = unique(beePop$state),
                    popN = N,
                    stringsAsFactors = F))
}

# Argentina bees
all_seq_AR <- do.call(rbind, lapply(1:nrow(n_seq_AR), function(j)
  find_set_toSeq(d = filter(no_dup_AR, !dried), 
                 N = n_seq_AR[j, "popN"],
                 n = n_seq_AR[j, "n_seq"],
                 minDist = 1.6, maxDist = 30)))
#unique(all_seq_AR$Bee_ID)
#group_by(all_seq_AR, popN) %>% summarise(mean(n_too_close)) %>% print(., n = 100)
#group_by(all_seq_AR, popN) %>% summarise(sum(toWing)) %>% View(.)
#group_by(all_seq_AR, popN) %>% summarise(mean(meanDist)) %>% View(.)
# CA bees
all_seq_CA <- do.call(rbind, lapply(1:(nrow(n_seq_CA)-1), function(j)
  find_set_toSeq(d = filter(no_dup_CA, !dried), 
                 N = n_seq_CA[j, "popN"],
                 n = n_seq_CA[j, "n_seq"],
                 minDist = 3.2, maxDist = 30)))
#group_by(all_seq_CA, popN) %>% summarise(mean(n_too_close)) %>% View(.)
#group_by(all_seq_CA, popN) %>% summarise(mean(maxDist)) %>% View(.)
#group_by(all_seq_CA, popN) %>% summarise(mean(minDist)) %>% View(.)
#group_by(all_seq_CA, popN) %>% summarise(mean(meanDist)) %>% View(.)
#group_by(all_seq_CA, popN) %>% summarise(sum(toWing)) %>% View(.)

# add labels to main bee dataframe
bees$duplicate <- bees$Bee_ID %in% 
  filter(rbind(bees_random_order_CA, bees_random_order_AR), 
         duplicate)$Bee_ID
bees$toWing <- bees$Bee_ID %in% c(wing_to_do_CA, all_seq_CA$Bee_ID, all_seq_AR$Bee_ID)
bees$toSequence <- bees$Bee_ID %in% c(all_seq_CA$Bee_ID, all_seq_AR$Bee_ID)

# lanes to do next (lanes 2-5 w/ 54 bees each)
# w/ populations evenly distributed across lanes
# and assigned in extraction groups to do next (sets of 24, groups #5-13)
set.seed(360)
bees_with_lanes <- bees[sample(1:nrow(bees)), ] %>% # randomize order
  filter(!(Bee_ID %in% c(extract1, extract2, extract3, extract4)) & toSequence) %>% # only bees to be sequenced & not already done in lane 1
  arrange(., popN) %>%
  arrange(., state) %>%
  mutate(lane = rep(2:5, 54)) %>%
  .[sample(1:nrow(.)), ] %>% # randomize order again
  arrange(., lane) %>%
  mutate(extraction = unlist(lapply(5:13, function(i) rep(i, 24)))) %>%
  select(Bee_ID, popN, state, lane, extraction) %>%
  left_join(bees, ., by = c("Bee_ID", "state", "popN"))

# add in extraction # and lane for first set of 63 bees:
bees_with_lanes$extraction[bees_with_lanes$Bee_ID %in% extract1] <- 1
bees_with_lanes$extraction[bees_with_lanes$Bee_ID %in% extract2] <- 2
bees_with_lanes$extraction[bees_with_lanes$Bee_ID %in% extract3] <- 3
bees_with_lanes$extraction[bees_with_lanes$Bee_ID %in% extract4] <- 4

bees_with_lanes$lane[bees_with_lanes$extraction %in% 1:4] <- 1 # lane one had bee extractions 1-4


# ok bees are evenly distributed over remaining lanes:
bees_with_lanes %>%
  group_by(lane, state) %>%
  summarise(n())
bees_with_lanes %>%
  group_by(lane, state) %>%
  summarise(mean(lat))
bees_with_lanes %>%
  group_by(extraction, state, lane) %>%
  summarise(n()) %>%
  View(.)

# add extraction and lane number to main bees database
bees <- left_join(bees,
                  select(bees_with_lanes, c("Bee_ID", "lane", "extraction")),
                  by = "Bee_ID")

# write file with bees to do for wing analysis:
# CA:
wing_print_CA2 <- do.call(cbind, lapply(unique(filter(bees, Bee_ID %in% wing_to_do_CA)$popN),
                                        function(i)
                                          filter(bees, Bee_ID %in% wing_to_do_CA & popN == i)$Bee_ID))
write.table(wing_print_CA2, "../labwork/CA_wing_bee_ids2.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)
wing_print_CA3 <- filter(bees, toWing, state == "CA", !(Bee_ID %in% wing_to_do_CA)) %>%
  arrange(popN, Bee_ID) %>%
  select(Bee_ID)
write.table(wing_print_CA3, "../labwork/CA_wing_bee_ids3.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)
# AR:
wing_print_AR <- filter(bees, toWing, state == "AR") %>%
  arrange(popN, lane != 1, Bee_ID) %>%
  select(Bee_ID, popN)
write.table(wing_print_AR$Bee_ID, "../labwork/AR_wing_bee_ids.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

# make some simple csv files of bee coordinates:
bees %>%
  filter(., state == "AR") %>%
  select(., Bee_ID, lat, long, popN, indN, Date, Time, toSequence, toWing, enjambre) %>%
  write.csv(., "Arg_bee_lat_long_4_google_maps.csv",
            quote = F, row.names = F)
bees %>%
  filter(., state == "CA") %>%
  select(., Bee_ID, lat, long, popN, indN, Date, Time, toSequence, toWing, enjambre) %>%
  write.csv(., "CA_bee_lat_long_4_google_maps.csv",
            quote = F, row.names = F)
bees %>%
  filter(., state == "MX") %>%
  mutate(., toSequence = T) %>%
  mutate(., toWing = T) %>%
  mutate(., enjambre = NA) %>% # unknown status
  select(., Bee_ID, lat, long, popN, indN, Date, Time, toSequence, toWing, enjambre) %>%
  write.csv(., "MX_bee_lat_long_4_google_maps.csv",
            quote = F, row.names = F)

# print files for labwork extractions 9-13:
set.seed(100)
for (x in 5:13){
  bees %>%
    filter(., extraction == x) %>%
    rename(extract = extraction) %>%
    .[sample(1:nrow(.)), ] %>% # randomize order
    arrange(lane) %>% # but keep lanes together
    mutate(n = 1:24) %>%
    select(., lane, extract, n, Bee_ID) %>%
    write.table(., 
                paste0("../labwork/extract_", x, ".txt"),
                quote = F, row.names = F,
                sep = "\t")
}

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
map <- get_map(location = 'Santa Fe, Argentina', 
               zoom =  6)
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
       function(N) mean(distm(bees[bees$toSequence & bees$state == "CA" & 
                                bees$popN == N, c("long", "lat")], 
                         fun = distGeo)/1000))
sapply(sort(unique(bees[bees$state == "CA" & bees$toSequence, "popN"])), 
       function(N) {d <- distm(bees[bees$toSequence & bees$state == "CA" & 
                                      bees$popN == N, c("long", "lat")], 
                               fun = distGeo)/1000
       diag(d) <- NA
       min(d, na.rm = T)})

# distances all AR bees to be sequenced
sapply(sort(unique(bees[bees$state == "AR" & bees$toSequence, "popN"])), 
       function(N) max(distm(bees[bees$toSequence & bees$state == "AR" & 
                                bees$popN == N, c("long", "lat")], 
                         fun = distGeo)/1000))
sapply(sort(unique(bees[bees$state == "AR" & bees$toSequence, "popN"])), 
       function(N) {d <- distm(bees[bees$toSequence & bees$state == "AR" & 
                                    bees$popN == N, c("long", "lat")], 
                             fun = distGeo)/1000
       diag(d) <- NA
       min(d, na.rm = T)})
distm(x = bees[bees$state == "AR" & bees$toSequence & bees$popN %in% c("01", "02", "27", "29"), 
               c("long", "lat")], fun = distGeo)/1000

# distance between all Argentina samples:
ar_dist <- distm(x = bees[bees$state == "AR", 
               c("long", "lat")], fun = distGeo)/1000
ar_dist_1 <- distm(x = bees[bees$state == "AR" & bees$popN == "01",
                            c("long", "lat")], 
                   fun = distGeo)/1000
heatmap(ar_dist_1)
hist(ar_dist_1)

ar_dist_29 <- distm(x = bees[bees$state == "AR" & bees$popN == "29",
                            c("long", "lat")], 
                   fun = distGeo)/1000
heatmap(ar_dist_29)
hist(ar_dist_29)
# attempt to make a non-duplicate GPS location
# set of bees w/ minimum distance for sequencing:
set.seed(8293)
no_dup_AR_dist <- distm(x = no_dup_AR[ ,
                             c("long", "lat")], 
                    fun = distGeo)/1000
heatmap(no_dup_AR_dist, scale = "none")

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

# zoomed in maps of argentina
# for marcelo
arg_map29 <- get_map(location = 'Villa Ocampo, Argentina', 
                     maptype = "terrain",
                     source = "stamen",
                     zoom =  11)
ggmap(arg_map29) + 
  geom_point(aes(x = long, 
                 y = lat,
                 col = popN), 
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
  ggsave(filename = "CA_mtDNA_Lin_2018.png", device = "png",
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
  ggsave(filename = "samples_world_watercolor.png", device = "png",
         plot = last_plot(), 
         width = 10, height = 10, units = "in",
         dpi = 300)

# zoom in samples already collected socal
ggmap(get_map(location = 'El Cajon, California', 
             maptype = "roadmap",
             # source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat), 
             data = bees[bees$state == "CA",],
             cex = 1,
             pch = 5,
             col = "black",
             alpha = .5)

# sample 10 bees per Northern CA pop for initial morphology analysis
#morph_list10 <- sapply(9:15, function(i) sort(sample(bees[bees$state == "CA" & as.integer(bees$popN) == i, "Bee_ID"], size = 12, replace = F)))
# note: In my 12 should be the bees that I chose for sequencing..
#write.table(morph_list10, "../labwork/ids2morph_10N_CA", col.names = F, row.names = F, quote = F)


# plot just the bees successfully sequenced:
# load population ancestry frequencies:
pops.seq <- read.table("../bee_samples_listed/byPop/pops_included.list", stringsAsFactors = F)$V1
bees.seq <- do.call(rbind, 
                lapply(pops.seq, function(p) data.frame(Bee_ID = read.table(paste0("../bee_samples_listed/byPop/", p, ".list"),
                                                                        stringsAsFactors = F)$V1, population = p, stringsAsFactors = F)))
meta.seq <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  left_join(bees.seq, ., by = c("Bee_ID", "population"))
CA_pops_included <- unique(meta.seq$population[meta.seq$group %in% c("N_CA", "S_CA", "CA_2018") & meta.seq$year >= 2014])


ggmap(get_map(location = 'Hollister, California', 
              maptype = "toner-lite",
              source = "stamen",
              zoom =  9)) +
  geom_point(aes(x = long, 
                 y = lat,
                 col = population,
                 shape = factor(year)), 
             data = meta.seq %>%
               filter(population %in% CA_pops_included)) +
  scale_color_distiller(palette = "RdPu")

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
map_CA_wide <- ggmap(get_map(location = 'Palmo, California',
  maptype = "toner-lite",
  source = "stamen",
  zoom = 7))
map_CA_wider <- ggmap(get_map(location = 'California',
                              maptype = "toner-lite",
                              source = "stamen",
                              zoom = 6))
map_CA_wide +
  geom_point(data = filter(meta.seq, population %in% CA_pops_included) %>%
               mutate(Collection = factor(year)),
             aes(x = long, 
                 y = lat,
                 color = population,
                 shape = Collection),
             alpha = .75) +
  guides(color = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Map of Bees Sequenced from California")
ggsave(filename = "plots/map_bees_sequenced_CA.png", 
       device = "png",
       width = 6, height = 6, 
       units = "in",
       dpi = 300)
# now plot on same scale as Argentina's map:
map_CA_wider +
  geom_point(data = filter(meta.seq, population %in% CA_pops_included) %>%
               mutate(Collection = factor(year)),
             aes(x = long, 
                 y = lat,
                 color = population,
                 shape = Collection),
             alpha = .75) +
  guides(color = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Map of Bees Sequenced from California")
ggsave(filename = "plots/map_bees_sequenced_CA_scaled2AR.png", 
       device = "png",
       width = 6, height = 6, 
       units = "in",
       dpi = 300)

# Argentina samples
map_AR_wide <- ggmap(get_map(location = 'Colon, Argentina',
                              maptype = "toner-lite",
                              source = "stamen",
                              zoom = 6))

map_AR_wide +
  geom_point(data = filter(meta.seq, geographic_location == "Argentina") %>%
               mutate(Collection = factor(year)),
             aes(x = long, 
                 y = lat,
                 color = population,
                 shape = Collection),
             alpha = .75) +
  guides(color = FALSE) +
  scale_shape_manual(values=c(17)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Map of Bees Sequenced from Argentina")
ggsave(filename = "plots/map_bees_sequenced_AR.png", 
       device = "png",
       width = 6, height = 6, 
       units = "in",
       dpi = 300)


