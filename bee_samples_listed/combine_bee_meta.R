library(dplyr)
library(tidyr)
# this script combines metadata for already sequenced ca bees (source: Ramirez),
# Kenya bees (source Sheppard NCBI) and A/C/M/Cerana bees (source: Harpur NCBI)

# harpur A C M group bees & one Cerana reference ind.
harpur <- read.table("../data/Harpur_2014_NCBI/Harpur_SraRunTable.txt",
                     sep = "\t", stringsAsFactors = F,
                     header = T)
harpur_incl <- harpur %>%
  dplyr::select(Run_s, geographic_location_s, strain_s) %>%
  mutate(population = ifelse(strain_s == "Apis mellifera mellifera" | strain_s == "Apis mellifera iberiensis", "M",
                        ifelse(strain_s == "Apis mellifera scutellata", "A",
                               ifelse(strain_s == "Apis mellifera yemenitica", "Y",
                                      ifelse(strain_s == "Apis mellifera carnica", "C",
                                             ifelse(strain_s == "Apis cerana", "Cerana", NA)))))) %>%
  filter(population != "Y") %>%
  rename(Bee_ID = Run_s) %>%
  rename(geographic_location = geographic_location_s) %>%
  rename(strain = strain_s) %>%
  mutate(year = 2013) %>%
  mutate(source = "Harpur") %>%
  mutate(group = population) # group is A, M, C, Cerana
# more A bees from Kenya, Sheppard collection
kenya <- read.table("../data/Kenya_Sheppard_NCBI/Kenya_SraRunTable.txt",
                    sep = "\t", stringsAsFactors = F,
                    header = T)
kenya_incl <- kenya %>%
  dplyr::select(Run_s, Organism_s, collection_date_s) %>%
  mutate(geographic_location = "Kenya") %>%
  rename(strain = Organism_s) %>%
  mutate(population = "A") %>%
  rename(Bee_ID = Run_s) %>%
  mutate(year = as.integer(collection_date_s)) %>%
  dplyr::select(-collection_date_s) %>%
  mutate(source = "Sheppard") %>%
  mutate(group = "A")

# CA bees from Ramirez lab
ca_bees <- read.table("../bee_filter.txt",
                      header = F, stringsAsFactors = F) %>%
  .[,1:2] # note: ap50 is a duplicated record -- I fix this below for the .meta file!!
colnames(ca_bees) <- c("population", "ID")
ca_list <- read.table("../data/CA_Bee/samples.list", header = F, stringsAsFactors = F)$V1 # has 2 fewer because Davis1968 didn't transfer and ap50 is duplicated below
ca_bees_incl <- ca_bees %>%
  dplyr::filter(!(ca_bees$population %in% c("M", "A", "C", "Cerana"))) %>%
  separate(., population, c("geographic_location", "year"), remove = F) %>%
  mutate(strain = "unknown") %>%
  mutate(year = as.integer(year)) %>%
  mutate(source = "Ramirez") %>%
  rename(Bee_ID = ID) %>%
  mutate(group = ifelse(geographic_location %in% c("Berkeley", "Stanislaus", "Humboldt", "Stebbins", "Davis", "Domestic"), "N_CA", "S_CA")) %>%
  filter(!duplicated(.))
ca_bees_gps <- read.table("Cridland_CA_bee_GPS.txt",
                          sep = "\t",
                          header = T, stringsAsFactors = F) %>%
  mutate(city = ifelse(city == "", "unknown", city))
riverside_subpops <- read.table("Riverside_subpops.txt",
                               sep = "\t",
                               header = T, stringsAsFactors = F)
ca_bees_incl_gps <- 
  left_join(ca_bees_incl, riverside_subpops, 
          by = c("geographic_location", "Bee_ID")) %>%
  mutate(city = ifelse(is.na(city), "unknown", city)) %>%
  left_join(., ca_bees_gps, by = c("geographic_location", "city")) %>%
  dplyr::select(-city)

# metadata for all reference panel bees from prior studies
# combine and sort into pop. groups before printing to file
all <- bind_rows(harpur_incl, kenya_incl, ca_bees_incl_gps) %>%
.[order(.$population), ] %>%
#.[order(.$year), ] %>%
.[order(.$strain), ]
write.table(all, "harpur_kenya_ca_bees.meta",
            row.names = F, col.names = T, sep = "\t", quote = F)


# make a list of the reference individuals + ca_bees post 1994 (CA introduction of Afr. honeybees)
post_1994 <- all[all$population %in% c("A", "C", "M") |
                   (all$strain == "unknown" & all$year >= 1994), ]
write.table(post_1994, "post_1994.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
write.table(post_1994$Bee_ID, "post_1994.list",
            row.names = F, col.names = F,
            quote = F)

# add in CA and AR bees newly collected:
newCA <- read.table("../maps/CA_bee_lat_long_4_google_maps.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(year = 2018) %>%
  mutate(strain = "unknown") %>%
  mutate(source = "Calfee") %>%
  mutate(population = substr(Bee_ID, 1, 4)) %>%
  mutate(geographic_location = "California") %>%
  mutate(group = "CA_2018")
newAR <- read.table("../maps/Arg_bee_lat_long_4_google_maps.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(year = 2018) %>%
  mutate(strain = "unknown") %>%
  mutate(source = "Calfee") %>%
  mutate(population = substr(Bee_ID, 1, 4)) %>%
  mutate(geographic_location = "Argentina") %>%
  mutate(group = "AR_2018")
newMX <- read.table("../maps/MX_bee_lat_long_4_google_maps.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(year = 2018) %>%
  mutate(strain = "unknown") %>%
  mutate(source = "Calfee") %>%
  mutate(population = substr(Bee_ID, 1, 4)) %>%
  mutate(geographic_location = "Mexico") %>%
  mutate(group = "MX_2019")
write.table(newCA, "CA_2018.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
write.table(newAR, "AR_2018.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
write.table(newMX, "MX_2019.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)


# quick kohn metadata:
kohn_meta <- data.frame(Bee_ID = c("SanDiego001", "SanDiego002", "Mexico001"),
                        source = "Kohn",
                        strain = "unknown",
                        year = 2015,
                        group = "Kohn",
                        geographic_location = c("San Diego", "San Diego", "Mexico"),
                        population = c("SanDiego_2015", "SanDiego_2015", "Mexico_2015"),
                        stringsAsFactors = F)
write.table(kohn_meta, "kohn.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
# get wallberg meta data
wallberg_ACMO = data.frame(geographic_location = c("Italy", "Austria", 
                                                   "Jordan", "Turkey",
                                                   "South Africa", "Nigeria",
                                                   "Norway", "Sweden", "Spain"),
                           group = c("C", "C",
                                     "O", "O",
                                     "A", "A",
                                     "M", "M", "M"),
                           stringsAsFactors = F)
wallberg_meta <- read.table("../bee_samples_listed/Wallberg_2014_SRA.info",
                            header = T, sep = "\t", stringsAsFactors = F) %>%
  unique(.) %>%
  dplyr::select(c("SRA_Sample", "geo_loc_name", "strain")) %>%
  rename(Bee_ID = SRA_Sample,
         geographic_location = geo_loc_name) %>%
  mutate(source = "Wallberg",
         year = 2014) %>%
  left_join(., wallberg_ACMO, by = "geographic_location")
write.table(wallberg_meta, "Wallberg_2014.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)

# merge all meta data
meta <- bind_rows(all, newCA, newAR, newMX, kohn_meta, wallberg_meta) %>%
  mutate(geographic_location_short = sapply(geographic_location, function(x)
    ifelse(stringr::str_count(x, ",") == 0,
           x,
           trimws(stringr::str_split_fixed(x, ",", n = 2)[2]))))


# write file with all metadata together
write.table(meta, "all.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
