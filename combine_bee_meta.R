library(dplyr)
library(tidyr)
# this script combines metadata for already sequenced ca bees (source: Ramirez),
# Kenya bees (source Sheppard NCBI) and A/C/M/Cerana bees (source: Harpur NCBI)

# harpur A C M group bees & one Cerana reference ind.
harpur <- read.table("data/Harpur_2014_NCBI/Harpur_SraRunTable.txt",
                     sep = "\t", stringsAsFactors = F,
                     header = T)
harpur_incl <- harpur %>%
  select(Run_s, geographic_location_s, strain_s) %>%
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
kenya <- read.table("data/Kenya_Sheppard_NCBI/Kenya_SraRunTable.txt",
                    sep = "\t", stringsAsFactors = F,
                    header = T)
kenya_incl <- kenya %>%
  select(Run_s, Organism_s, collection_date_s) %>%
  mutate(geographic_location = "Kenya") %>%
  rename(strain = Organism_s) %>%
  mutate(population = "A") %>%
  rename(Bee_ID = Run_s) %>%
  mutate(year = as.integer(collection_date_s)) %>%
  select(-collection_date_s) %>%
  mutate(source = "Sheppard") %>%
  mutate(group = "A")

# CA bees from Ramirez lab
ca_bees <- read.table("bee_filter.txt",
                      header = F, stringsAsFactors = F) %>%
  .[,1:2]
colnames(ca_bees) <- c("population", "ID")
ca_list <- read.table("data/CA_Bee/samples.list", header = F, stringsAsFactors = F)$V1 # has 2 fewer because Davis1968 didn't transfer and ap50 is duplicated below
ca_bees_incl <- ca_bees %>%
  dplyr::filter(!(ca_bees$population %in% c("M", "A", "C", "Cerana"))) %>%
  separate(., population, c("geographic_location", "year"), remove = F) %>%
  mutate(strain = "unknown") %>%
  mutate(year = as.integer(year)) %>%
  mutate(source = "Ramirez") %>%
  rename(Bee_ID = ID) %>%
  mutate(group = ifelse(geographic_location %in% c("Berkeley", "Stanislaus", "Humboldt", "Stebbins", "Davis", "Domestic"), "N_CA", "S_CA"))
# note: ap50 is a duplicated record -- I will need to fix this (later)!!

# metadata for all reference panel bees from prior studies
# combine and sort into pop. groups before printing to file
all <- bind_rows(harpur_incl, kenya_incl, ca_bees_incl) %>%
.[order(.$population), ] %>%
#.[order(.$year), ] %>%
.[order(.$strain), ]
write.table(all, "bee_samples_listed/harpur_kenya_ca_bees.meta",
            row.names = F, col.names = T, sep = "\t", quote = F)


# make a list of the reference individuals + ca_bees post 1994 (CA introduction of Afr. honeybees)
post_1994 <- all[all$population %in% c("A", "C", "M") |
                   (all$strain == "unknown" & all$year >= 1994), ]
write.table(post_1994, "bee_samples_listed/post_1994.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
write.table(post_1994$Bee_ID, "bee_samples_listed/post_1994.list",
            row.names = F, col.names = F,
            quote = F)

# add in CA and AR bees newly collected:
newCA <- read.table("maps/CA_bee_lat_long_4_google_maps.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(year = 2018) %>%
  mutate(strain = "unknown") %>%
  mutate(source = "Calfee") %>%
  mutate(population = substr(Bee_ID, 1, 4)) %>%
  mutate(geographic_location = "California") %>%
  mutate(group = "CA_2018")
newAR <- read.table("maps/Arg_bee_lat_long_4_google_maps.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(year = 2018) %>%
  mutate(strain = "unknown") %>%
  mutate(source = "Calfee") %>%
  mutate(population = substr(Bee_ID, 1, 4)) %>%
  mutate(geographic_location = "Argentina") %>%
  mutate(group = "AR_2018")

w_new_bees <- bind_rows(all, newCA, newAR)
w_new_bees <- w_new_bees[!duplicated(w_new_bees),] # temporary patch to remove duplicated ap50 entry
write.table(w_new_bees, "bee_samples_listed/all.meta",
            row.names = F, col.names = T, sep = "\t",
            quote = F)
