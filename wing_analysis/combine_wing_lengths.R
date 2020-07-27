library(dplyr)
library(tidyr)
library(ggplot2)

# load measurements and metadata for reference bee populations:
acm.meta <- read.table("results/reference_samples_scut_car_lig_iber_mell.csv",
                       header = T, sep = ",", stringsAsFactors = F)
acm.measurements <- read.table("results/results_museum_bees_1_per_pop_wing_lengths_1.29.20.csv",
                               header = T, sep = ",", stringsAsFactors = F) %>%
  rename(length_pixels = Length) %>%
  mutate(id = stringr::str_replace(Label, ".bmp", "")) %>%
  mutate(type = ifelse(Label == "Masstab_1107-Sud-Afrika-21.bmp", "ruler", "forewing"))
acm.ruler <- filter(acm.measurements, type == "ruler") %>%
  mutate(pixels_per_cm = length_pixels)
subspecies_codes <- data.frame(subspecies_code = c(1, 2, 6, 9, 21), 
                               subspecies = c("mellifera", "iberiensis", "carnica", "ligustica", "scutellata"),
                               ACM_group = c("M", "M", "C", "C", "A"),
                               stringsAsFactors = F)
acm.wings <- filter(acm.measurements, type == "forewing") %>%
  separate(., col = "id", into = c("pop_id", "bee_id", "subspecies_code")) %>%
  mutate(pixels_per_cm = acm.ruler$pixels_per_cm,
         wing_cm = length_pixels/pixels_per_cm,
         pop_id = as.integer(pop_id),
         bee_id = as.integer(bee_id),
         subspecies_code = as.integer(subspecies_code)) %>%
  left_join(., subspecies_codes, by = "subspecies_code") %>%
  left_join(., acm.meta, by = c("pop_id"="CASEID")) %>%
  filter(!is.na(ACM))


# load unscaled wing measurements
measurements <- read.table("results/results_wings_rulers_1.28.20.csv", 
                           sep = ",", header = T) %>%
  mutate(id = stringr::str_replace(Label, "[.+]JPG", "")) %>%
  #id = stringr::str_replace(id)) %>%
  tidyr::separate(., id, c("wing", "letter", "Bee_ID", "month", "day", "year"), 
                  convert = F, sep = "-") %>%
  tidyr::unite(., "date", c("month", "day", "year"), sep = ".") %>%
  mutate(Bee_ID = stringr::str_replace(Bee_ID, "[_]", ""),
         type = ifelse(X <= 293, "forewing", "ruler")) %>%
  rename(length_pixels = Length) %>%
  filter(., Label != "UNK_RW-G-AR1103-05-30-19.JPG") # label error; it's CA1103 but I already have a good image of this bee wing
# extract ruler/scale measurements
rulers <- measurements %>%
  filter(type == "ruler") %>%
  mutate(length_cm = ifelse(length_pixels < 900, 0.5, 1), # all images have 0.5cm or 1cm length measurement
         pixels_per_cm = length_pixels/length_cm) %>%
  arrange(desc(X)) %>%
  filter(!duplicated(Label)) # take last measurement only -- some rulers were remeasured
rulers %>%
  ggplot(., aes(x = pixels_per_cm)) +
  geom_histogram() # some variation in microscope magnification across measurements means different scales pixels per cm

wings <- measurements %>%
  filter(type == "forewing") %>%
  filter(!duplicated(Label)) %>%
  left_join(., rulers[ , c("Label", "pixels_per_cm")], by = "Label") %>%
  mutate(wing_cm = length_pixels/pixels_per_cm) # scale wing measurements
wings %>% # one bee was imaged twice
  filter(Bee_ID == "AR0605") # same bee diff images, .003cm diff. in measurement error
wings <- filter(wings, Label != "RW-D-AR0605-01-25-20.JPG") # just use the original higher res. image
rulers <- filter(rulers, Label != "RW-D-AR0605-01-25-20.JPG") # duplicated pictures, just use other higher res one

# write data
write.table(acm.wings, "results/ACM_bee_wing_lengths_cm.txt",
            col.names = T, row.names = F, sep = "\t")

write.table(wings[ , c("Bee_ID", "Label", "type", "wing_cm", "length_pixels", "pixels_per_cm", "wing", "date")],
            file = "results/bee_wing_lengths_cm_1.28.20.txt",
            col.names = T, row.names = F, sep = "\t")

rulers %>%
  dplyr::select(Bee_ID, Label, pixels_per_cm) %>%
  write.table(., file = "results/bee_scale_1cm_rulers_1.28.20.txt", col.names = T, row.names = F, sep = "\t")
# note: there are more rulers than wing lengths because some wings were broken for full length but are otherwise ok for vein measurements etc.

# load bee metadata
PREFIX = "combined_sept19"
# bees
pops <- read.table("../bee_samples_listed/byPop/combined_sept19_pops.list",
                   stringsAsFactors = F)$V1
# ancestries
ACM <- c("A", "C", "M")
bees <- do.call(rbind, 
                lapply(pops, function(p) data.frame(Bee_ID = read.table(paste0("../bee_samples_listed/byPop/", p, ".list"),
                                                                        stringsAsFactors = F)$V1, population = p, stringsAsFactors = F)))

# metadata
load("../local_ancestry/results/meta.RData")
# get admixture data
load("../global_ancestry/results/NGSAdmix/ACM_K3_combined_sept19_chr_prunedBy250.rData") 

# combine data
wings.meta <- wings %>%
  left_join(., dplyr::select(meta.ind, -c(date, time)), by = "Bee_ID") %>%
  left_join(., d_admix_ACM_labelled[ , c("Bee_ID", ACM)], by = "Bee_ID")

# make output file of all samples for supplementary information:
# translate names from Oberursel wing data
german_names <- data.frame(LANDNAM = c("Italien", "Tanzania", "England", "Norwegen", "Rhodesien", 
                                       "Frankreich", "Kenya", "Osterreich", "JugoslawienSlow", 
                                       "JugoslawienSerb", "Rumanien", "Ungarn", "Burundi", 
                                       "Rusland", "JugoslawienKroa", "Spanien", "Sud-Afrika"), 
                           location = c("Italy", "Tanzania", "England", "Norway", "Rhodesia (Zimbabwe)", 
                                        "France", "Kenya", "Austria", "Yugoslavia (Slovenia)", 
                                        "Yugoslavia (Serbia)", "Romania", "Hungary", "Burundi", 
                                        "Russia", "Yugoslavia (Croatia)", "Spain", "South Africa"))
sources <- data.frame(source = c("Ramirez", "Sheppard", "Harpur", "Calfee"),
                      publication = c("Published in Cridland et al. 2018, NCBI PRJNA385500",
                                      "Published in Cridland et al. 2017, collected by Sheppard et al., NCBI PRJNA294105", 
                                      "Published in Harpur et al. 2014, NCBI PRJNA216922",
                                      "This study"))


# add notes about bees collected near observed feral colonies
enjambre <- data.frame(Bee_ID = 
                         c("AR1002", 
                           "AR1403", 
                           "AR1511", 
                           "AR2002", 
                           "AR2201",
                           "AR2512",
                           "AR2603",
                           "AR2913"),
                       collection_note =
                         c("at entrance of a feral colony in an old gas tank",
                           "foraging on water within 5m of a feral colony",
                           "on ground within 5m of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a tree cavity"),
                       stringsAsFactors = F)

bees_all <- d_admix_ACM_labelled %>%
  left_join(., meta.ind[ , c("Bee_ID", "date", "time")], by = "Bee_ID") %>%
  dplyr::select(Bee_ID, geographic_location_short, population, date, time, enjambre, lat, long, est_coverage, A, C, M, source) %>%
  left_join(., full_join(dplyr::select(wings, -date), 
                         dplyr::select(rulers, c("Bee_ID", "Label")), 
                         by = "Bee_ID") %>% # get label from ruler if there is no measurement (e.g. broken wing has no length but has an image + ruler)
               mutate(Label = ifelse(is.na(Label.x), as.character(Label.y), as.character(Label.x))) %>% # keep wing image label if exists, if not take ruler label instead
               dplyr::select(., Bee_ID, Label, wing_cm),
            by = "Bee_ID") %>%
  mutate(location = ifelse(source == "Ramirez", "California", 
                           ifelse(geographic_location_short == "Pag Island, Croatia", "Croatia", geographic_location_short))) %>%
  left_join(., sources, by = "source") %>%
  bind_rows(., acm.wings %>%
              mutate(Bee_ID = paste0("Pop", pop_id, "_Ind", bee_id)) %>% 
              rename(population = ACM) %>%
              left_join(., german_names, by = "LANDNAM") %>%
              dplyr::select(Bee_ID, Label, wing_cm, location, population) %>%
              mutate(publication = "Published in Ruttner 1988, Oberursel Collection, Germany")) %>%
  rename(wing_image_file = Label) %>%
  rename(wing_length_cm = wing_cm) %>%
  rename(collection_date = date,
         collection_time = time,
         feral_nest = enjambre) %>%
  left_join(., enjambre, by = "Bee_ID") %>%
  mutate(collection_note = ifelse(is.na(collection_note) & publication == "This study", 
                                  "foraging on vegetation", 
                                  collection_note)) %>%
  rename(., latitude = lat, longitude = long) %>%
  dplyr::select(Bee_ID, location, population, latitude, longitude, collection_date, collection_time, collection_note,
                feral_nest, publication, wing_image_file, wing_length_cm, est_coverage, A, C, M) %>%
  arrange(latitude, population, location) %>%
  arrange(publication == "Published in Ruttner 1988, Oberursel Collection, Germany",
          publication != "This study")
#View(bees_all)

write.table(bees_all, "../../bee_manuscript/files_supp/S1_Table_Sample_information.txt",
            col.names = T, row.names = F, sep = "\t", quote = F) # write out supplementary table with sample information
save(file = "results/wings_bees_all.RData", list = c("bees_all"))
