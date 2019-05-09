library(dplyr)
library(tidyr)
library(ggplot2)
d <- do.call(rbind, 
                lapply(5:13, 
                       function(x) 
                         read.table(paste0("extract_", x, ".txt"),
                   header = T, 
                   stringsAsFactors = F,
                   sep = "\t"))) %>%
  mutate(pop = substr(Bee_ID, 1, 2)) %>%
  arrange(n) %>%
  arrange(extract) %>%
  arrange(pop == "MX") %>%
  arrange(lane) %>%
  mutate(plate_row = rep(rep(LETTERS[1:7], rep(8, 7)), 4)) %>%
  mutate(plate_col = rep(rep(1:8, 7), 4)) %>%
  mutate(P7 = c(rep(paste0(rep(letters[1:7], c(rep(8, 7))), # lanes 2-4 use newer 56 barcodes A-G, 1-8
                         rep(1:8, 7)), 3), 
         paste0(rep(letters[2:8], c(rep(8, 7))), 
                rep(1:8, 7)))) # lane 5 uses older subset from A-H, 1-8, but where h8 doesn't exist
# and some of the other barcodes are lower because of repeated use, 
# so I switched some of them out for a1-8 barcodes.
d[d$lane == 5 & d$P7 %in% c("f6", "g1", "g7", "h5", "h7", "h8"), "P7"] <- paste0("a", c(4:8, 1)) # final

test <- d[1:56, ]
test %>%
  mutate(x = paste(extract, n, sep = "-")) %>%
  select(x, plate_row, plate_col) %>%
  tidyr::spread(., plate_col, x)
wells <- lapply(1:5, function(l)
  filter(d, lane == l) %>%
  mutate(x = paste(extract, n, sep = "-")) %>%
  select(x, plate_row, plate_col) %>%
  tidyr::spread(., plate_col, x))
P7s <- lapply(1:5, function(l)
  filter(d, lane == l) %>%
    select(P7, plate_row, plate_col) %>%
    tidyr::spread(., plate_col, P7))
lapply(2:5, function(l) write.table(wells[[l]],
  paste0("lane_", l, "_plate_layout.txt"),
                        quote = F, row.names = F,
                        sep = "\t"))
lapply(2:5, function(l) write.table(P7s[[l]],
                                    paste0("lane_", l, "_P7_layout.txt"),
                                    quote = F, row.names = F,
                                    sep = "\t"))
# make list for novogene of libraries, sample ID's, and barcodes
barcodes <- read.table("i7_barcodes.txt", stringsAsFactors = F)
colnames(barcodes) <- c("index", "p7_id")
d_novogene <- barcodes %>%
  mutate(P7 = tolower(p7_id)) %>% # make lower vs. upper case letters
  left_join(d, ., by = "P7") %>%
  mutate(., library = paste0("EC_E2_L", lane))
d_novogene %>%
  select(c("library", "Bee_ID", "index", "lane", "extract", "n", "plate_row", "plate_col", "P7")) %>%
  write.table(., "Novogene_libraries2-5_barcodes.txt",
              quote = F, row.names = F)

