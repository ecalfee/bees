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
  mutate(P7 = rep(paste0(rep(letters[2:8], c(rep(8, 7))), 
                         rep(1:8, 7)), 4))
# h8 doesn't exist and some of the other barcodes are lower because of repeated use, 
# so I switch some of them out for a1-8 barcodes.
d[d$P7 == "h8" | (d$lane == 5 & P7 %in% c("f6", "g1", "g7", "h5", "h7")), "P7"] <- paste0("a", c(1:8, 1))
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


