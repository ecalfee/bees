# for bootstrap, divide genome into 0.2cM windows (with non-uniform recomb rates within windows)
# use 10kb recombination map to start:
rmap <- read.table("../data/recomb_map/Wallberg_HAv3.1/map_rates_extended_10kb.bed", header = T, sep = "\t")
# add cM positions to rmap:
rmap <- rmap %>%
  group_by(chr) %>%
  mutate(rpos_end = cumsum((end-start)*cM_Mb/10^6))
# interpolate 0.2cM breakpoints:
rmap2 <- do.call(rbind, # spans all chromosomes at 0.2cM windows 
                 lapply(unique(rmap$chr), function(chromosome)
                   data.frame(approx(x = rmap$rpos_end[rmap$chr == chromosome], 
                                     y = rmap$end[rmap$chr == chromosome], 
                                     xout = seq(0, max(rmap$rpos_end[rmap$chr == chromosome]), by = 0.2), method = "linear")) %>%
                     mutate(., y = ifelse(x == 0, 0, y)) %>%
                     mutate(., chr = chromosome) %>%
                     mutate(., rpos_start = x, start = floor(y), 
                            end = c(sapply(2:nrow(.), function(i) start[i]), 
                                    max(rmap$end[rmap$chr == chromosome]))))) %>%
  mutate(., window_0.2cM = 1:nrow(.)) %>%
  dplyr::select(., chr, start, end, window_0.2cM, rpos_start)
write.table(rmap2, 
            "../data/recomb_map/Wallberg_HAv3.1/map_10kb_in_0.2cM_windows.bed", 
            quote = F,
            col.names = T, row.names = F, sep = "\t")
