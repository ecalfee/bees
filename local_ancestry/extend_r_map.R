library(dplyr)
library(ggplot2)

# looking at recombination map for new genome HAv3.1:
rmaps <- lapply(c("100", "10000", "1000000"),
                function(x) read.table(paste0("../data/recomb_map/Wallberg_HAv3.1/apis_ALL.rates.bp1.w_headers.csv.windows.",
                                              x,
                                              ".csv"), sep = "\t", header = T))

lapply(rmaps, function(x) mean(x$RATE*100/(3*442757), na.rm = T))
lapply(rmaps, function(x) var(x$RATE, na.rm = T))
lapply(rmaps, function(x) max(x$WIN_NR, na.rm = T))
lapply(rmaps, function(x) hist(x$RATE))
lapply(rmaps, function(x) min(x$RATE, na.rm = T))
par(mfrow=c(3, 1))
lapply(rmaps, function(x) with(x, plot(WIN_NR, RATE)))
par(mfrow=c(1, 1))
with(rmaps[[3]], plot(WIN_NR, RATE))
rmaps[[3]] %>%
  mutate(cum_pos = cumsum(100000*RATE)) %>%
  ggplot(., aes(x = WIN_NR, y = cum_pos, color = CHROM)) +
  geom_line()
rmaps[[1]] %>%
  mutate(cum_pos = cumsum(100000*RATE)) %>%
  ggplot(., aes(x = WIN_NR, y = cum_pos, color = CHROM)) +
  geom_line()
rmaps[[2]] %>%
  mutate(gain = ifelse(is.na(RATE), 0, 1000*RATE)) %>%
  mutate(cum_pos = cumsum(gain)) %>%
  ggplot(., aes(x = WIN_NR, y = cum_pos, color = CHROM)) +
  geom_line() +
  facet_wrap(~CHROM, scales = "free")

rmaps[[1]] %>% # takes a bit to plot
  mutate(gain = ifelse(is.na(RATE), 0, 100*RATE)) %>%
  mutate(cum_pos = cumsum(gain)) %>%
  ggplot(., aes(x = WIN_NR, y = cum_pos, color = CHROM)) +
  geom_line() +
  facet_wrap(~CHROM, scales = "free")

# A few windows at 10kb scale do not have rates:
table(is.na(rmaps[[2]]))
rmaps[[2]][is.na(rmaps[[2]][ , "RATE"]), ]

# extend the recombination map to the end of the chromosome
# by using adjacent recombination bin as the r value:
max_wind_10kb <- rmaps[[2]] %>%
  group_by(CHROM) %>%
  summarise(max_win = max(WIN_NR))
rmap_10kb_ext <- rmaps[[2]] %>%
  left_join(., max_wind_10kb, by = "CHROM") %>%
  mutate(rate_1_prev = c(NA, RATE[1:(nrow(.) - 1)])) %>%
  mutate(rate_1_after = c(RATE[2:(nrow(.))], NA)) %>%
  mutate(rate_ext = ifelse(!is.na(RATE), RATE,
                           ifelse(WIN_NR == 0, 
                                  rate_1_after,
                                  ifelse(WIN_NR == max_win,
                                         rate_1_prev,
                                         NA))))
rmap_10kb_ext[is.na(rmap_10kb_ext$RATE), ]
lapply(rmaps, function(x) mean(x$RATE*100000/(3*442757), na.rm = T)) # convert from rho/kb -> cM/Mb
lapply(rmaps, function(x) median(x$RATE*10^6*100/1000/(3*442757), na.rm = T))
rmap_10kb_ext %>%
  rename(chr = CHROM) %>%
  mutate(start = START - 1) %>% # make zero indexed for bed file
  mutate(end = start + 10000) %>%
  mutate(cM_Mb = rate_ext*100000/(3*442757)) %>%
  mutate(off_end_map = is.na(RATE)) %>%
  dplyr::select(c("chr", "start", "end", "cM_Mb", off_end_map)) %>%
  write.table(.,
              paste0("../data/recomb_map/Wallberg_HAv3.1/map_rates_extended_10kb.bed"),
              sep = "\t",
              quote = F,
              row.names = F)

