library(dplyr)
library(ggplot2)
file_in <- "coverage/backup_pass1_plus_kohn_and_wallberg_random_100_regions.coverage"
coverage <- read.table(file_in, stringsAsFactors = F)
colnames(coverage) <- c("scaffold", "start", "end", "x_coverage", "n_bp", "region_length", "perc_bp")
by_region <- coverage %>%
  group_by(paste(scaffold, start, end, sep = "_")) %>%
  summarize(mean_coverage = sum(x_coverage*perc_bp)) 
summary(by_region$mean_coverage)
hist(by_region$mean_coverage)
coverage %>%
  filter(x_coverage <= 900) %>%
  summarize(count_low = sum(n_bp))/sum(coverage$n_bp)
coverage %>%
  filter(x_coverage >= 4000)%>%
  summarize(count_high = sum(n_bp))/sum(coverage$n_bp)
summary(coverage$x_coverage)
