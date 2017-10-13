# I will be using the observed crosses in haploid drones from Liu 2015 to create a smoothed recombination map 
# to get estimated genetic positions for physical markers

# Load recombination map from Liu 2015
r <- read.csv("data/recomb_map/Liu_2015/S5_rates.csv")
colnames(r) <- c("linkage_group", "pos_Mb", "n_crossovers_per_drone", "rate")
# there were 43 drones, but the rates per drone must be heavily rounded
unique(43*r$n_crossovers_per_drone)

r$CO_counts <- round(r$rate/min(r[r$rate != 0, "rate"]))
r$CO_counts <- round(r$rate*43/10000) #equivalent. Need to round because of shortened decimal places when file was converted to .csv


sum(r$CO_counts) #2177 not expected 3505 -- are some missing? 
# Yes, emailed authors about the discrepancy and some are omitted because they fall in genome assembly gaps of unknown size.
# these gaps are not included in the genetic length either, but increase uncertainty about local recomb. rate in those regions.
unique(r$linkage_group)
# I could make a map based on the posterior probability based on # of crossovers observed:
# prior on local recombination rate (lambda) centered around mean for the whole genome -- gamma? not sure.
# calculate posterior
genome_length = dim(r)[1] #in .01 Mb units, this is the total length included in Liu 2015 (skips gaps)
mean_r = sum(r$CO_counts)/genome_length
#dpois(x = unique(r$CO_counts), lambda = mean_r)
plot(unique(r$CO_counts), 
     log(dpois(x = unique(r$CO_counts), lambda = mean_r)),
     main = paste("Prob. observations under a fixed genomewide rate ", round(mean_r, digits = 2)/.01, "r/Mb"),
     ylab = "log probability",
     xlab = "counted recombination events in .01 Mb intervals")
# how different are the LD based map and the posterior distribution based on the observed meioses?