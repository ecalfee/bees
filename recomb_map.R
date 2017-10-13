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
unique(r$linkage_group)
# Emailed authors about the problem.