library(rethinking)
library(MonoPoly)
# I will be using the observed crosses in haploid drones from Liu 2015 to create a smoothed recombination map 
# to get estimated genetic positions for physical markers

# Load recombination map from Liu 2015
r <- read.csv("data/recomb_map/Liu_2015/S5_rates.csv", stringsAsFactors = FALSE)
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
# prior on local recombination rate (lambda) centered around mean for the whole genome -- gamma? not sure; would estimate both
# calculate posterior -- moved models down to end of file.
genome_length = dim(r)[1] #in .01 Mb units, this is the total length included in Liu 2015 (skips gaps)
mean_r = sum(r$CO_counts)/genome_length
r$id = 1:length(r$CO_counts) # add an id column for the segment of the genome of interest
r$lg = as.numeric(as.factor(r$linkage_group))

# THe modeling above (now at end of file) doesn't make a lot of sense -- instead I will try to fit a monotone increasing polynomial
# as a smoothing technique
r1 = r[r$linkage_group == "LG1",][1:10,]
r1$posMb = sapply(strsplit(r1$pos_Mb, split = "~"), function(x) as.numeric(x[2])) # get physical Mb position
r1$posM = cumsum(r1$n_crossovers_per_drone)
pos_byLG = lapply(unique(r$linkage_group), # for each linkage group (of 16)
             function(g) sapply(strsplit(r[r$linkage_group == g, "pos_Mb"], split = "~"), # get the end of each interval
                               function(x) as.numeric(x[2]))) # as a number
r_byLG = lapply(unique(r$linkage_group),
                function(g) cumsum(r[r$linkage_group == g, "n_crossovers_per_drone"]))
par(mfrow=c(2,2))
lapply(1:16, function(i) plot(pos_byLG[[i]], r_byLG[[i]], 
                              xlab = "position (omit gaps)",
                              ylab = "r (Morgans)",
                              main = unique(r$linkage_group)[i]))
#approxfun
#monopoly -smooths over larger windows
d = data.frame(list(pos_byLG = pos_byLG[[15]], r_byLG = r_byLG[[15]]))
par(mfrow=c(1,1))
m1 = monpol(r_byLG ~ pos_byLG, data = d,
       degree = 9, plot.it = T) # will need to plot on my own just final fit (not all)

#dpois(x = unique(r$CO_counts), lambda = mean_r)
plot(unique(r$CO_counts), 
     log(dpois(x = unique(r$CO_counts), lambda = mean_r)),
     main = paste("Prob. observations under a fixed genomewide rate ", round(mean_r, digits = 2)/.01, "r/Mb"),
     ylab = "log probability",
     xlab = "counted recombination events in .01 Mb intervals")
# how different are the LD based map and the posterior distribution based on the observed meioses?

# because this genetic map concatenates scaffolds within linkage groups together, and I have SNP position relative
# to the scaffold not whole LG, I need to know how long each scaffold is to find where they connect
start = read.csv(file = "data/honeybee_genome/start_line_n.txt", header = F, stringsAsFactors = F)
last = read.csv(file = "data/honeybee_genome/last_line_prev.txt", header = F, stringsAsFactors = F)
LG = data.frame(id = last$V1[c(T,F,F)], tail = c("", last$V1[c(F,F,T)]), 
                start_line = sapply(start$V1, function(x) as.numeric(strsplit(x, split = "[:>]")[[1]][1])),
                stringsAsFactors = F)
LG$tail_n = nchar(LG$tail)
# there are 70 characters in a normal line
max(LG$tail_n) # maximum of last lines is also 70
nchar("AACGAGCACATGACAAGTAGTACATGACTTCCCATCCTTGGATCCAGACAAGATCCAAAGATGAGAAACC") # line taken from file
LG$lg = sapply(LG$id, function(x) strsplit(x, split = "[>.]")[[1]][2]) # e.g. Group1
LG$scaffold = as.numeric(sapply(LG$id, function(x) strsplit(x, split = "[>.]")[[1]][3])) # e.g. 3 for scaffold
# finally, I calculate the start position for each scaffold (relative to the whole linkage group)
LG$pos = NA # placeholder
LG[LG$scaffold ==1, "pos"][1] = 1 # first scaffold always starts at position 1
for (g in unique(LG$lg)){
  n = unique(LG[LG$lg == g, "scaffold"])
  for (i in 2:length(n)){
    LG[(LG$lg==g & LG$scaffold==n[i]), "pos"] = LG[(LG$lg==g & LG$scaffold==n[i]), "tail_n"] + 
      70*(LG[(LG$lg==g & LG$scaffold==n[i]), "start_line"] - 
            LG[(LG$lg==g & LG$scaffold==n[i-1]), "start_line"] - 2) + 
      LG[(LG$lg==g & LG$scaffold==n[i-1]), "pos"] # length of prior segment + start pos prior segment is new position (start of current scaffold)
   # minus 2 because 1 line is taken up by the scaffold label and the other is a partial line accounted for in tail_n
  }
}
LG$id = sapply(LG$id, function(x) substr(x, 2, nchar(x))) # get rid of leading ">"


# also sum up all the lengths to make sure I get the correct total length!!!!
rownames(LG) = LG$tail = LG$start_line = LG$tail_n = NULL
write.csv(LG, file = "data/recomb_map/Liu_2015/scaffold_start_positions.csv", quote = F)








# Some models I tried but didn't use:
dont_run = function(x){
# mean of gamma = shape*scale (scale = theta = inverse of beta)
m0.r <- map2stan(
  alist(
    CO_counts ~ dpois(mu), # observed crossovers are poisson with mean = mean rate for whole genome
    mu ~ dgamma(shape = .01, scale = 1) # estimate a mean rate
    # using a hyperprior to estimate 1 mu doesn't make any sense
    #k ~ dgamma(shape = .01, sd = 1), # 1 is close to the mean observation; too close to zero to put a normal prior
    #theta ~ dexp(1)
  ),
  data = r,
  warmup = 1000, iter = 2000, chains = 4, cores = 4, start = list(mu = 0.1, k = 0.1, theta = 1)
)
m0LG.r <- map2stan( # can I first get estimates by linkage group
  alist(
    CO_counts ~ dpois(lambda = lambda), # observed crossovers are poisson with mean = mean rate for whole genome
    lambda <- mu[lg],
    mu[lg] ~ dgamma(shape = k, scale = theta), # estimate a mean rate for each linkage group, with hyperprior over chromosomes
    k ~ dgamma(shape = 0.1, scale = 1), # 1 is close to the mean observation; too close to zero to put a normal prior
    theta ~ dexp(1)
  ),
  data = r,
  warmup = 1000, iter = 2000, chains = 4, cores = 4, start = list(mu = rep(1, 16), k = 0.1, theta = 1)
)
post0LG.r <- extract.samples( m0LG.r )
par(mfrow=c(1,1))
dens(post0LG.r$k)
dens(post0LG.r$theta)
apply(post0LG.r$mu, 2, dens)

r_short = as.list(r[1:100,])
# I thought this is what I wanted to do but it doesn't make a lot of sense --- and would save very large # of chains (21970 of them)
m1.r <- map2stan(
  alist(
    CO_counts ~ dpois(lambda = lambda), # observed crossovers are poisson with mean = mean rate for that local region
    lambda <- mu[id],
    mu[id] ~ dgamma(shape = k, scale = theta), ## estimate a mean rate per segment using a hyperprior on the mean # COs
    k ~ dgamma(shape = 0.1, scale = 1), # 1 is close to the mean observation; too close to zero to put a normal prior
    theta ~ dexp(1) 
  ),
  data = r_short,
  warmup = 10, iter = 10, chains = 4, cores = 4, start = list(mu = rep(1, 21970), k = 0.1, theta = 1)
)
}
