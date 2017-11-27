library(rethinking)
library(MonoPoly)
library(dplyr)
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
# calculate posterior
genome_length = dim(r)[1] #in .01 Mb units, this is the total length included in Liu 2015 (skips gaps)
mean_r = sum(r$CO_counts)/genome_length
r$id = 1:length(r$CO_counts) # add an id column for the segment of the genome of interest
r$lg = as.numeric(as.factor(r$linkage_group))

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


# THe modeling above doesn't make a lot of sense -- instead I will try to fit a monotone increasing polynomial
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
# 9 degrees of freedom appears to capture all major variation across the chromosome
for (i in 1:16){
d = data.frame(list(pos = pos_byLG[[i]], r = r_byLG[[i]]))
par(mfrow=c(1,1))
m1 = monpol(r ~ pos, data = d, 
            ptype= "SOS", monotone = "increasing",
       degree = 9, plot.it = T) 
#plot(seq(-1,1,by=.01), fitted(m1))
# will need to plot on my own just final fit (not all)

dist = seq(0, max(d$pos), by = .1)
plot(x = dist, 
    y = sapply(dist, function(x) 
      coef(m1) %*% c(1, x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9)),
    xlab = "distance (Mbp)",
    ylab = "distance (Morgans)",
    main = paste("chr", i, "fit"))
plot(d$r ~ d$pos, main = paste("chr", i, "raw"))
}

for (i in 1:16){ # only 5 instead of 9 degree polynomial
  d = data.frame(list(pos = pos_byLG[[i]], r = r_byLG[[i]]))
  par(mfrow=c(1,1))
  m1 = monpol(r ~ pos, data = d, 
              ptype= "SOS", monotone = "increasing",
              degree = 5, plot.it = T) 
  #plot(seq(-1,1,by=.01), fitted(m1))
  # will need to plot on my own just final fit (not all)
  
  dist = seq(0, max(d$pos), by = .1)
  plot(x = dist, 
       y = sapply(dist, function(x) 
         coef(m1)%*%c(1, x, x^2, x^3, x^4, x^5)),
       xlab = "distance (Mbp)",
       ylab = "distance (Morgans)",
       main = paste("chr", i, "fit"))
  plot(d$r ~ d$pos, main = paste("chr", i, "raw"))
}



#dpois(x = unique(r$CO_counts), lambda = mean_r)
plot(unique(r$CO_counts), 
     log(dpois(x = unique(r$CO_counts), lambda = mean_r)),
     main = paste("Prob. observations under a fixed genomewide rate ", round(mean_r, digits = 2)/.01, "r/Mb"),
     ylab = "log probability",
     xlab = "counted recombination events in .01 Mb intervals")
# how different are the LD based map and the posterior distribution based on the observed meioses?
# simple test case with 2 scaffolds, 1 linkage group and each scaffold is 100bp long
test1 = data.frame(LG=1, rel_pos = c(12, 16, 5), scaffold_start = c(1, 1, 101),
                   scaffold=c("Group1.1", "Group1.1", "Group1.2"))
test2 = test1 %>%
  mutate(abs_pos = rel_pos + scaffold_start - 1)


#get_fit_pos(data_pos = pos_byLG[[1]], data_r = r_byLG[[1]]),
#abs_pos = test1$abs_pos, degree = 5)
get_fit_pos = function(data_pos, data_r, abs_pos, degree = 9){
  # where degree is the degrees for the polynomial to fit
  # could be chosen by dividing data in half and finding ML of observations
  # using a training and test set
  d = data.frame(pos = data_pos, r = data_r)
  m1 = monpol(r ~ pos, data = d, 
              ptype= "SOS", monotone = "increasing",
              degree = 5, plot.it = F) 
 abs_r_pos = sapply(abs_pos, function(x) 
         coef(m1)%*%c(1, x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9))
  return(abs_r_pos) # return genetic distance from beginning of linkage group
}
# steps to get input for HMM Corbett-Detig:
# 1) Start with bp positions from beginning of scaffold, e.g. Group1.6
# 2) Using scaffold start positions, translate to bp positions 
# from beginning of chromosome, all scaffolds stuck together (no gaps)
# 3) Use monopoly interpolation to get cM position from beginning of chromosome
# 4) Subtract to get cM distance between markers.
# 5) Decide what to do with markers that span a scaffold junction
# 5a) nothing - ignore gap has any recombination
# 5b) treat like start of independent chromosome (extreme independence)
# 5c) assume gaps are all = length and add in extra recomb. distance constant to reflect this
# I can check sensitivity to these choices a-c and also get a sense of effect by just calculating what constant c would be under chromosome mean recomb. rate
