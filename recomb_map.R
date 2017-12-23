#library(rethinking)
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
# calculate posterior -- moved models down to end of file.
genome_length = dim(r)[1] #in .01 Mb units, this is the total length included in Liu 2015 (skips gaps)
mean_r = sum(r$CO_counts)/genome_length
r$id = 1:length(r$CO_counts) # add an id column for the segment of the genome of interest
r$lg = as.numeric(as.factor(r$linkage_group))

# THe modeling above (now at end of file) doesn't make a lot of sense -- instead I will try to fit a monotone increasing polynomial
# as a smoothing technique
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
       degree = 3, plot.it = T) # will need to plot on my own just final fit (not all)

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
starts = sapply(start$V1, function(x) as.numeric(strsplit(x, split = "[:>]")[[1]][1]))
LG = data.frame(id = last$V1[c(T,F,F)][1:340], tail = last$V1[c(F,F,T)], 
                start_line = starts[1:(length(starts)-1)], # ignore last
                next_start_line = starts[-1], # ignore first
                stringsAsFactors = F)

# extract scaffold ids
LG$id = substr(LG$id, 2, nchar(LG$id))
LG$lg = as.integer(sapply(LG$id, function(x) substr(strsplit(x, split = "[>.]")[[1]][1], 6, 100000))) # e.g. Group1
LG$scaffold = as.numeric(sapply(LG$id, function(x) strsplit(x, split = "[>.]")[[1]][2])) # e.g. 3 for scaffold

# there are 70 characters = bp in all but the last line
max(sapply(LG$tail, nchar)) # maximum of last lines is also 70
nchar("AACGAGCACATGACAAGTAGTACATGACTTCCCATCCTTGGATCCAGACAAGATCCAAAGATGAGAAACC") # line taken from file

# calculate length of each scaffold (in bp)
LG$len = nchar(LG$tail) + 70*(LG$next_start_line - LG$start_line - 2)

# finally, I calculate the start position for each scaffold (relative to the whole linkage group)
LG$pos = NA # placeholder
LG[LG$scaffold == 1, "pos"] = 1 # first scaffold always starts at position 1 bp

# then arrange scaffolds within chromosomes in numerical order and get starting position (because 11 originally follows 10 follows 1)
LG = arrange(LG, scaffold) %>% 
  arrange(., lg)
for (lg in unique(LG$lg)){
  for (i in 2:length(unique(LG$scaffold[LG$lg == lg]))){
    LG[(LG$lg == lg & LG$scaffold == i), "pos"] = 
      LG[(LG$lg == lg & LG$scaffold == (i-1)), "pos"] + 
      LG[(LG$lg == lg & LG$scaffold == (i-1)), "len"] +
      50000 # see below but 50000 bp is the asummed length in Amel_4.5 for any unspanned gap
  }
}
# check that it worked visually

plot(LG[LG$lg == "Group1", "pos"], 
     LG[LG$lg == "Group1", "scaffold"], main = "plot chr 1 scaffold positions",
     ylab = "Group 1 scaffold #",
     xlab = "starting position of chromosome 1")
scaff_lengths = tapply(LG$len, LG$lg, sum)
sum(tapply(LG$len, LG$lg, sum)) # length of concatenated chromosomes matches positions, good
sapply(unique(LG$lg), function(x) 
  (max(LG[LG$lg == x, "pos"]) + LG$len[which(LG$lg == x & LG$pos == max(LG[LG$lg == x, "pos"]))] - 1)/1000000)
# So it's internally consistent and matches Amel_4.5 (see below),and with addition of 50000 bp gaps matches Liu 2015 
# and 'total length' in NCBI. Note that 'ungapped length' takes out medium-sized runs of NNNNs that are spanned gaps.
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000002195.4/#/st

# I have also confirmed using command line that counting all ACTGN's (minus all G's in Group1 header etc)
# equals the total genome length I calculate above.
# I do not understand how the 'total length' is determined or how the gaps are distributed inbetween scaffolds
# As a rough approximation for now, I will just treat the difference between true length and length I calculate 
# divided by total unspanned gaps as the expected unspanned gap length for each unspanned gap on a chromosome
# gunzip < Amel_4.5_scaffolds.fa.bgz | head -2906642 | tr -cd ACTGN | wc -c # counted all ACTGN's before first GroupUn (unmapped)
# subtract all G's from Group1 etc. not actual nucleotides
# gunzip < Amel_4.5_scaffolds.fa.bgz | head -2906642 | grep "^>Group" | wc -l
203429952-341 # [1] 203429611
# can see gaps placed here too:
# https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000002195.4

genome_meta = read.csv("data/honeybee_genome/ncbi_genome4.5_scaffold_stats.csv")[-17,] # get rid of last row (unmapped scaffolds)
gap_length = (genome_meta$total_length - scaff_lengths)/genome_meta$unspanned.gaps
# ok so what I found out is that all the unspanned gaps are put in as 50000 bp gaps as a proxy, which I can do too.

rownames(LG) = NULL
  #LG$tail = LG$start_line = LG$tail_n = NULL
write.csv(LG[,c("id", "lg", "scaffold", "pos", "len")], file = "data/honeybee_genome/Amel_4.5_scaffold_start_positions.csv", quote = F)

s <- read.csv("data/honeybee_genome/Amel_4.5_scaffold_start_positions.csv", stringsAsFactors = F)
# so now I just need to add the scaffold-relative position to their start position to get an absolute chromosome position per marker


# create ELAI input - positions file
# start with just 1 scaffold, e.g. Group1.1
raw_bp = read.csv("data/geno_AC/geno_AC_Group1.1.mafs.gz", sep = "\t", stringsAsFactors = F, header = T)
raw_bp$snp_id = paste0(raw_bp$chromo, "_", raw_bp$pos)
length(raw_bp$knownEM)/max(raw_bp$position)
table(diff(raw_bp$position))
mean(diff(raw_bp$position))
median(diff(raw_bp$position))
filter_bp = raw_bp[raw_bp$knownEM > .02,]
thres_bp = raw_bp[raw_bp$pK.EM < .000001,]
dim(filter_bp)
dim(thres_bp)
length(filter_bp$knownEM)/max(filter_bp$position)
table(diff(filter_bp$position))
length(raw_bp$knownEM)



# which column goes with which bee/pop?
meta = read.csv("Bee_bam_include_geno_AC_12_metadata.txt", stringsAsFactors = F, header = T)
meta$pop = ifelse(is.na(meta$year), meta$pop, paste0(meta$pop, "_", meta$year)) 
meta$pop = ifelse(meta$pop == "S", "A", meta$pop) # set S and A equivalent (both A ancestry)

# genotype data
raw_geno = read.csv("data/geno_AC/geno_AC_Group1.1.geno.gz", sep = "\t", stringsAsFactors = F, header = F)
raw_geno$V1 <- NULL
raw_geno$V2 <- NULL
colnames(raw_geno) = meta$id

# make ELAI files for one short scaffold Group1.1

# function below formats a genotype input file for ELAI based on the full set of raw genotypes and snp names
write_elai_geno_input = function(columns, snp_ids, raw_geno, output_file){ #columns = subset of columns from genotype file to use
  write.table(length(columns), output_file, sep = ",", quote = F, append = F, row.names = F, col.names = F)
  write.table(length(snp_ids), output_file, sep = ",", quote = F, append = T, row.names = F, col.names = F)
  write.table(data.frame(IND = snp_ids, raw_geno[,columns]), output_file, sep = ",", 
              quote = F, append = T, row.names = F, col.names = T)
}
# make input genotype files for ancestral and 2 admixed populations
for (p in c("A", "M", "C", "Placerita_2014", "Riverside_2014")){
  write_elai_geno_input(columns = which(meta$pop == p), snp_ids = raw_bp$snp_id, raw_geno = raw_geno, 
                        output_file = paste0("data/ELAI_input/Group1.1/", p, ".txt"))
}
# combined Plac_Riv_2014
write_elai_geno_input(columns = which(meta$pop == "Placerita_2014" | meta$pop == "Riverside_2014"), snp_ids = raw_bp$snp_id, raw_geno = raw_geno, 
                      output_file = paste0("data/ELAI_input/Group1.1/", "Plac_Riv_2014", ".txt"))

# function below formats a position input file for ELAI based on a set of snps
# format: (no header) snp_id, position on chromosome (in bp), chromosome name
# ELAI assumption is 1,000,000 bp (humans) = 1 cM, or 1Mbp per cM. So I will 1) put scaffolds together to make chromosomes in bp
# then 2) convert bp to cM for honeybee genome markers and finally 
# 3) multiply my cM by 1,000,000 to get a position on a human-bp scale
# Initially I will use just the chromosome avg recombination rate for 2) and skip 1) for just 1 scaffold
avg_r_per_chr = read.csv(file = "data/recomb_map/Liu_2015/avg_r_per_chr_Liu_2015.csv", header = T, stringsAsFactors = F)
avg_r = avg_r_per_chr[avg_r_per_chr$chr_n == 1, "avg_r_cMperMb"] # cM/Mbp
raw_bp$ELAI_pos = round(avg_r*raw_bp$position, digits = 0) # round to nearest whole bp
write.table(raw_bp[ , c("snp_id", "ELAI_pos", "chromo")], "data/ELAI_input/Group1.1/pos.txt", sep = ",", 
            quote = F, append = F, row.names = F, col.names = F)

# make all input files for ELAI
for (chr in 1:16){
  raw_geno = read.csv(paste0("data/geno_AC/All/All_AC_", chr, ".geno"), sep = "\t", stringsAsFactors = F, header = F)
  snp_id = paste0(raw_geno$V1, "_", raw_geno$V2)
  snp_raw_pos = raw_geno$V2 
  snp_lg = raw_geno$V1
  raw_geno$V1 <- NULL
  raw_geno$V2 <- NULL
  colnames(raw_geno) = meta$id
  # write genotype files
  for (p in c("A", "M", "C")){
    write_elai_geno_input(columns = which(meta$pop == p), snp_ids = snp_id, raw_geno = raw_geno, 
                          output_file = paste0("data/ELAI_input/chr", chr, "/", p, ".txt"))
  }
  # combined Plac_Riv_2014
  write_elai_geno_input(columns = which(meta$pop == "Placerita_2014" | meta$pop == "Riverside_2014"), 
                        snp_ids = snp_id, raw_geno = raw_geno, 
                        output_file = paste0("data/ELAI_input/chr", chr, "/", "Plac_Riv_2014", ".txt"))
  
  # write snp position files
  avg_r = avg_r_per_chr[avg_r_per_chr$chr_n == i, "avg_r_cMperMb"] # cM/Mbp
  
  snp_new_pos = sapply(1:length(snp_raw_pos), function(p) 
    round(avg_r*(snp_raw_pos[p] + s[s$id == snp_lg[p], "pos"]), 
          # add starting position of the scaffold to raw relative position & multiply by chromosome recomb. rate 
          digits = 0)) # then round to nearest whole bp
  write.table(data.frame(snp_id, snp_new_pos, rep(chr, length(snp_id))),
              paste0("data/ELAI_input/chr", chr, "/pos.txt"), sep = ",", 
              quote = F, append = F, row.names = F, col.names = F)
  
}



# visualizing elai test results
elai_test_out = read.csv("data/ELAI_input/Group1.1/output/test_Group1.1.ps21.txt",
                           header = F, sep = " ")
elai_test_out2 = elai_test_out[,!is.na(elai_test_out[1,])] # out2 excludes MAF <.1 and in the future I'll use .05 (saves time on less useful snps)
elai_test_out_snps = read.csv("data/ELAI_input/Group1.1/output/test_Group1.1.snpinfo.txt",
                              header = T, sep = "\t")
elai_test = list(A = t(elai_test_out2[, c(T, F, F)]),
                       C = t(elai_test_out2[, c(F, T, F)]),
                       M = t(elai_test_out2[, c(F, F, T)]),
                 snp = elai_test_out_snps$rs,     
                 pos = elai_test_out_snps$pos)
par(mfrow=c(1,1))
for(i in 1:14){
  plot(elai_test$pos, elai_test$A[,i], type = "p", col = "blue")
}





# I have way too high a density of SNPs! (ok for ELAI but not for other software with assumptions about no within-ancestry LD)

# Calculate LD within ancestral population 
frq1 = read.csv("data/geno_12/geno_12_Group1.1.geno.gz", sep = "\t", stringsAsFactors = F, header = F)
frq1[, 1:4] <- NULL
frq1[,c(T,F)] <- NULL # get rid of columns not related to ind. bees

# The modeling at the end doesn't make a lot of sense -- instead I will try to fit a monotone increasing polynomial
# as a smoothing technique
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
# visualize with 9 degrees of freedom
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
# visualize with 5 degrees of freedom
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

# visualize error with 1 degree of freedom (linear - mean is the constant recomb. rate)
#dpois(x = unique(r$CO_counts), lambda = mean_r)
plot(unique(r$CO_counts), 
     log(dpois(x = unique(r$CO_counts), lambda = mean_r)),
     main = paste("Prob. observations under a fixed genomewide rate ", round(mean_r, digits = 2)/.01, "r/Mb"),
     ylab = "log probability",
     xlab = "counted recombination events in .01 Mb intervals")
# how different are the LD based map and the posterior distribution based on the observed meioses?



# for each chromosome, load data into r to manipulate, make elai input files and add recomb. distances
chr = 15
g <- read.csv(paste0("data/geno_AC/All/All_AC_", chr, ".geno"), sep = "\t", stringsAsFactors = F, header = F)

colnames(g) = c("scaffold", "rel_pos", meta$id, "other")
g$other = NULL
g$IND = paste0(g$scaffold, "_", g$rel_pos)

g$start_pos = sapply(g$scaffold, function(i) s[s$id == i, "pos"])
g$abs_pos = g$start_pos + g$rel_pos - 1 # absolute bp position from start of chromosome
# get recombination data for that chromosome
d = data.frame(pos = pos_byLG[[chr]], r = r_byLG[[chr]])
# fit the recombination data to a polynomial degree 5
m5 = monpol(r ~ pos, data = d, 
            ptype= "SOS", monotone = "increasing",
            degree = 5, plot.it = F) 
g$abs_r = sapply(g$abs_pos, function(x) 
  coef(m5)%*%c(1, x, x^2, x^3, x^4, x^5))











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
    warmup = 1000, iter = 2000, chains = 4, cores = 4, start = list(mu = 0.1) #, k = 0.1, theta = 1)
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
  
  # an alternative would be just to use simple prior conjugacy between the gamma and poisson
  # r | lambda ~ pois(lambda*t) where t is the observed meioses and lambda is the recombination rate for an interval = window size
  # the prior for lambda is based on total number of recombinations previously observed in total # of genomic windows*meioses
  # I would have to choose how strong to make the prior based on how much genomic evidence there was
  # lambda ~ gamma(shape, rate) # note I use rate not scale parameterization
  # then the posterior after x observed recombinations in w windows*meioses to observe them =
  # lambda ~ gamma(shape + x, rate + w)
  # unspanned gaps would have no observations and therefore just reflect the prior
  # If I just want the expectation out of this posterior, however,
  # I think this is the equivalent of just taking a weighted average p*(prior mean) + (1-p)*(observed mean)
  # where p is how much you want to rely on prior beliefs .. so I can justify an average basically as the best point approximation
}
