#!/usr/bin/env Rscript

# this script takes in an input file for ancestry_hmm, with sample genotypes,
# and simulates binomially sampled reads to ?x coverage across the genome
# using a user-specified random seed 
# and either "fixed" or "poisson" or "negbinom" to indicate whether coverage should be fixed (even across every bp)
# or sampled from a poisson distribution or even more dispersed negative binomial dist (with set variance-to-mean ratio > 1)
# to run: Rscript sim_sample_x_coverage.R ../data/TEST/head_ACM_Riv2014.counts 4 1353323 fixed ../data/TEST/head_ACM_Riv2014.counts

# to test
#args <- c("../data/TEST/head_ACM_Riv2014.counts", "4", "1353323", "fixed", "../data/TEST/head_ACM_Riv2014.counts")
#args <- c("../data/TEST/head_ACM_Riv2014.counts", "4", "1353323", "poisson", "../data/TEST/head_ACM_Riv2014.counts")
#args <- c("../data/TEST/head_ACM_Riv2014.counts", "4", "1353323", "negbinom", "../data/TEST/head_ACM_Riv2014.counts")
# user input:
args <- commandArgs(trailingOnly=TRUE)
INPUT_FILE = args[1]
X = as.numeric(args[2]) # desired coverage to simulate
SEED = as.integer(args[3])
set.seed(SEED) # set random seed
MODEL = args[4]
if (!(MODEL %in% c("poisson", "fixed", "negbinom"))){ # produce error if incorrect input
  stop("4th argument must be poisson or fixed or negbinom")
}
VMR = 3 # used in negative binomial, this sets the variance to mean ratio at 3 
# 3 comes from ReadDepth paper w/ Illumina sequenced human genomes by Christopher Miller 2011: 
# https://doi.org/10.1371/journal.pone.0016327

OBS_ERROR = 0.01 # observed sequencing error rate (1% is standard for illumina)
# (e.g. poisson has VMR = 1; fixed has no variance or VMR = 0
OUTPUT = args[5]
print(paste(INPUT_FILE, "simulating", X, "x coverage, random seed", SEED, args[4], "to output:", OUTPUT))

ignoreLeft = 9 # 9 columns before the genotypes (2 per ancestry + 3 for snp position)
infile <- file(INPUT_FILE, 'r') 
outfile <- file(paste0(OUTPUT, ".", X, "x.", SEED, ".", args[4]), 'w') 
while (length(oneLine <- readLines(infile, n=1)) > 0){ # read one line at a time
  d = strsplit(x = oneLine, split = "\t")[[1]]
  other = d[1:(ignoreLeft)] # other (SNP & reference pops) information columns
  geno = as.integer(d[(ignoreLeft+1):length(d)]) # sample genotypes
  byInd = matrix(geno, ncol = 2, byrow = T)
  p = byInd[,1]/apply(byInd, 1, sum) # binomial probability based on genotype
  # note: there are some NAs for ind's without called genotypes at a position
  if (MODEL == "poisson"){
    # total coverage is taken from a poisson distribution
    tot = rpois(n = length(p), lambda = X)
  }
  if (MODEL == "fixed"){
    # total coverage is fixed (even) across genome for each bp = X cov.
    tot = rep(X, length(p))
  }
  if (MODEL == "negbinom"){
    # total coverage is taken from a negative binomial distribution (with set variance-to-mean ratio VMR)
    tot = rnbinom(n = length(p), size = X/(VMR - 1), mu = X) 
    # mean (mu) = X coverage; dispersion parameter (size) is parameterized s.t. var = mu(mu/size + 1)
    # so variance-to-mean ratio VMR = mu(mu/size + 1)/mu
  }
  readA1 = sapply(1:length(p), function(i)
    ifelse(is.na(p[i]),
           NA,
           rbinom(n = 1, size = tot[i], prob = p[i])))
  readA2 = tot - readA1
  # low sequencing error rates should have minimal effect, but are included
  errorsA1 = sapply(1:length(p), function(i)
    ifelse(is.na(p[i]),
           NA,
           rbinom(n = 1, size = readA1[i], prob = OBS_ERROR))
  )
  errorsA2 = sapply(1:length(p), function(i)
    ifelse(is.na(p[i]),
           NA,
           rbinom(n = 1, size = readA2[i], prob = OBS_ERROR))
  )
  readA1_rmError = readA1 - errorsA1 + # error reads don't get counted
    sapply(1:length(p), function(i) # plus on avg. one third of A2 errors become A1 reads
      ifelse(is.na(p[i]),
             NA,
             rbinom(n = 1, size = errorsA2[i], prob = 1/3))
    )
  readA2_rmError = readA2 - errorsA2 + # error reads don't get counted
    sapply(1:length(p), function(i) # plus one third of A1 errors become A2 reads
      ifelse(is.na(p[i]),
             NA,
             rbinom(n = 1, size = errorsA1[i], prob = 1/3))
    )  
  readRow = as.vector(t(cbind(readA1_rmError, readA2_rmError))) # puts in vector sample1_A1, sample1_A2, sample2_A1, sample2_A2 etc.
  readRow[is.na(readRow)] <- 0 # turn NAs into zeros (for both reads)
  writeLines(paste(c(other, readRow), collapse = "\t"), sep = "\n", con=outfile) 
}
close(infile)
close(outfile)
