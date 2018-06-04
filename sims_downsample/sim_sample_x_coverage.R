#!/usr/bin/env Rscript

# this script takes in an input file for ancestry_hmm, with sample genotypes,
# and simulates binomially sampled reads to ?x coverage across the genome
# using a user-specified random seed 
# and either "fixed" or "poisson" to indicate whether coverage should be fixed (even across every bp)
# or sampled from a poisson distribution
# to run: Rscript sim_sample_x_coverage.R ../data/TEST/head_ACM_Riv2014.counts 4 1353323 fixed

# to test
#args <- c("../data/TEST/head_ACM_Riv2014.counts", "4", "1353323", "fixed")
#args <- c("../data/TEST/head_ACM_Riv2014.counts", "4", "1353323", "poisson")
# user input:
args <- commandArgs(trailingOnly=TRUE)
INPUT_FILE = args[1]
X = as.numeric(args[2]) # desired coverage to simulate
SEED = as.integer(args[3])
set.seed(SEED) # set random seed
POISSON = ifelse(args[4]=="poisson", T, ifelse(args[4]=="fixed", F, stop("4th argument must be poisson or fixed")))
print(paste(INPUT_FILE, "simulating", X, "x coverage, random seed", SEED, args[4]))

ignoreLeft = 9 # 9 columns before the genotypes (2 per ancestry + 3 for snp position)
infile <- file(INPUT_FILE, 'r') 
outfile <- file(paste0(INPUT_FILE, ".", X, "x.", SEED, ".", args[4]), 'w') 
while (length(oneLine <- readLines(infile, n=1)) > 0){ # read one line at a time
  d = strsplit(x = oneLine, split = "\t")[[1]]
  other = d[1:(ignoreLeft)] # other (SNP & reference pops) information columns
  geno = as.integer(d[(ignoreLeft+1):length(d)]) # sample genotypes
  byInd = matrix(geno, ncol = 2, byrow = T)
  p = byInd[,1]/apply(byInd, 1, sum) # binomial probability based on genotype
  # note: there are some NAs for ind's without called genotypes at a position
  if (POISSON){
    # total coverage is taken from a poisson distribution
    tot = rpois(n = length(p), lambda = X)
  } else{
    # total coverage is fixed (even) across genome for each bp = X cov.
    tot = rep(X, length(p))
  }
  readA1 = sapply(1:length(p), function(i)
    ifelse(is.na(p[i]),
           NA,
           rbinom(n = 1, size = tot[i], prob = p[i])))
  readA2 = tot - readA1
  readRow = as.vector(t(cbind(readA1, readA2))) # puts in vector sample1_A1, sample1_A2, sample2_A1, sample2_A2 etc.
  readRow[is.na(readRow)] <- 0 # turn NAs into zeros (for both reads)
  writeLines(paste(c(other, readRow), collapse = "\t"), sep = "\n", con=outfile) 
}
close(infile)
close(outfile)
