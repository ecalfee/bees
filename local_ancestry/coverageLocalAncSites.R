# what is the distribution of coverage across sites used for local ancestry inference?
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/AR0302.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/CA0502.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/AR1403.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/CA1120.counts.txt", header = F, stringsAsFactors = F)
a <- read.table("results/SNPs/thin1kb_common3/countsMajMin/AR2603.counts.txt", header = F, stringsAsFactors = F)
tot <- apply(a, 1, sum)
mean(tot)
hist(rpois(n=length(tot), lambda = mean(tot)), freq = F, col = "pink")
hist(rnbinom(n = length(tot), size = mean(tot)/(3 - 1), mu = mean(tot)), freq = F,
          col = "yellow", add = T)
hist(tot, freq = F, add = T, col = "blue")
hist(rnbinom(n = length(tot), size = mean(tot)/(2 - 1), mu = mean(tot)), freq = F,
     col = "green", add = T, bins = 20)
mean(tot)
var(tot) # true variance in coverage
mean(tot) + mean(tot)*(3-1) # expected variance with VMR=3
var(rnbinom(n = length(tot), size = mean(tot)/(3 - 1), mu = mean(tot))) # observed variance
mean(tot) + mean(tot)*(2-1) # expected variance with VMR=2
var(rnbinom(n = length(tot), size = mean(tot)/(2 - 1), mu = mean(tot))) # observed variance

# so it looks like at least for each of this small set of bees with varying coverage from ~1 to ~11x mean coverage,
# their distribution of coverage across included SNPs has mean variance and distribution similar to a negative binomial with VMR=2 
# (slightly less variable that what I modeled for neg binom but more variable than a poisson distribution of coverage)
    