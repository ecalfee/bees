# snp density
chr_length_tot <- sum(read.table("../data/honeybee_genome/chr.lengths")$V2)
gaps_tot <- read.table("../data/honeybee_genome/gaps.bed") %>%
  mutate(length = V3 - V2) %>%
  summarise(tot = sum(length)) %>%
  unlist(.)
n_snps <- 3510834 # from wc -l chr.var.sites
# genome_size <- chr_length_tot - gaps_tot #236*10^6 # genome size
# actually, use all sites passing quality filters for the denominator of pi:
all_sites <- data.frame(
  bin_name = read.table("../geno_lik_and_SNPs/results/1cM_bins.names", header = F, stringsAsFactors = F)$V1,
  n_sites = read.table("../geno_lik_and_SNPs/results/1cM_bins_total_counts.txt",
                       header = F, sep = "\t")$V2)
genome_size <- sum(all_sites$n_sites)
frac_snps <- n_snps/genome_size # multiplier for any snp-based heterozygosity estimate
# fraction of the genome that is variable, i.e. snps


# some functions:
# heterozygosity corrected for small sample size:
# by default also drop any snps with fewer than 2 individuals with data (<= 2 alleles observed, i.e. 1 ind.)
het_small_sample_correction_nind <- function(p, n_ind, min_ind = 2){ # assumes diploid individuals
  ifelse(n_ind < min_ind, NA, 2*(p - p^2*(2*n_ind/(2*n_ind-1)) + p*(1/(2*n_ind-1))))}

het_small_sample_correction <- function(p, n){ # here n is the number of haplotypes observed
  2*(p - p^2*(n/(n-1)) + p*(1/(n-1)))}


# Hudson pairwise Fst estimator from Bhatia 2013 
# (I implement equation 10 from the methods) 
hudson_Fst <- function(p1, p2, n_ind1, n_ind2, min_ind = 2){ # allele frequencies and number of alleles sampled
  if (n_ind1 < min_ind | n_ind2 < min_ind){
    fst = NA
  } else{
    N = (p1 - p2)^2 - p1*(1 - p1)/(n_ind1*2 - 1) - p2*(1-p2)/(n_ind2*2 - 1)
    D = p1*(1 - p2) + (1 - p1)*p2
    fst = N/D
  }
  return(fst)
}

hudson_Fst_window <- function(p1, p2, n_ind1, n_ind2, min_ind = 2, min_snps = 10){ # VECTORS of allele frequencies and number of alleles sampled
  x <- data.frame(p1 = p1, p2 = p2, n_ind1 = n_ind1, n_ind2 = n_ind2, stringsAsFactors = F) %>%
    filter(complete.cases(.)) %>%
    filter(n_ind1 >= min_ind & n_ind2 >= min_ind)
  if (nrow(x) < min_snps){
    fst = NA
  } else{
    N = (x$p1 - x$p2)^2 - x$p1*(1 - x$p1)/(x$n_ind1*2 - 1) - x$p2*(1-x$p2)/(x$n_ind2*2 - 1)
    D = x$p1*(1 - x$p2) + (1 - x$p1)*x$p2
    fst = mean(N)/mean(D)
  }
  return(fst)
}
hudson_Fst_window2 <- function(p1, p2, n_ind1, n_ind2, min_ind = 2, min_snps = 10){ # VECTORS of allele frequencies and number of alleles sampled
  x <- data.frame(p1 = p1, p2 = p2, n_ind1 = n_ind1, n_ind2 = n_ind2, stringsAsFactors = F) %>%
    filter(complete.cases(.)) %>%
    filter(n_ind1 >= min_ind & n_ind2 >= min_ind) %>%
    filter((p1 < 1 & p1 > 0) & (p2 < 1 & p2 > 0))
  if (nrow(x) < min_snps){
    fst = NA
  } else{
    N = (x$p1 - x$p2)^2 - x$p1*(1 - x$p1)/(x$n_ind1*2 - 1) - x$p2*(1-x$p2)/(x$n_ind2*2 - 1)
    D = x$p1*(1 - x$p2) + (1 - x$p1)*x$p2
    fst = mean(N)/mean(D)
  }
  return(fst)
}
#hudson_Fst_window(p1 = c(.5, .5, .5), p2 = c(1, .5, 1), n_ind1 = rep(5, 3), n_ind2 = rep(5, 3), min_snps = 1)
#hudson_Fst_window(p1 = c(.5, .5, .5), p2 = c(1, 1, 1), n_ind1 = rep(5, 3), n_ind2 = c(5, 5, 1), min_snps = 3)

# For genomewide mean Fst, I take the ratio of the avg. numerator & avg. denominator) 
calc_avg_hudson_Fst <- function(v1, v2, n1, n2){ # allele frequencies and number of alleles sampled
  exclude = is.na(v1) | is.na(v2)
  p_1 = v1[!exclude]
  p_2 = v2[!exclude]
  n_1 = n1[!exclude]
  n_2 = n2[!exclude]
  
  N = (p_1 - p_2)^2 - p_1*(1 - p_1)/(n_1 - 1) - p_2*(1-p_2)/(n_2 - 1)
  D = p_1*(1 - p_2) + (1 - p_1)*p_2
  
  fst = mean(N)/mean(D)
  return(fst)
}