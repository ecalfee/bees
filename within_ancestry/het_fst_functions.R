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

# # For genomewide mean Fst, I take the ratio of the avg. numerator & avg. denominator) 
# calc_avg_hudson_Fst <- function(v1, v2, n1, n2){ # allele frequencies and number of alleles sampled
#   exclude = is.na(v1) | is.na(v2)
#   p_1 = v1[!exclude]
#   p_2 = v2[!exclude]
#   n_1 = n1[!exclude]
#   n_2 = n2[!exclude]
#   
#   N = (p_1 - p_2)^2 - p_1*(1 - p_1)/(n_1 - 1) - p_2*(1-p_2)/(n_2 - 1)
#   D = p_1*(1 - p_2) + (1 - p_1)*p_2
#   
#   fst = mean(N)/mean(D)
#   return(fst)
# }