FDR_values = c(.1, .05, .01)
# false discovery rate calculations
fdr_shared_high <- function(a1, a2, pop1, pop2, sims1, sims2){# takes in an ancestry, test for high ancestry in both zones
  obs <- sum(pop1 >= a1 & pop2 >= a2) # observations exceeding threshold
  null <- sum(sims1 >= a1 & sims2 >= a2)/length(sims1)*length(pop1) # expected number of neutral loci exceeding threshold in data set size length(pop1)
  null/obs
}
fdr_shared_low <- function(a1, a2, pop1, pop2, sims1, sims2){# takes in an ancestry, test for low ancestry in both zones
  obs <- sum(pop1 <= a1 & pop2 <= a2)
  null <- sum(sims1 <= a1 & sims2 <= a2)/length(sims1)*length(pop1)
  null/obs
}
fdr_1pop_high <- function(a, pop, sims){# takes in an ancestry, test for high ancestry in 1 zone
  obs <- sum(pop >= a)
  null <- sum(sims >= a)/length(sims)*length(pop)
  null/obs
}
fdr_1pop_low <- function(a, pop, sims){# takes in an ancestry, test for low ancestry in 1 zone
  obs <- sum(pop <= a)
  null <- sum(sims <= a)/length(sims)*length(pop)
  null/obs
}