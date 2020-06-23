# Kinship matrix using empirical pairwise ancestry covariance 
# (using individual genome-wide ancestry proportion alpha[i] as the 'mean' from which ancestry deviates)
# k[ij] = mean over loci l of (x[il] - alpha[i])*(x[jl] - alpha[j])

calcK = function(ancFreqMatrix, alpha){ # for all markers and # pops = length(alpha)
  ((ancFreqMatrix - alpha) %*% t(ancFreqMatrix - alpha))/ncol(ancFreqMatrix) #matrix multiply
  #above matrix operation is equivalent to the original code below
  #K <- matrix(nrow = sample_size, ncol = sample_size)
  #for (i in 1:sample_size) {
  #for (j in 1:sample_size) {
  #K[i,j] <- mean((anc[i,] - alpha[i])*(anc[j,] - alpha[j]))
  #}
  #}
}


calcInvL = function(K){
  solve(t(chol(K)))  
  #calculate the inverse Cholesky. t(chol(K)) is a lower left triangular matrix such that LL' = K
}

# function for calculating statistic:
make_K_calcs = function(pop_anc){
  # calculate mean population frequency across all snps
  pop_alpha = apply(pop_anc, 1, mean)
  # calculate K matrix for populations
  pop_K = calcK(ancFreqMatrix = pop_anc, alpha = pop_alpha)
  pop_InvL = calcInvL(pop_K)
  return(list(alpha = pop_alpha, K = pop_K, InvL = pop_InvL))
}

get_mean_from_K <- function(K, m = meta.pop){ # for covariance use K, for correlation set K = cov2cor(K)
# need meta.pop information to group populations
  lower_tri <- K
  lower_tri[lower.tri(lower_tri, diag = T)] <- NA # ignore variances and repeats
  AR_pops_S <- m$population[m$zone == "S. America" & m$region == "Low A"]
  AR_pops_N <- m$population[m$zone == "S. America" & m$region == "High A"]
  CA_pops <- m$population[m$zone == "N. America"]
  k_lower_tri <- melt(lower_tri) %>%
    filter(!is.na(value)) %>%
    mutate(type = ifelse(Var1 %in% CA_pops & Var2 %in% CA_pops, "CA_CA",
                         ifelse(Var1 %in% AR_pops_S & Var2 %in% AR_pops_S, "ARS_ARS",
                                ifelse(Var1 %in% AR_pops_N & Var2 %in% AR_pops_N, "ARN_ARN",
                                       ifelse((Var1 %in% AR_pops_N & Var2 %in% AR_pops_S) | (Var1 %in% AR_pops_S & Var2 %in% AR_pops_N), "ARN_ARS",
                                              ifelse((Var1 %in% AR_pops_S & Var2 %in% CA_pops) | (Var1 %in% CA_pops & Var2 %in% AR_pops_S), "CA_ARS",
                                                     ifelse((Var1 %in% AR_pops_N & Var2 %in% CA_pops) | (Var1 %in% CA_pops & Var2 %in% AR_pops_N), "CA_ARN", 
                                                            NA)))))))
  mean_corr <- k_lower_tri %>%
    group_by(type) %>%
    summarise(mean_anc_corr = mean(value))
  return(mean_corr)
}
