# expected drift in ancestry based on Long 1991, with continuous migration:

# some parameter estimates:
# generations of admixture
t = 60
# starting admixture proportion
alpha_0 = 0.84
# ending admixture proportion
alpha_t = 0.4
# approximate Ne
Ne = 1000

# plot expectation if started at alpha_t for whole time,
# vs. equal influx migrants to get to alpha_t at time t
# classic
fst_noMig <- function(Ne, t){
  1 - exp(-t/(2*Ne))
}
# calculate migration rate per generation
calc_m <- function(alpha_0, alpha_t, t){
  #m = -log(alpha_0 - alpha_t)/t
  m = 1 - (alpha_t/alpha_0)^(1/t)
  return(m)
}

# add migration
fst_withMig <- function(Ne, t, alpha_0, m){
  # proportion A ancestry diluted by 
  # migration m each gen. from C/M ancestry
  # alpha_t = alpha_0 - (1 - m)^t
  
  # expected A ancestry over generations, with m mig.
  if (t == 0){
    fst = 0
  } else{

    
    # eq. 24 Long 1991 (? derivation)
    j = 1:t
    v = ((1 - 1/(2*Ne))^(j-1))*((1 - m)^2)^j
    alpha_j = alpha_0*(1 - m)^j
    fst = sum(j)/(2*Ne)
    #var_alpha = sum(j*alpha_j*(1-alpha_j))/(2*Ne)
    #fst = var_alpha/(alpha_j[t]*(1-alpha_j[t]))
  }
  return(fst)
}
# asymptote with migration
fst_asymp = function(m, Ne){
  fst = (1-m)^2/(2*Ne - (2*Ne - 1)*(1 - m)^2)
  return(fst)
}
t_max = 60
m = calc_m(alpha_0 = 0.84, alpha_t = 0.4, t = t_max + 1)
alphas = numeric(length = (t_max + 1))
alphas[1] <- 0.84
for (i in 1:t_max) {
  alphas[i + 1] = alphas[i]*(1-m)
}
plot(0:t_max, alphas)
plot(0:1000, sapply(0:1000, function(t) 
  fst_withMig(Ne = 1000, t = t, alpha_0 = 0.84, m = m)),
  col = "blue", xlab = "gen", ylab = "Fst")
curve(fst_noMig(Ne = 1000, t = x), from = 0, to = t_max,
      add = T)
fst_asymp(m = m, Ne = 1000)

# simulate the drift process?
sim_drift <- function(Ne, alpha_0, t, m = 0){
   As = numeric(length = t) # ancestry freq at a locus
   As[1] <- alpha_0
   for (i in 1:(t-1)){
     As[i + 1] = rbinom(n = 1, size = 2*Ne, prob = As[i])/(2*Ne)*(1 - m) # all migrants have other ancestry
   }
   return(As)
}

# number of simulations - precision
n_sim = 10^5

# fst
#fst1 = 1 - mean(A_ts*(1-A_ts))/(0.84*(1-0.84))
fst1 = 1 - apply(sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = 0)), 
                 1, function(r) mean(r*(1-r)))/(0.84*(1-0.84))
fst1_compare <- sapply(1:60, function(t) fst_noMig(Ne = 1000, t = t))
plot(fst1_compare, fst1)


fst0 = 1 - apply(sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.4, t = 60, m = 0)), 
                 1, function(r) mean(r*(1-r)))/(0.4*(1-0.4))
fst0_compare <- sapply(1:60, function(t) fst_noMig(Ne = 1000, t = t))
plot(fst0_compare, fst0)
plot(1:60, fst0, main = "Ast no mig. alpha = 0.4, t = 60",
     xlab = "generation t", ylab = "Ancestry Fst")
points(1:60, fst0)
# migration per gen
m = calc_m(alpha_0 = 0.84, alpha_t = 0.4, t = 60)
m
# expected ancestry per gen
alphas = 0.84*(1-m)^(1:60)
fst2 = 1 - apply(sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = m)), 
                 1, function(r) mean(r*(1-r)))/(alphas*(1-alphas))
plot(1:60, fst2, col = "blue")
points(1:60, fst2, col = "blue")

fst3 = 1 - apply(sapply(1:n_sim, function(x) sim_drift(Ne = 100, alpha_0 = 0.5, t = 60, m = m)), 
                 1, function(r) mean(r*(1-r)))/(0.5*(1-m)^(1:60)*(1-0.5*(1-m)^(1:60)))
points(1:60, fst3, col = "orange")


#sim2 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = m))
sim2_a = apply(sim2, 1, function(r) mean(2*r*(1-r)))
sim2_b = apply(sim2, 1, var)
sim2_c = sapply(1:nrow(sim2), function(i) mean((sim2[i, ]-alphas[i])^2))
plot(1:60, sim2_a)
plot(1:60, sim2_b)
plot(sim2_b, sim2_c)
abline(0, 1, col = "blue")
sim0 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = 0))
sim2 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = m))
sim4 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = m*2))
sim5 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.84, t = 60, m = m/2))
sim6 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.5, t = 60, m = 0))
sim3 = sapply(1:n_sim, function(x) sim_drift(Ne = 1000, alpha_0 = 0.5, t = 60, m = m))

sim3_b = apply(sim3, 1, var)
sim4_b = apply(sim4, 1, var)
sim5_b = apply(sim5, 1, var)
sim0_b = apply(sim0, 1, var)
sim6_b = apply(sim6, 1, var)
ms = c(0, m, m*2, m/2, 0, m)
a0s = c(rep(0.84, 4), rep(0.5, 2))
sims = list(sim0, sim2, sim4, sim5, sim6, sim3)
sim_cs = lapply(1:length(ms), function(s)
                sapply(1:60, function(t) 
                         mean((sims[[s]][t, ] - 
                                 (a0s[s]*(1-ms[s])^t))^2)))
png("plots/non-intuitive_relationship_drift_ancestry_variance.png")
plot(1:60, sim0_b, type = "l", ylab = "ancestry var", xlab = "t (gen)")
for (i in 1:5){
  lines(1:60, list(sim2_b, sim4_b, sim5_b, sim6_b, sim3_b)[[i]], 
        type = "l", col = c(rainbow(5)[[i]]))
}
for (i in 1:6){
  lines(1:60, sim_cs[[i]], 
        type = "l", lty = 2, col = c("black", rainbow(5))[[i]])
}
legend("topleft", legend = c("0.84, m=0", "0.84, m", "0.84, 2m", "0.84, m/2",
                             "0.5, m=0", "0.5, m"),
       title = "a0, m = 0.12",
       col = c("black", rainbow(5)), lty = 1)
dev.off()

plot(0:t_max, sapply(0:t_max, function(t) 
  fst_withMig(Ne = 10000000, t = t, alpha_0 = 0.84, alpha_t = 0.4)),
  col = "blue", xlab = "gen", ylab = "Fst")
curve(fst_noMig(Ne = 1000, t = x), from = 0, to = t_max,
      add = T)
fst_withMig(Ne = 1000, t = 0, alpha_0 = 0.84, alpha_t = 0.4)
m = -log(0.84 - 0.4)/t_max
fst_asymp(m = 0.0008, Ne = 1000)

# how do I correct for sampling variance?

# actual population ancestry at a locus:
n_large = 100000
anc1 <- rnorm(n = n_large, mean = 0.5, sd = 0.1)
# samples:
# e.g. sample 10 haplotypes
samples10 <- rbinom(n = n_large, size = 10, prob = anc1)
samples10_genomewide <- rbinom(n = n_large, size = 10, prob = 0.5)
samples20 <- rbinom(n = n_large, size = 20, prob = anc1)
samples20_genomewide <- rbinom(n = n_large, size = 20, prob = 0.5)
var(anc1)
var(samples10) - var(samples10_genomewide) # two processes aren't independent; can't sum variances.
# but, can still test against a null model var(samples10) == var(samples10_genomewide)
# because variance in anc1 will increase var(samples10)
var(samples10_genomewide)
10*.5*.5
mean((samples10_genomewide - mean(samples10_genomewide))^2) 
var(samples20_genomewide)
20*.5*.5
var(rbinom(n = n_large, size = 3, prob = 0.5))
3*.5*.5

load("../local_ancestry/results/zAnc.RData")
load("../local_ancestry/results/A.RData")
load("../local_ancestry/results/meta.RData")
load("../local_ancestry/results/sites_r.RData")
load("results/bootstrap_pi.RData") # bootstrap_pi for within-A ancestry diversity
colnames(zAnc_bees$K)
colnames(A)
meta.pop$population == colnames(A)
colnames(A) == colnames(zAnc_bees$K)
het_small_sample_correction <- function(p, n, filter_under_2 = T){ # this simplifies to just 2pq*n/(n-1)
  ifelse(filter_under_2 & n <= 2, NA, 2*(p - p^2*(n/(n-1)) + p*(1/(n-1))))}
ancestry_het <- sapply(1:nrow(meta.pop), function(i) 
  mean(sapply(A[ , i], function(x) 
    het_small_sample_correction(p = x, n = 2*meta.pop$n_bees[i]))))
Ast = 1 - ancestry_het/(2*zAnc_bees$alpha*(1-zAnc_bees$alpha))
binomial_var_anc = zAnc_bees$alpha*(1 - zAnc_bees$alpha)/(2*meta.pop$n_bees)
observed_var_anc = diag(zAnc_bees$K)
# observed ancestry variance vs. expected
plot(binomial_var_anc, # expected variance from binomial
     # do small sample size correction on to get sample variance
     observed_var_anc*(2*meta.pop$n_bees)/(2*meta.pop$n_bees - 1), 
     main = "Variance in ancestry expected vs. observed")
abline(a = 0, b = 1, col = "purple")
plot(abs(meta.pop$lat), Ast, main = "Ast across latitude")
abline(h = 0, col = "purple")
plot(2*zAnc_bees$alpha*(1-zAnc_bees$alpha), ancestry_het)
abline(a = 0, b = 1, col = "purple")

# observed vs. expected ancestry variance:
plot(binomial_var_anc, observed_var_anc, 
     xlab = "Expected ancestry variance (alpha*(1-alpha))/(2N)", 
     ylab = "Observed ancestry variance (diag K matrix)",
     main = "Observed ancestry variance is slightly more than binomial expectation")
points(binomial_var_anc, 
     2*meta.pop$n_bees/(2*meta.pop$n_bees - 1)*
       observed_var_anc, col = "blue")
abline(a = 0, b = 1, col = "blue")
colnames(meta.pop)

# sanity check on calculation of observed ancestry variance:
observed_var_anc3 = apply(apply(A, 1, function(row) (row - zAnc_bees$alpha)^2),
                          1, mean)
plot(observed_var_anc, observed_var_anc3)
abline(a = 0, b = 1, col = "green")

# let p = 0.5 in the true pop
ps_sample8 <- rbinom(n = 100, size = 8, p = 0.5)/8
2*0.5*0.5 # true het. at all loci
mean(2*ps_sample8*(1-ps_sample8)) # small sample het
mean(2*(ps_sample8 - ps_sample8^2*(8/(8-1)) + ps_sample8*(1/(8-1))))
test_small_sample <- function(l, n, p){ # l loci, n haplotypes, freq in pop p
  ps = rbinom(n = l, size = n, p = p)/n
  return(c(obs = mean(2*ps*(1-ps)),
           n_corr = mean(2*ps*(1-ps)*n/(n-1)),
           corrected = mean(2*(ps - ps^2*(n/(n-1)) + ps*(1/(n-1)))),
           true = 2*p*(1-p)))
}
test_small_sample(100, 8, 0.5)
test_small_sample(1000000, 4, 0.25)
hist(sapply(1:100, function(x) test_small_sample(10000, 3, 0.4)[2]))
abline(v = test_small_sample(10000, 3, 0.4)[3], col = "red")

# Ast vs. Fst from within-A ancestry pi estimates:
data.frame(population = colnames(zAnc_bees$K),
           Ast = Ast,
           stringsAsFactors = F) %>%
  left_join(meta.pop, ., by = "population") %>%
  left_join(., filter(bootstrap_pi, ancestry == "A"), by = c("population"="pop")) %>%
  mutate(pi_A_ref = bootstrap_pi[bootstrap_pi$pop == "A" & bootstrap_pi$ancestry == "A", "estimate"]) %>%
  mutate(Fst_within_A = 1 - estimate/pi_A_ref,
         Fst_within_A_lower = 1 - upper/pi_A_ref,
         Fst_within_A_upper = 1 - lower/pi_A_ref) %>%
  ggplot(., aes(x = Ast, y = Fst_within_A, color = exclude)) +
  geom_point() +
  geom_pointrange(aes(ymin = Fst_within_A_lower, ymax = Fst_within_A_upper)) +
  facet_grid(.~zone) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_fixed() +
  xlab("Ancestry Ast = 1 - ancestry het/(2*a*(1-a))") +
  ylab("Fst within A ancestry = 1 - pop het within A homozyg. / ref. A pop het") +
  ggtitle("Ancestry drift (Ast) vs. within A ancestry drift (Fst)") +
  scale_color_manual(values = c("blue", "orange"), labels = c("good pi estimate", "unreliable, low data"), name = element_blank())
ggsave("plots/Ast_vs_Fst.png", device = "png", height = 12, width = 7, units = "in", dpi = 300)  
  
plot(Ast, observed_var_anc/(zAnc_bees$alpha*(1-zAnc_bees$alpha)))
plot(2*zAnc_bees$alpha*(1-zAnc_bees$alpha)-ancestry_het, observed_var_anc)
abline(0,1)

# what heterozygosity in ancestry do I observe and what do I expect under simple binomial sampling?
het_obs_uncorrected = apply(2*A*(1-A), 2, mean)
het_binomial_only = sapply(1:nrow(meta.pop), function(i){
  a = rbinom(n = 100000, size = meta.pop$n_bees[i]*2, prob = zAnc_bees$alpha[i])/(2*meta.pop$n_bees[i])
  mean(2*a*(1-a))
})
png("plots/het_obs_expected_binomial_only.png", height = 4, width = 5, units = "in", res = 300)
plot(het_binomial_only, het_obs_uncorrected, main = "Ancestry het. observed vs. binomial only",
     xlab = "Het. expected from binomial sampling", ylab = "Het. observed (uncorrected for sample size)",
     pch = ifelse(meta.pop$population == "Avalon_2014", 2, 1))
abline(0, 1, col = "purple")
dev.off()
