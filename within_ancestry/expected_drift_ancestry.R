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





