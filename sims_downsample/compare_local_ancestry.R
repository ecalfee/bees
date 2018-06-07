# compares local ancestry calls for those inferred from high coverage
# sequence data with inferred genotypes to those inferred from
# low coverage data (simulated binomial sampling of alleles -> reads)

id = "SRCD47A"
id = "SRCD10B"
id = "SRCD9A"
# high coverage 'true ancestry'
d = read.table(paste0("../data/sims_downsample/ancestry_hmm/est_t2/",
id, ".posterior"), 
               stringsAsFactors = F, header = T)

s = 1001
time = "fix_t" #est_t2
colnames(d)[3:8] 
ancTypes = c("AA", "AC", "AM", "CC", "CM", "MM") 
postH = .9 # posterior minimum for high confidence ancestry calls
postT = .9 # posterior to call a true call
postF = .5 # posterior to call a false call
# 70% genome is called confidently for ancestry
sum(d$X2.0.0>=postH | d$X1.1.0>=postH |
      d$X0.2.0>=postH | d$X0.1.1>=postH |
      d$X0.0.2>=postH | d$X1.0.1>=postH)/nrow(d)

calcPwr = function(dTrue, dSim, ancN){
  colN = ancN+2
  # true positives
  tp = sum(dSim[dTrue[ , colN] >= postH, colN] >= postT)/sum(dTrue[ , colN]>=postH)
  # false positives
  fp = (sum(apply(dSim[dTrue[ , colN]>=postH, 3:8], 1, max) >= postF) - sum(dSim[dTrue[ , colN]>=postH, colN] >= postF))/sum(dTrue[ , colN]>=postH)
  return(c(tp=tp, fp=fp))  
}
poisX = lapply(1:10, function(x)
  read.table(paste0("../data/sims_downsample/ancestry_hmm/", x, 
                                       "x/", "poisson", "/seed", s, "/", time, "/", id, ".posterior"), 
                                stringsAsFactors = F, header = T))

poisPwr = lapply(poisX, function(pois) sapply(1:6, function(i) 
  calcPwr(dSim = pois, dTrue = d, ancN = i)))

# calculate power for fixed coverage results
fixedX = lapply(1:10, function(x)
  read.table(paste0("../data/sims_downsample/ancestry_hmm/", x, 
                    "x/", "fixed", "/seed", s, "/", time, "/", id, ".posterior"), 
             stringsAsFactors = F, header = T))

fixedPwr = lapply(fixedX, function(fixed) sapply(1:6, function(i) 
  calcPwr(dSim = fixed, dTrue = d, ancN = i)))

# plot & save results for poisson simulation
png(paste0("../plots/sim_downsample_pwr_", id, ".png"),
     height = 6, width = 8, units = "in", res = 120)
plot( # make empty plot frame
    1:10, rep(0, 10), col = NULL,
    ylim = c(0,1), ylab = paste0("Proportion correct ancestry (post>", postT,")"),
    xlab = paste0("simulated coverage"),
    main = paste0("Proportion high conf. (post>", postH, ") ancestry calls replicated in low-coverage"),
    sub = paste("False-positive rate (post>", postF, ") for incorrect ancestry plotted at bottom") 
) # add lines 
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) poisPwr[[x]][1, i]), 
        col = rainbow(6)[i], lwd = 2) 
}
legend("right", legend = ancTypes, col = rainbow(6), pch = 20, title = "True Ancestry")

# add lines for fixed coverage
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) fixedPwr[[x]][1, i]), 
        col = rainbow(6)[i], lwd = 2, lty = 2) 
}
legend("bottomright", legend = c("Poisson", "Fixed at mean"), 
       title = "Coverage distribution", lwd = 2, lty = c(1, 2), col = "black")

# optional - add lines for false positives
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) poisPwr[[x]][2, i]), 
        col = rainbow(6)[i], lwd = 1) 
}
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) fixedPwr[[x]][2, i]), 
        col = rainbow(6)[i], lwd = 1, lty = 2) 
}
dev.off()
