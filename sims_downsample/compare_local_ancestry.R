library(ggplot2)
library(tidyr)
library(dplyr)

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
continents = c("AfrAfr", "AfrEu", "AfrEu", "EuEu", "EuEu", "EuEu")
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


# make 'confusion histogram' for specific x coverage
# plots for 4x coverage
conf4x_list = lapply(1:length(ancTypes), function(i) poisX[[4]][d[ , i+2] >= postH, ])
for (i in 1:length(ancTypes)){
  colnames(conf4x_list[[i]]) = c("chrom", "position", ancTypes)
  conf4x_list[[i]] = tidyr::gather(conf4x_list[[i]], key = "ancestryPred", value = "posterior", ancTypes)
  conf4x_list[[i]]$ancestryTrue = ancTypes[i]
}
conf4x = do.call(rbind, conf4x_list)
filter(conf4x, ancestryTrue == ancestryPred) %>%
  ggplot(., aes(posterior, fill = ancestryPred)) + geom_density(alpha = 0.5) +
  ggtitle("posterior in low-cov data for true ancestry")
conf4x[conf4x$posterior >= .001,] %>% # only plot over .001 posterior because too many at low post.
  ggplot(., aes(posterior, fill = ancestryPred)) + geom_density(alpha = 0.5) +
  facet_wrap(~ ancestryTrue)

# identify true and predicted continents of Ancestry
# i.e. African or European homozygous or heterozygous
conf4x$continentPred = sapply(conf4x$ancestryPred, function(anc)
  continents[which(ancTypes == anc)])
conf4x$continentTrue = sapply(conf4x$ancestryTrue, function(anc)
  continents[which(ancTypes == anc)])
# group data by continent-level ancestry, summing across ancestry within continental groupings
conf4x_continent = group_by(conf4x, continentTrue, continentPred, chrom, position) %>%
  summarise(., contPost = sum(posterior))

filter(conf4x_continent, continentTrue == continentPred) %>%
  ggplot(., aes(contPost, fill = continentPred)) + geom_density(alpha = 0.5) +
  ggtitle("posterior in low-cov data for true continental ancestry")
conf4x_continent[conf4x_continent$contPost >= .001,] %>% # only plot over .001 posterior because too many at low post.
  ggplot(., aes(contPost, fill = continentPred)) + geom_density(alpha = 0.5) +
  facet_wrap(~ continentTrue)
# this is hard to visualize, so I'll plot a 'confusion matrix'
# where the low-confidence classification is the one with the highest posterior
# and the high-confidence classification is > 0.9
# and there is an 'unclassified' category for loci that don't meet the threshold for any ancestry classification




