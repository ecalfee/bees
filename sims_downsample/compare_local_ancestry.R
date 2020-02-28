library(ggplot2)
library(tidyr)
library(dplyr)

# compares local ancestry calls for those inferred from high coverage
# sequence data with inferred genotypes to those inferred from
# low coverage data (simulated binomial sampling of alleles -> reads)

#id = "SRCD47A"
#id = "SRCD10B"
#id = "SRCD9A"
id = "Riv2014"
ids = paste0("SRCD", c("6A", "9A", "10B", "47A", "49A", "51B", "51C", "54C"))
# high coverage 'true ancestry'
d = do.call(rbind, lapply(ids, function(id) 
  read.table(paste0("../data/sims_downsample/ancestry_hmm/CMA/highx/", 
                   id, ".posterior"), stringsAsFactors = F, header = T) %>%
    mutate(Bee_ID = id)))

s = 100
#time = "fix_t" 
time = "est_t2"
colnames(d)[3:8] 
ancTypes = c("CC", "CM", "CA", "MM", "MA", "AA")
posterior_columns <- c("chrom", "position", ancTypes, "Bee_ID")
continents = c("AfrAfr", "AfrEu", "EuEu")
continents_long = c("EuEu", "EuEu", "AfrEu", "EuEu", "AfrEu", "AfrAfr")
postH = .9 # posterior minimum for high confidence ancestry calls
postT = .9 # posterior to call a true call
#postH = 0.8
#postT = 0.8
postF = .5 # posterior to call a false call
dir_prefix = "../data/sims_downsample/ancestry_hmm/CMA/"


# get and summarise data for high coverage
d_high <- d %>%
  data.table::setnames(posterior_columns) %>%
  mutate(A = AA + .5*CA + .5*MA) %>% # best estimate of A ancestry at a locus
  pivot_longer(., cols = ancTypes, names_to = "ancestry", values_to = "posterior") %>%
  group_by(chrom, position, Bee_ID, A) %>%
  summarise(ancestry_state = ancestry[which.max(posterior)],
            posterior = max(posterior)) %>%
  arrange(chrom, position, Bee_ID)

# and lower coverage results
d_10x <- lapply(1:10, function(x)
  do.call(rbind, lapply(ids, function(id) 
    read.table(paste0(dir_prefix, x, 
                      "x/negbinom/seed100/est_t2/", id, ".posterior"), 
               stringsAsFactors = F, header = T) %>%
      mutate(Bee_ID = id) %>%
      data.table::setnames(posterior_columns) %>%
      mutate(A = AA + .5*CA + .5*MA) %>%
      pivot_longer(., cols = ancTypes, names_to = "ancestry", values_to = "posterior") %>%
      group_by(chrom, position, Bee_ID, A) %>%
      summarise(ancestry_state = ancestry[which.max(posterior)],
                posterior = max(posterior)))) %>%
      arrange(chrom, position, Bee_ID))

save(file = "results/d_1-10x_high_CMA.RData", list = c("d_high", "d_10x"))

# % high confidence sites:
sum(d_high$posterior >= .9)/nrow(d_high)
sum(d_high$posterior >= .8)/nrow(d_high) # > 81.6% of the genome has a high conf. anc. call
d_high %>%
  group_by(ancestry_state) %>%
  summarise(min = min(A),
            max = max(A),
            mean = mean(A))

# A ancestry correlations:
d_corr <- sapply(d_10x, function(x) cor(d_high$A, x$A))
d_corr
# % replicated:
d_high_conf <- table(d_high$ancestry_state[d_high$posterior >= 0.8])
d_correct <- lapply(d_10x, function(x) table(d_high$ancestry_state[x$ancestry_state == d_high$ancestry_state & x$posterior >= 0.8 & d_high$posterior >= 0.8]))
d_incorrect <- lapply(d_10x, function(x) table(d_high$ancestry_state[x$ancestry_state != d_high$ancestry_state & x$posterior >= 0.8 & d_high$posterior >= 0.8]))

results <- do.call(rbind, lapply(1:10,
                                 function(i) d_correct[[i]]/d_high_conf)) %>%
  data.frame() %>%
  mutate(type = "correct",
         x = 1:10) %>%
  bind_rows(., do.call(rbind, lapply(1:10, 
                                     function(i) d_incorrect[[i]]/d_high_conf)) %>%
              data.frame() %>%
              mutate(type = "error",
                     x = 1:10)) %>%
  pivot_longer(cols = ancTypes, names_to = "ancestry", values_to = "p") 

# plot
labels_type <- c("Correct\n(same ancestry)", "Incorrect\n(diff. ancestry)")
names(labels_type) = c("correct", "error")
pr <- ggplot(results, aes(x = x, y = p, color = ancestry, lty = ancestry)) +
  geom_line() +
  facet_wrap(~type, labeller = labeller(type = labels_type)) +
  theme_light() +
  scale_color_viridis_d(option = "viridis", name = "Ancestry") +
  scale_linetype_manual(values = c(1, 2, 1, 2, 2, 1), name = "Ancestry") +
  scale_x_continuous(breaks = seq(2, 10, by = 2), labels = seq(2, 10, by = 2)) +
  xlab("Coverage") +
  ylab("% High confidence ancestry calls") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
pr
pc <- data.frame(x = 1:10,
           p = d_corr,
           type = "Correlation\n") %>%
  ggplot(., aes(x = x, y = p)) +
  geom_line() +
  theme_light() +
  facet_wrap(~type)+
  scale_x_continuous(breaks = seq(2, 10, by = 2), labels = seq(2, 10, by = 2)) +
  xlab("Coverage") +
  ylab("African ancestry correlation") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
pc  

pcr <- ggarrange(pc + ggtitle("A"), 
          pr + ggtitle("B") + theme(strip.background = element_rect(color = NA)), 
          nrow = 1, widths = c(1, 2))
pcr
ggsave("plots/sims_downsample.png",
       plot = pcr,
       width = 5.2, height = 4, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures/sims_downsample.png",
       plot = pcr,
       width = 5.2, height = 4, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/sims_downsample.tiff",
       plot = pcr,
       width = 5.2, height = 4, 
       units = "in", dpi = 600, device = "tiff")



# ~72% genome is called confidently for ancestry
sum(d$X2.0.0>=postH | d$X1.1.0>=postH |
      d$X0.2.0>=postH | d$X0.1.1>=postH |
      d$X0.0.2>=postH | d$X1.0.1>=postH)/nrow(d)
# functions
# calculate power to recover in low coverage data ancestry for positions with posteriors over postH in high-cov. data
calcPwr = function(dTrue, dSim, ancN){
  colN = ancN+2
  # true positives
  tp = sum(dSim[dTrue[ , colN] >= postH, colN] >= postT)/sum(dTrue[ , colN]>=postH)
  # false positives
  fp = (sum(apply(dSim[dTrue[ , colN]>=postH, 3:ncol(dSim)], 1, max) >= postF) - sum(dSim[dTrue[ , colN]>=postH, colN] >= postF))/sum(dTrue[ , colN]>=postH)
  return(c(tp=tp, fp=fp))  
}
calcPwr = function(dTrue, dSim, ancestry){ # e.g ancestry = 'AA'
  # true positives
  tp = sum(dSim[dTrue[ , colN] >= postH, colN] >= postT)/sum(dTrue[ , ancestry] >= postH)
  # false positives
  fp = (sum(apply(dSim[dTrue[ , colN]>=postH, 3:ncol(dSim)], 1, max) >= postF) - sum(dSim[dTrue[ , colN]>=postH, colN] >= postF))/sum(dTrue[ , colN]>=postH)
  return(c(tp=tp, fp=fp))  
}
# convert population level ancestry posteriors to continent level ancestry posteriors:
ACM_to_continent = function(posteriors, col_names = posterior_columns){
  colnames(posteriors) <- col_names
  post_continent = mutate(posteriors, AfrAfr = AA) %>%
    mutate(., AfrEu = CA + MA) %>%
    mutate(., EuEu = CC + CM + MM) %>%
    select(., chrom, position, AfrAfr, AfrEu, EuEu)
  return(post_continent)
}


# calculate power for negative binomial coverage results
negbinX = lapply(1:10, function(x)
  do.call(rbind, lapply(ids, function(id) 
    read.table(paste0(dir_prefix, x, 
                      "x/", "negbinom", "/seed", s, "/", time, "/", id, ".posterior"), 
               stringsAsFactors = F, header = T))))

negbinPwr = lapply(negbinX, function(x) sapply(1:6, function(i) 
  calcPwr(dSim = x, dTrue = d, ancN = i)))






################# OLD ###################################################


poisX = lapply(1:10, function(x)
  do.call(rbind, lapply(ids, function(id) 
    read.table(paste0(dir_prefix, x, 
                      "x/", "poisson", "/seed", s, "/", time, "/", id, ".posterior"), 
               stringsAsFactors = F, header = T))))

poisPwr = lapply(poisX, function(x) sapply(1:6, function(i) 
  calcPwr(dSim = x, dTrue = d, ancN = i)))

# calculate power for fixed coverage results
fixedX = lapply(1:10, function(x)
  do.call(rbind, lapply(ids, function(id) 
    read.table(paste0(dir_prefix, x, 
                      "x/", "fixed", "/seed", s, "/", time, "/", id, ".posterior"), 
               stringsAsFactors = F, header = T))))

fixedPwr = lapply(fixedX, function(x) sapply(1:6, function(i) 
  calcPwr(dSim = x, dTrue = d, ancN = i)))

# calculate power for negative binomial coverage results
negbinX = lapply(1:10, function(x)
  do.call(rbind, lapply(ids, function(id) 
    read.table(paste0(dir_prefix, x, 
                      "x/", "negbinom", "/seed", s, "/", time, "/", id, ".posterior"), 
               stringsAsFactors = F, header = T))))

negbinPwr = lapply(negbinX, function(x) sapply(1:6, function(i) 
  calcPwr(dSim = x, dTrue = d, ancN = i)))


# plot & save results for poisson simulation
png(paste0("../plots/sims_downsample/wSeqError/pwr_ACM_", id, ".png"),
     height = 6, width = 10, units = "in", res = 150)
plot( # make empty plot frame
    1:10, rep(0, 10), col = NULL,
    ylim = c(0,1), ylab = paste0("Proportion correct ancestry (post>", postT,")"),
    xlab = paste0("simulated coverage"),
    main = paste0("Proportion high conf. (post>", postH, ") ancestry calls replicated in low-coverage"),
    sub = paste("False-positive rate (post>", postF, ") for incorrect ancestry plotted at bottom") 
) 
# add lines for poisson 
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) poisPwr[[x]][1, i]), 
        col = rainbow(6)[i], lwd = 2) 
}

# add lines for fixed coverage
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) fixedPwr[[x]][1, i]), 
        col = rainbow(6)[i], lwd = 2, lty = 2) 
}

# add lines for neg. binomial
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) negbinPwr[[x]][1, i]), 
        col = rainbow(6)[i], lwd = 2, lty = 3) 
}

# add lines with estimated times:
#for (i in 1:6){
#  lines(1:10, sapply(1:10, function(x) negbinPwr_est_t2[[x]][1, i]), 
#        col = rainbow(6)[i], lwd = 1) 
#}


# add legend
legend("right", legend = ancTypes, col = rainbow(6), pch = 20, title = "True Ancestry")
legend("bottomright", legend = c("Poisson", "Fixed at mean", "Neg. Binom."), 
       title = "Coverage distribution", lwd = 2, lty = c(1, 2, 3), col = "black")

# optional - add lines for false positives
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) poisPwr[[x]][2, i]), 
        col = rainbow(6)[i], lwd = 1) 
}
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) fixedPwr[[x]][2, i]), 
        col = rainbow(6)[i], lwd = 1, lty = 2) 
}
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) negbinPwr[[x]][2, i]), 
        col = rainbow(6)[i], lwd = 1, lty = 3) 
}
dev.off()

# now a similar analysis at the continent (Afr-Eu) level:
# 'true ancestry' coded for continent-level analysis (groups C & M together as 'European')
dCont = ACM_to_continent(d)

c_poisX = lapply(poisX, ACM_to_continent)

c_poisPwr = lapply(c_poisX, function(x) sapply(1:3, function(i) 
  calcPwr(dSim = x, dTrue = dCont, ancN = i)))

# calculate power for fixed coverage results
c_fixedX = lapply(fixedX, ACM_to_continent)

c_fixedPwr = lapply(c_fixedX, function(x) sapply(1:3, function(i) 
  calcPwr(dSim = x, dTrue = dCont, ancN = i)))

# calculate power for negative binomial coverage results
c_negbinX = lapply(negbinX, ACM_to_continent)

c_negbinPwr = lapply(c_negbinX, function(x) sapply(1:3, function(i) 
  calcPwr(dSim = x, dTrue = dCont, ancN = i)))

# plot resulting power for continent-level:
png(paste0("../plots/sims_downsample/wSeqError/pwr_AfrEu_", id, ".png"),
    height = 6, width = 10, units = "in", res = 150)
plot( # make empty plot frame
  1:10, rep(0, 10), col = NULL,
  ylim = c(0,1), ylab = paste0("Proportion correct ancestry (post>", postT,")"),
  xlab = paste0("simulated coverage"),
  main = paste0("Proportion high conf. (post>", postH, ") ancestry calls replicated in low-coverage"),
  sub = paste("False-positive rate (post>", postF, ") for incorrect ancestry plotted at bottom") 
) 
# add lines for poisson 
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_poisPwr[[x]][1, i]), 
        col = rainbow(3)[i], lwd = 2) 
}

# add lines for fixed coverage
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_fixedPwr[[x]][1, i]), 
        col = rainbow(3)[i], lwd = 2, lty = 2) 
}

# add lines for neg. binomial
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_negbinPwr[[x]][1, i]), 
        col = rainbow(3)[i], lwd = 2, lty = 3) 
}

# add legend
legend("right", legend = continents, col = rainbow(3), pch = 20, title = "True Ancestry")
legend("bottomright", legend = c("Poisson", "Fixed at mean", "Neg. Binom."), 
       title = "Coverage distribution", lwd = 2, lty = c(1, 2, 3), col = "black")

# optional - add lines for false positives
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_poisPwr[[x]][2, i]), 
        col = rainbow(3)[i], lwd = 1) 
}
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_fixedPwr[[x]][2, i]), 
        col = rainbow(3)[i], lwd = 1, lty = 2) 
}
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_negbinPwr[[x]][2, i]), 
        col = rainbow(3)[i], lwd = 1, lty = 3) 
}
dev.off()


# just negative binomial (for higher readability):
# plot resulting power for continent-level:
png(paste0("../plots/sims_downsample/wSeqError/pwr_negbin_", id, ".png"),
    height = 6, width = 10, units = "in", res = 150)
plot( # make empty plot frame
  1:10, rep(0, 10), col = NULL,
  ylim = c(0,1), ylab = paste0("Proportion correct ancestry (post>", postT,")"),
  xlab = paste0("simulated coverage"),
  main = paste0("Proportion high conf. (post>", postH, ") ancestry calls replicated w/ negBin low-coverage"),
  sub = paste("False-positive rate (post>", postF, ") for incorrect ancestry plotted at bottom") 
) 

# add lines for neg. binomial - continent level
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_negbinPwr[[x]][1, i]), 
        col = rainbow(9)[i+6], lwd = 2, lty = 1) 
}

# add lines for neg. binomial - ancestry level
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) negbinPwr[[x]][1, i]), 
        col = rainbow(9)[i], lwd = 2, lty = 2) 
}

# add lines for false positives
for (i in 1:3){
  lines(1:10, sapply(1:10, function(x) c_negbinPwr[[x]][2, i]), 
        col = rainbow(9)[i+6], lwd = 1, lty = 3) 
}
for (i in 1:6){
  lines(1:10, sapply(1:10, function(x) negbinPwr[[x]][2, i]), 
        col = rainbow(9)[i], lwd = 1, lty = 3) 
}


# add legend
legend("right", legend = c(ancTypes, continents, "(false pos.)"), col = c(rainbow(9), "black"), pch = 20, title = "True Ancestry",
       lty = c(rep(2, 6), rep(1, 3), 3))


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
  continents_long[which(ancTypes == anc)])
conf4x$continentTrue = sapply(conf4x$ancestryTrue, function(anc)
  continents_long[which(ancTypes == anc)])
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




