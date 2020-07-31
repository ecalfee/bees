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
#load("results/d_1-10x_high_CMA.RData")

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
  geom_line(alpha = .8) +
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
  ylab("A ancestry correlation") +
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
       units = "in", dpi = 600, device = "tiff",
       compression = "lzw", type = "cairo")



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
