# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)

bees <- read.table("results/SNPs/thin1kb_common3/pass1_2018.ploidy", stringsAsFactors = F, 
                                     header = F, sep = "\t")$V1
prior <- c(0.4, 0.4, 0.2) # prior on global ancestry proportions

#dir_post = "results/ancestry_hmm/thin1kb_common3/pass1_2018_0.4_0.4_0.2/fixed_t_60_30"
#dir_post = "results/ancestry_hmm/thin1kb_common3/pass1_2018_0.05_0.05_0.9"
dir_post = "results/ancestry_hmm/thin1kb_common3/pass1_2018_0.4_0.4_0.2/reverse_order_CMA"
genotypes = c("CC", "CM", "AC", "MM", "AM", "AA") # used in calc_anc_from_post

# helper function reads posterior files
read_post <- function(id, dir = dir_post){
  post <- read.table(paste0(dir, "/", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., CC = X2.0.0) %>%
    rename(., CM = X1.1.0) %>%
    rename(., AC = X1.0.1) %>%
    rename(., MM = X0.2.0) %>%
    rename(., AM = X0.1.1) %>%
    rename(., AA = X0.0.2) %>%
    mutate(., Bee_ID = id)
}
# get posterior for individuals
get_post_simple = function(id, dir = dir_post){ 
  post <- read_post(id, dir)
  post$max_p = apply(post[ , genotypes], 1, max) # which ancestry state has the highest probability?
  small1 <- tidyr::gather(post, "anc", "p", genotypes) %>%
    filter(., p == max_p) %>% # only keep highest posterior prob. ancestry; and filter to one in every nSkip + 1
    arrange(., chrom, position)
  return(small1)
}


# ancestry proportions
calc_anc_from_post <- function(post, genos = genotypes){
  anc1 <- c(1, .5, .5, 0, 0, 0)
  anc2 <- c(0, .5, 0, 1, .5, 0)
  anc3 <- c(0, 0, .5, 0, .5, 1)
  anc_matrix <- cbind(anc1, anc2, anc3)
  anc <- as.data.frame(t(apply(post[ , genos], 1, function(row) row %*% anc_matrix)), stringsAsFactors = F)
  colnames(anc) <- c("C", "M", "A")
  return(anc)
}

# get ancestry proportions across the genome
# ancestry proportions
get_anc <- function(id, dir = dir_post){
  post <- read_post(id, dir)
  anc <- calc_anc_from_post(post)
  return(bind_cols(dplyr::select(post, c("chrom", "position", "Bee_ID")), anc))
}

# moderate A ancestry
postCA0401 = read_post("CA0401")
postSmallCA0401 = get_post_simple(id = "CA0401")
apply(calc_anc_from_post(postCA0401), 2, mean)

# high A ancestry
postAR2704 = read_post("AR2704")
apply(calc_anc_from_post(postAR2704), 2, mean)

# low A ancestry
postCA1410 = read_post("CA1410")
postSmallCA1410 = get_post_simple(id = "CA1410")
table(postSmallCA1410$anc)/nrow(postSmallCA1410)
apply(calc_anc_from_post(postCA1410), 2, mean)

for (i in c("CA0401", "AR2711", "CA1410", "AR0302")){
  posti <- read.table(paste0(dir_post, "/", i, ".posterior"), 
                      stringsAsFactors = F, header = T)
  apply(calc_anc_from_post(posti[ , 3:8]), 2, mean)
}
bee_subset <- c("CA0401", "AR2711", "CA1410", "AR0302")

# for each bee, calculate mean ancestry
mean_anc_all <- as.data.frame(t(sapply(bees, function(id)
  apply(calc_anc_from_post( # get mean ancestry from posterior
    read.table(paste0(dir_post, "/", id, ".posterior"), 
               stringsAsFactors = F, header = T)[ , 3:8]), 
    2, mean))))
colnames(mean_anc_all) <- c("C", "M", "A")
mean_anc_all$Bee_ID <- rownames(mean_anc_all)

# plot ancestry one way
mean_anc_all %>%
  gather(., "anc", "mean", c("A", "C", "M")) %>%
  ggplot(., aes(x = Bee_ID, y = mean)) +
  geom_point(aes(color = anc))
# plot ancestry similar to NGSAdmix results to compare:
mean_anc_all %>%
  gather(., "anc", "mean", c("A", "C", "M")) %>%
  ggplot(., aes(fill=anc, y=mean, x=Bee_ID)) +
  geom_bar(stat = "identity", position = "fill") + facet_wrap(~substr(Bee_ID, 1, 2))



# get all the posterior data
post_all = do.call(rbind,
                lapply(bees, function(i) 
  get_post_simple(id = i)))

# plot all ind's ancestry (with uncertainty) over some region
post_all %>%
  filter(chrom == "Group1") %>%
  ggplot(aes(x=position, y = Bee_ID)) +
  geom_point(aes(color = anc, alpha = p), size = .5)

# plot all ind's ancestry for the whole genome, chromosome by chromosome
for (chr_i in paste0("Group", 1:16)){
  p_chrom <- post_all %>%
    filter(chrom == chr_i) %>%
    ggplot(aes(x=position, y = Bee_ID)) +
    geom_point(aes(color = anc, alpha = p), size = .5) + 
    ggtitle(paste0("Ancestry 2018 bees -- ", chr_i))
  ggsave(paste0("plots/local_ancestry_tracts_pass1_2018_", chrom, ".png"), 
         plot = p_chrom,
         device = "png", 
         width = 20, height = 16, units = "in",
         dpi = 200)
}
table(small_all$chrom)
p_by_chrom <- small_all %>%
  ggplot(aes(x=position, y = Bee_ID)) +
  geom_point(aes(color = anc, alpha = p), size = .5) +
  facet_wrap(~chrom) +
  ggtitle(paste0("Ancestry 2018 bees by chromosome"))
plot(p_by_chrom)
ggsave("plots/local_ancestry_tracts_pass1_2018_by_chrom.png", 
       plot = p_by_chrom,
       device = "png", 
       width = 20, height = 16, units = "in",
       dpi = 200)
  
# what confidence does the HMM have in the calls it makes?
hist(post_all$p)
summary(post_all$p)
# does the confidence vary by ancestry?
post_all %>%
  group_by(., anc) %>%
  summarize(mean_p = mean(p))
# very confident in all calls, but especially homozygous ancestry calls.

# now compare tracts from different runs of ancestry_hmm for a subset of bees:
# how sensitive is the result to different times of admixture (t)?
est_ts <- c("fix_t_100_60", "fix_t_60_30", "est_t2")
dirs <- paste0("results/ancestry_hmm/thin1kb_common3/pass1_2018_0.4_0.4_0.2_scaffolds/reverse_order_CMA_",
               est_ts)
post_subset = lapply(dirs, function(x)
  do.call(rbind,
                   lapply(bee_subset, function(b) 
                     get_post_simple(id = b, dir = x))))
str(post_subset)
anc_subset = post_subset[[3]] %>%
  rename(anc_35.18 = anc) %>%
  mutate(anc_60.30 = post_subset[[2]]$anc) %>%
  mutate(anc_100.60 = post_subset[[1]]$anc)
# does setting different admixture timing change estimates for mean A ancestry?
mean_anc_subset <- lapply(bee_subset, function(b)
          do.call(rbind,
                  lapply(dirs, function(x) 
            apply(calc_anc_from_post( # get mean ancestry from posterior
              read.table(paste0(x, "/", b, ".posterior"), 
                         stringsAsFactors = F, header = T)[ , 3:8]), 2, mean)
          )))
for (i in 1:length(bee_subset)){
  colnames(mean_anc_subset[[i]]) <- c("C", "M", "A")
  rownames(mean_anc_subset[[i]]) <- est_ts
}
# no more than 1-2% difference, but larger times lead to a little more minor ancestry


# what is the rate of different top ancestry calls? generally around 1% ish
anc_subset[, c("anc_35.18", "anc_60.30", "Bee_ID")] %>%
  table(.)

# plot ancestry tracts for a few segments:
anc_subset %>%
  gather(., "est_t", "top_anc", c("anc_35.18", "anc_60.30", "anc_100.60")) %>%
  filter(chrom == "Group10.1") %>%
  filter(position < 1000000) %>%
  ggplot(aes(x = position, y = est_t, alpha = max_p, color = top_anc)) +
  geom_point() +
  facet_wrap(~Bee_ID) +
  ggtitle(paste0("Ancestry tracts sensitivity to time of admixture LG10.1"))
ggsave("plots/example_4_bees_LG10.1_ancestry_calls_diff_t_est.png", 
       device = "png", 
       width = 10, height = 8, units = "in",
       dpi = 200)
# plot second example set of tracts - some additional short tracts with longer times but no big differences
anc_subset %>%
  gather(., "est_t", "top_anc", c("anc_35.18", "anc_60.30", "anc_100.60")) %>%
  filter(chrom == "Group1.5") %>%
  filter(position < 1000000) %>%
  ggplot(aes(x = position, y = est_t, alpha = max_p, color = top_anc)) +
  geom_point() +
  facet_wrap(~Bee_ID) +
  ggtitle(paste0("Ancestry tracts sensitivity to time of admixture LG1.5"))
ggsave("plots/example_4_bees_LG1.5_ancestry_calls_diff_t_est.png", 
       device = "png", 
       width = 10, height = 8, units = "in",
       dpi = 200)
# look at spot where Nelson 2017 paper found high european ancestry
anc_subset %>%
  gather(., "est_t", "top_anc", c("anc_35.18", "anc_60.30", "anc_100.60")) %>%
  filter(chrom == "Group1.5") %>%
  filter(position < 1000000) %>%
  ggplot(aes(x = position, y = est_t, alpha = max_p, color = top_anc)) +
  geom_point() +
  facet_wrap(~Bee_ID) +
  ggtitle(paste0("Ancestry tracts sensitivity to time of admixture LG1.5"))
ggsave("plots/example_4_bees_LG1.5_ancestry_calls_diff_t_est.png", 
       device = "png", 
       width = 10, height = 8, units = "in",
       dpi = 200)

head(post_subset[[1]])


# is ancestry roughly in HW? I don't necessarily expect that it will be, but probably not way off for most individuals..
print("Observed ancestry homozygosity:")
table(filter(post_subset[[1]], anc %in% c("CC", "MM", "AA"))[ , c("anc", "Bee_ID")])/(nrow(post_subset[[1]])/4) 
# expected
print("HW expectation:")
apply(table(post_subset[[1]][ , c("anc", "Bee_ID")])/(nrow(post_subset[[1]])/4), 
      2, 
      function(x) c(C = x["CC"] + .5*(x["CM"] + x["AC"]), 
                    M = x["MM"] + .5*(x["CM"] + x["AM"]),
                    A = x["AA"] + .5*(x["AC"] + x["AM"])))^2
# the ancestry calls are pretty close to HW expectation; seems good
# closest to HW for CA0401 which is close to the ancestry prior
# For some other bees with mean ancestry further from the prior, 
# the ancestry homozygosity is higher than expected from HW

# what do ancestry blocks look like?


# are there clear peaks of high (or low) African ancestry?
str(post_subset[[1]])
# what is the mean A ancestry?
getAncFreq = function(id, dir = dir_post){ 
  post <- read.table(paste0(dir, "/", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., CC = X2.0.0) %>%
    rename(., CM = X1.1.0) %>%
    rename(., AC = X1.0.1) %>%
    rename(., MM = X0.2.0) %>%
    rename(., AA = X0.0.2) %>%
    rename(., AM = X0.1.1)
  anc_matrix = cbind(C = c(1, .5, .5, 0, 0, 0),
                     M = c(0, .5, 0, 1, .5, 0),
                     A = c(0, 0, .5, 0, .5, 1))
  anc <- as.data.frame(t(apply(post[ , c("CC", "CM", "AC", "MM", "AA", "AM")], 1, 
               function(row) row %*% anc_matrix)))
  colnames(anc) <- c("C", "M", "A")
  return(anc)
}

test <- getAncFreq(id = "CA0401", dir = dirs[1])
anc_all <- lapply(bees, function(b) getAncFreq(id = b, dir = dirs[1]))
# what is the ancestry in the total population at each site?
anc_total <- Reduce('+', anc_all)/length(anc_all) # sum all C's, M's, and A's contribution across individuals
apply(anc_total, 2, hist) # make histogram of the distribution of % C/M/A across loci
summary(anc_total)
# split it into Argentina vs. California samples -- maybe not a lot of selection is shared?
listy = list(1:5, 1:3, 3:4)
anc_CA <- Reduce('+', anc_all[substr(bees, 1, 2) == "CA"])/sum(substr(bees, 1, 2) == "CA")
anc_AR <- Reduce('+', anc_all[substr(bees, 1, 2) == "AR"])/sum(substr(bees, 1, 2) == "AR")
summary(anc_CA)
summary(anc_AR)
plot(anc_CA$A ~ anc_AR$A)
summary(lm(anc_CA$A ~ anc_AR$A)) # not very correlated in general. but what about across diff. frequency bins?
summary(lm(anc_CA$A[anc_CA$A > .5 | anc_AR$A > .5] ~ anc_AR$A[anc_CA$A > .5 | anc_AR$A > .5]))
summary(lm(anc_CA$A[anc_CA$A < .3 | anc_AR$A < .3] ~ anc_AR$A[anc_CA$A < .3 | anc_AR$A < .3]))
# what if I divide CA and AR bees into two groups? Are CA bees more correlated with each other that with AR bees?
mean_anc_all$location <- substr(mean_anc_all$Bee_ID, 1, 2)
mean_anc_all$popN <- as.integer(substr(mean_anc_all$Bee_ID, 3, 4))

anc_CA1 <- Reduce('+', anc_all[mean_anc_all$location == "CA"][1:12])/12
anc_CA2 <- Reduce('+', anc_all[mean_anc_all$location == "CA"][13:24])/12
anc_AR1 <- Reduce('+', anc_all[mean_anc_all$location == "AR"][1:20])/20
anc_AR2 <- Reduce('+', anc_all[mean_anc_all$location == "AR"][21:39])/19
plot(anc_CA1$A ~ anc_CA2$A)
summary(lm(anc_CA1$A ~ anc_CA2$A))
plot(anc_AR1$A ~ anc_AR2$A)
summary(lm(anc_AR1$A ~ anc_AR2$A)) # less correlated than the CA bees, but perhaps unfair to include whole zone
anc_AR1 <- Reduce('+', anc_all[substr(bees, 1, 2) == "AR"][1:20])/20



# I should plot clines across CA and Argentina for allele frequencies. 
# Overall freq may be less impressive than cline shape.

# ancestry estimates for all bees at all positions:
anc_all = do.call(rbind,
        lapply(bees, function(i) get_anc(id = i, dir = dir_post)))
head(anc_all)
anc_all_plus_meta <- anc_all %>%
  left_join(., dplyr::select(meta, c("Bee_ID", "geographic_location", "population", "group", "lat", "long", "popN", "indN", "enjambre")),
            by = "Bee_ID")
anc_pops <- anc_all_plus_meta %>% 
  group_by(group, chrom, position) %>%
  summarise(pop_A = mean(A))
anc_pops %>%
  group_by(group) %>%
  summarise(max(pop_A))
anc_pops %>%
  group_by(group) %>%
  summarise(min(pop_A))
anc_pops %>%
  ggplot(aes(x = position, y = pop_A, color = group)) +
  geom_point(size = .1) +
  facet_wrap(~chrom)
# look at euro ancestry peak region from Nelson 2017:
# looks like the peak in European ancestry is present in Argentinian pops 
# (lowest A ancestry region in that data)
# but not starkly in California pops
# though both CA has a more localized dip right before 13Mb, 
# that doesn't quite line up with their candidate gene but could indicate a diff. candidate
anc_pops %>%
  filter(chrom == "Group11") %>%
  filter(position > 10000000) %>%
  filter(position < 15000000) %>%
  ggplot(aes(x = position, y = pop_A, color = group)) +
  geom_point()
# correlation is fairly low between CA and Argentina -- what is variation within?
cov(filter(anc_pops, group == "AR_2018")$pop_A, filter(anc_pops, group == "CA_2018")$pop_A)
var(filter(anc_pops, group == "AR_2018")$pop_A)
var(filter(anc_pops, group == "CA_2018")$pop_A) # somewhat higher in CA than AR, though not huge

anc_all_plus_meta %>%
  filter(chrom == "Group11") %>%
  filter(position > 10000000) %>%
  filter(position < 15000000) %>%
  ggplot(aes(x = position, y = Bee_ID, color = A)) +
  geom_point(size = 1) + 
  facet_wrap(~group) +
  ggtitle("high EU region from Nelson 2017")

anc_all_plus_meta %>%
  filter(chrom == "Group1") %>%
  filter(position > 20000000) %>%
  filter(position < 25000000) %>%
  ggplot(aes(x = position, y = Bee_ID, color = A)) +
  geom_point(size = 1) +
  facet_wrap(~group) +
  ggtitle("random region LG1")

# plot a frequency across latitude
anc_all_plus_meta %>%
  filter(chrom == "Group11") %>%
  filter(position > 13000000) %>%
  filter(position < 13100000) %>%
  group_by(popN, group) %>%
  summarise(locus_A = mean(A),
            mean_lat = abs(mean(lat))) %>%
  ggplot(aes(x = mean_lat, y = locus_A, color = group)) +
  geom_point()
