# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(reshape2)
library(viridis)
library(LaplacesDemon)
library(emdbook)
library(betareg)
library(gridExtra)
library(MASS) # for mvrnorm

source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions


bees1 <- read.table("results/SNPs/thin1kb_common3/pass1_2018.ploidy", stringsAsFactors = F, 
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
mean_anc_all <- as.data.frame(t(sapply(bees1, function(id)
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
                lapply(bees1, function(i) 
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
    ggtitle(paste0("Ancestry 2018 bees1 -- ", chr_i))
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
  ggtitle(paste0("Ancestry 2018 bees1 by chromosome"))
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

# now compare tracts from different runs of ancestry_hmm for a subset of bees1:
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
anc_all <- lapply(bees1, function(b) getAncFreq(id = b, dir = dirs[1]))
# what is the ancestry in the total population at each site?
anc_total <- Reduce('+', anc_all)/length(anc_all) # sum all C's, M's, and A's contribution across individuals
apply(anc_total, 2, hist) # make histogram of the distribution of % C/M/A across loci
summary(anc_total)
# split it into Argentina vs. California samples -- maybe not a lot of selection is shared?
listy = list(1:5, 1:3, 3:4)
anc_CA <- Reduce('+', anc_all[substr(bees1, 1, 2) == "CA"])/sum(substr(bees1, 1, 2) == "CA")
anc_AR <- Reduce('+', anc_all[substr(bees1, 1, 2) == "AR"])/sum(substr(bees1, 1, 2) == "AR")
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
        lapply(bees1, function(i) get_anc(id = i, dir = dir_post)))
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

#************************************************************************************************************#

# Newest bee sequences
# load metadata
# load population ancestry frequencies:
pops <- read.table("../bee_samples_listed/byPop/pops_included.list", stringsAsFactors = F)$V1
bees <- do.call(rbind, 
                lapply(pops, function(p) data.frame(Bee_ID = read.table(paste0("../bee_samples_listed/byPop/", p, ".list"),
                                                    stringsAsFactors = F)$V1, population = p, stringsAsFactors = F)))

# get metadata
meta.ind <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  left_join(bees, ., by = c("Bee_ID", "population")) 
meta.pop <- meta.ind %>%
  dplyr::select(c("population", "source", "year", "group", "lat", "long")) %>%
  dplyr::group_by(population, source, year, group) %>%
  dplyr::summarise(n_bees = n(),
                   lat = mean(lat),
                   long = mean(long)) %>%
  left_join(data.frame(population = pops, stringsAsFactors = F),
            ., by = "population") %>%
  mutate(lat = ifelse(population == "Riverside_1999", .[.$population == "Riverside_2014", "lat"], lat)) %>%
  mutate(zone = ifelse(group == "AR_2018", "S. America", "N. America"))


# get ancestry frequencies for each population across the genome
dir_results <- "results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot"

popA <- lapply(pops, function(p) read.table(paste0(dir_results, "/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A <- do.call(cbind, popA)
colnames(A) <- pops
popM <- lapply(pops, function(p) read.table(paste0(dir_results, "/anc/", p, ".M.anc"),
                                            stringsAsFactors = F))
M <- do.call(cbind, popM)
colnames(M) <- pops
popC <- lapply(pops, function(p) read.table(paste0(dir_results, "/anc/", p, ".C.anc"),
                                            stringsAsFactors = F))
C <- do.call(cbind, popC)
colnames(C) <- pops

# mean ancestry across populations
meanA0 <- apply(A, 1, mean)
meanC0 <- apply(C, 1, mean)
meanM0 <- apply(M, 1, mean)

# mean ancestry across the genome for each pop
admix_proportions <- data.frame(population = pops,
                                A = apply(A, 2, mean),
                                C = apply(C, 2, mean),
                                M = apply(M, 2, mean),
                                stringsAsFactors = F)

# get time of admixture estimates
time_pops <- read.table(paste0(dir_results, "/", "time_pops.txt"), 
                        stringsAsFactors = F, header = F) %>%
  data.table::setnames("population")
admix_times <- do.call(rbind, lapply(c("A", "C", "M"), function(anc) read.table(paste0(dir_results, "/", "time_", anc, ".txt"), 
                     stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("ancestry_n", "time", "proportion")) %>%
    mutate(ancestry = anc) %>% 
    cbind(time_pops, .)))

# the ancestry proportion values in the log rom ancestry_hmm are just the priors from NGSAdmix:
# d_ACM and d_pop_ACM are loaded from plot_NGSadmix.R script
left_join(d_pop_ACM, filter(admix_times, ancestry == "A"), by = "population") %>%
  with(., plot(A, proportion, xlim = 0:1, ylim = 0:1, col = "blue"))


# plot ancestry times vs. mean ancestry proportion:
admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  filter(., ancestry != "C") %>% # no time, first ancestry
  ggplot(., aes(x = abs(lat), y = time, color = zone, shape = factor(year))) +
  geom_point() +
  facet_grid(. ~ ancestry) +
  ggtitle("Inferred time of admixture pulses from HMM") +
  ylab("Time (generations)") +
  xlab("Degrees latitude from the equator") +
  labs(color = "Hybrid Zone", shape = "Collection")
ggsave("plots/time_of_admixture_vs_latitude.png",
       height = 5, width = 14, units = "in")
ggsave("../../bee_manuscript/figures/time_of_admixture_vs_latitude.png",
       height = 5, width = 14, units = "in")

# compare California and Argentina: 
# For the same ancestry proportions, do they have similar inferred times of admixture?
# This would be evidence for different amounts of ongoing gene flow in the two zones.
admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  #filter(., ancestry == "A") %>%
  filter(., ancestry != "C") %>%
  mutate(ancestry_labels = paste(ancestry, "ancestry admixture pulse")) %>%
  ggplot(., aes(x = proportion, y = time, color = zone, shape = factor(year))) +
  geom_point() +
  ggtitle("Time since admixture estimated from ancestry block lengths") +
  facet_grid(. ~ ancestry_labels) +
  ylab("Time (generations in the past)") +
  xlab("Admixture proportion (ancestry_hmm)") +
  labs(color = "Hybrid Zone", shape = "Collection")
ggsave("plots/California_has_shorter_ancestry_blocks.png",
       height = 3, width = 8, units = "in")
ggsave("../../bee_manuscript/figures/California_has_shorter_ancestry_blocks.png",
       height = 3, width = 8, units = "in")


sites0 <- read.table("results/SNPs/thin1kb_common3/included_scaffolds.pos", stringsAsFactors = F,
                    sep = "\t", header = F)
colnames(sites0) <- c("scaffold", "pos")
sites <- tidyr::separate(sites0, scaffold, c("chr", "scaffold_n"), remove = F) %>%
  mutate(scaffold_n = as.numeric(scaffold_n)) %>%
  mutate(chr_n = as.numeric(substr(chr, 6, 100))) %>%
  mutate(snp_id = paste0("snp", chr_n, ".", scaffold, ".", pos))
#test_id <- read.table("results/SNPs/thin1kb_common3/included.snplist", stringsAsFactors = F,
#                      sep = "\t", header = F)$V1
#table(test_id == sites$snp_id) # good

# get individual ancestries for just 
# CA_2018 and AR_2018 samples,
# same # individuals per pop (=8)
exclude_over8 <- c("AR0115", "AR0501", "AR0818",
                   "AR1019", "AR1401", "AR1604", 
                   "AR2002", "AR2303", "AR2608",
                   "AR2912", "CA0201", "CA0502",
                   "CA1010", "CA1302")
meta.ind %>% 
  filter(!(Bee_ID %in% exclude_over8)) %>%
  group_by(population) %>%
  summarise(n = n()) %>%
  filter(n != 8)
# get ID's for 8 bees per CA and AR pop (2018)
# divide into subgroups:
bees_all8_0 <- meta.ind %>%
  filter(!(Bee_ID %in% exclude_over8)) %>%
  filter(group %in% c("CA_2018", "AR_2018")) %>%
  #left_join(., subgroups, by = "population") %>%
  arrange(lat)
subgroups = data.frame(population = unique(bees_all8_0$population),
                       subgroup = c(rep("AR_S", 6),
                                    rep("AR_mid", 8),
                                    rep("AR_N", 7),
                                    rep("CA_S", 6),
                                    rep("CA_N", 6)),
                       stringsAsFactors = F)
bees_all8 <- bees_all8_0 %>%
  left_join(., subgroups, by = "population")

# get ancestries
A_bees_all8 <- lapply(bees_all8$Bee_ID, function(p) read.table(paste0("results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A_all8 <- do.call(cbind, A_bees_all8)
colnames(A_all8) <- bees_all8$Bee_ID

# make plots - all loci, mean ancestry across all pops
png("plots/histogram_A_ancestry_all_loci.png", height = 6, width = 8, units = "in", res = 300)
hist(meanA, main = "all loci: mean ancestry across populations")
abline(v = max(meanA), col = "blue")
abline(v = min(meanA), col = "blue")
abline(v = quantile(meanA, .99), col = "orange")
abline(v = quantile(meanA, .01), col = "orange")
dev.off()

CA_pops_included <- meta.pop[meta.pop$group %in% c("N_CA", "S_CA", "CA_2018") & meta.pop$year >= 2014, ]
CA_A <- A[, CA_pops_included$population]
CA_A_earlier <- A[, meta.pop$population[meta.pop$group %in% c("N_CA", "S_CA") & meta.pop$year < 2014]]
CA_A_2014 <- A[, meta.pop$population[meta.pop$group %in% c("N_CA", "S_CA") & meta.pop$year == 2014]]
CA_A_2018 <- A[, meta.pop$population[meta.pop$group %in% c("CA_2018")]]
plot(apply(CA_A_2014, 1, mean), apply(CA_A_2018, 1, mean))
plot(apply(CA_A_earlier, 1, mean), apply(CA_A_2018, 1, mean))
AR_pops_included <- meta.pop[meta.pop$group == "AR_2018", ]
AR_A <- A[, AR_pops_included$population]
# take mean across individuals, not across populations:
meanA_CA0 <- apply(CA_A, 1, mean)
meanA_CA <- apply(CA_A, 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanA_AR0 <- apply(AR_A, 1, mean)
meanA_AR <- apply(AR_A, 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))

# with only included pops, mean of all individuals:
meanA <- apply(A[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(CA_pops_included$n_bees, AR_pops_included$n_bees)))


sites %>%
  filter(meanA_CA > quantile(meanA_CA, .99) & 
            meanA_AR > quantile(meanA_AR, .99)) %>%
  dplyr::select(scaffold) %>%
  table() # outliers are on just a few scaffolds
# plot
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter(meanA_CA > quantile(meanA_CA, .99) & 
           meanA_AR > quantile(meanA_AR, .99)) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) # really just a few peaks
sites %>%
  filter(meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  dplyr::select(scaffold) %>%
  table() # low A outliers still only on 14 scaffolds
#plot
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter(meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) # really just a few peaks
sites %>%
  filter(meanA_CA > quantile(meanA_CA, .25) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  dplyr::select(scaffold) %>%
  table() # several regions, but at least one big hit on chr 11 - Group11.18
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>% # plot
  filter(meanA_CA > quantile(meanA_CA, .25) & 
           meanA_AR < quantile(meanA_AR, .01)) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) # several regions, but at least one big hit on chr 11 - Group11.18
sites %>%
  filter(meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR > quantile(meanA_AR, .25)) %>%
  dplyr::select(scaffold) %>%
  table() # more spread out but a few scaffolds with higher concentrations
sites %>%
  filter(meanA_CA < quantile(meanA_CA, .75) & 
           meanA_AR > quantile(meanA_AR, .99)) %>%
  dplyr::select(scaffold) %>%
  table()
sites %>%
  filter(meanA_CA > quantile(meanA_CA, .99) & 
           meanA_AR < quantile(meanA_AR, .75)) %>%
  dplyr::select(scaffold) %>%
  table()
# plot shared high and shared low on same axes:
#plot
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter((meanA_CA < quantile(meanA_CA, .01) & 
           meanA_AR < quantile(meanA_AR, .01)) |
           (meanA_CA > quantile(meanA_CA, .99) & 
              meanA_AR > quantile(meanA_AR, .99))) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point() +
  facet_wrap(~scaffold) + # really just a few peaks
  ggtitle("shared peaks high and low A ancestry")


# plot all data. oops I need to use position relative to the chromosome instead.
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>% # plot
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point(size = .1) +
  facet_wrap(~chr, scales = "free_x")

png("plots/scatterplot_CA_vs_AR_A_ancestry_all_loci.png", height = 6, width = 8, units = "in", res = 300)
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(meanA_CA, .99) & 
                                                   meanA_AR > quantile(meanA_AR, .99)) |
                                                  (meanA_CA < quantile(meanA_CA, .01) & 
                                                     meanA_AR < quantile(meanA_AR, .01)), 
          alpha("orange", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, orange = top 1% shared outliers")
abline(lm(meanA_AR~meanA_CA), col = "blue")
dev.off()
cor(meanA_CA, meanA_AR)

# color outliers by scaffold/chromosome:
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>% # plot
  ggplot(aes(x = CA, y = AR, color = scaffold_n)) +
  geom_point(size = .1) +
  facet_wrap(~chr) + 
  ggtitle("African ancestry frequencies in CA vs. AR, all loci")
ggsave("plots/mean_A_ancestry_CA_vs_AR_rainbow_by_chr.png", 
       device = "png", height = 10, width = 12, units = "in")
# zoom in on Group1.23 because it has both high and low outliers:
data.frame(sites, CA = meanA_CA, AR = meanA_AR) %>%
  filter(chr == "Group1" & scaffold_n %in% 20:23) %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point(size = .1) +
  facet_wrap(~scaffold_n) + 
  ggtitle("African ancestry frequencies in CA and AR")
ggsave("plots/mean_A_ancestry_CA_and_AR_Group1_scaffolds_20-23.png", 
       device = "png", height = 10, width = 12, units = "in")


# plot K matrix
zAnc_bees = make_K_calcs(t(A))


melt(zAnc_bees$K) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K covariance - bees")
ggsave("plots/k_matrix_all_pops.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
melt(cov2cor(zAnc_bees$K)) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K correlation - bees")
ggsave("plots/k_correlation_matrix_all_pops.png", 
       height = 6, width = 8, 
       units = "in", device = "png")


summary(meanA)
sd(meanA)

# simulate from MVN distribution:
AR_CA_K <- make_K_calcs(t(cbind(AR_A, CA_A)))
#sim_n <- 10^6 # slow
sim_n <- 10^5
set.seed(101)
MVNsim <- mvrnorm(n = sim_n, 
                  mu = AR_CA_K$alpha, 
                  Sigma = AR_CA_K$K, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)



# sets bounds at 0 and 1
MVNsim_bounded <- MVNsim 
MVNsim_bounded[MVNsim < 0] <- 0
MVNsim_bounded[MVNsim > 1] <- 1

MVNsim_AR <- MVNsim[ , AR_pops_included$population]
MVNsim_AR_bounded <- MVNsim_bounded[ , AR_pops_included$population]
MVNsim_CA <- MVNsim[ , CA_pops_included$population]
MVNsim_CA_bounded <- MVNsim_bounded[ , CA_pops_included$population]
# mean across populations
meanA_MVNsim_AR_bounded0 <- apply(MVNsim_AR_bounded, 1, mean)
meanA_MVNsim_CA_bounded0 <- apply(MVNsim_CA_bounded, 1, mean)
# mean across individuals
meanA_MVNsim_AR_bounded <- apply(MVNsim_AR_bounded[ , AR_pops_included$population], 1, 
                                  function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_bounded <- apply(MVNsim_CA_bounded[ , CA_pops_included$population], 1, 
                                  function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
# combined mean across individuals
meanA_MVNsim_bounded <- apply(cbind(MVNsim_AR_bounded, MVNsim_CA_bounded)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))
meanA_MVNsim <- apply(cbind(MVNsim_AR, MVNsim_CA)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))

# Add in constraint that all covariances must be 0 or positive
# (effectively I zero out negative covariances):
# simulate from MVN distribution:
AR_CA_K_zero <- AR_CA_K$K
summary(AR_CA_K_zero[AR_CA_K$K < 0])
sum(AR_CA_K$K < 0)/prod(dim(AR_CA_K$K)) # 11% of covariances < 0
summary(AR_CA_K$K[AR_CA_K$K > 0])
summary(AR_CA_K$K[T])
AR_CA_K_zero[AR_CA_K$K < 0] <- 0
set.seed(101)
MVNsim_zero <- mvrnorm(n = sim_n, 
                  mu = AR_CA_K$alpha, 
                  Sigma = AR_CA_K_zero, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)

# how many simulations exceed the bounds?
# overall few, but more than 15% for some populations
sum(MVNsim_zero < 0)/sum(table(MVNsim_zero<0))
sum(MVNsim_zero > 1)/sum(table(MVNsim_zero>1))
hist(MVNsim[MVNsim_zero<0])
table(MVNsim_bounded[MVNsim_zero < 1])
apply(MVNsim, 2, function(x) sum(x > 1))/(nrow(MVNsim_zero))
summary(apply(MVNsim, 2, function(x) sum(x <0))/(nrow(MVNsim_zero)))
data.frame(population = colnames(MVNsim_zero),
           percent_0 = apply(MVNsim_zero, 2, function(x) sum(x < 0))/nrow(MVNsim_zero),
           stringsAsFactors = F) %>%
  left_join(., admix_proportions, by = "population") %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = A, y = percent_0, color = zone)) +
  geom_point() +
  xlab("Mean A ancestry proportion (ancestry_hmm)") +
  ylab("Percent sims < 0 (MVN)") +
  ggtitle("Percent simulated population ancestry frequencies < 0 (MVN)")
ggsave("plots/percents_MVN_sims_need_truncation_low.png",
       height = 3, width = 8, units = "in")



MVNsim_zero_bounded <- MVNsim_zero # sets bounds at 0 and 1
MVNsim_zero_bounded[MVNsim_zero < 0] <- 0
MVNsim_zero_bounded[MVNsim_zero > 1] <- 1
MVNsim_AR_zero <- MVNsim_zero[ , AR_pops_included$population]
MVNsim_AR_zero_bounded <- MVNsim_zero_bounded[ , AR_pops_included$population]
MVNsim_CA_zero <- MVNsim_zero[ , CA_pops_included$population]
MVNsim_CA_zero_bounded <- MVNsim_zero_bounded[ , CA_pops_included$population]
# mean across individuals
meanA_MVNsim_AR_zero_bounded <- apply(MVNsim_AR_zero_bounded[ , AR_pops_included$population], 1, 
                                 function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_zero_bounded <- apply(MVNsim_CA_zero_bounded[ , CA_pops_included$population], 1, 
                                 function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
# combined mean across individuals
meanA_MVNsim_zero_bounded <- apply(cbind(MVNsim_AR_zero_bounded, MVNsim_CA_zero_bounded)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                              function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))
meanA_MVNsim_zero <- apply(cbind(MVNsim_AR_zero, MVNsim_CA_zero)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                      function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))




#summarise
summary(meanA_MVNsim)
summary(meanA_MVNsim_bounded) # makes little difference
summary(meanA_MVNsim_zero)
summary(meanA_MVNsim_zero_bounded)
# plots
hist(apply(MVNsim, 1, mean))
hist(meanA_MVNsim)
hist(meanA_MVNsim_bounded)
plot(meanA_MVNsim_AR_bounded, meanA_MVNsim_CA_bounded, xlim = c(0, 1), ylim = c(0, 1),
     main = "simulation CA vs. AR frequencies (zero-d K matrix MVN sim in blue)")
points(meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_CA_zero_bounded, col = "blue")
plot(meanA_MVNsim_bounded, meanA_MVNsim_zero_bounded, col = "orange", 
     xlim = c(0.1, 0.7), ylim = c(0.1, 0.7),
     main = "effect of zero-ing out small negative covariances")
abline(a = 0, b = 1)
var(meanA_MVNsim_bounded)
var(meanA_MVNsim_zero_bounded) # very slightly higher variance -- could just be simulation noise
# this makes sense because the covariances I zero-d out were very slight
summary(diag(AR_CA_K$K))
summary(AR_CA_K$K[lower.tri(AR_CA_K$K, diag = F)])

#QQ-plots:
png("plots/QQ_plots_against_MVN.png",
    height = 10, width = 15, res = 300, units = "in")
par(mfrow=c(2,2))
qqplot(meanA_MVNsim_bounded, meanA_MVNsim,
       main = "Effect of 0-1 bounds on MVN distribution")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_bounded, meanA,
       main = "QQ-plot: MVN dist. vs. data, combined mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_CA_bounded, meanA_CA,
       main = "QQ-plot: MVN dist. vs. data, CA mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_AR_bounded, meanA_AR,
       main = "QQ-plot: MVN dist. vs. data, AR mean")
abline(0, 1, col = "blue")
dev.off()

png("plots/QQ_plots_against_MVN_zero-ed_out_negK.png",
    height = 10, width = 15, res = 300, units = "in")
par(mfrow=c(2,2))
qqplot(meanA_MVNsim_zero_bounded, meanA_MVNsim_bounded,
       main = "Effect of non-neg K on MVN distribution")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_zero_bounded, meanA,
       main = "QQ-plot: MVN dist. non-neg K vs. data, combined mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_CA_bounded,
       main = "QQ-plot: Effect of non-neg K on MVN CA mean")
abline(0, 1, col = "blue")
qqplot(meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_AR_bounded,
       main = "QQ-plot: Effect of non-neg K on MVN AR mean")
abline(0, 1, col = "blue")
dev.off()


make_qqplot_lines <- function(x, legend = T){
  abline(0, 1, col = "black")
  abline(v = quantile(x, 
                      c(0.99, 0.95, 0.75, 0.5, 0.01, 0.05, 0.25)), 
         col = c("lightgreen", "salmon", "lightblue", "yellow", "lightgreen", "salmon", "lightblue"))
  if (legend) {legend("bottomright", 
         legend = c(0.99, 0.95, 0.75, 0.5), 
         title = "Quantile",
         lty = 1,
         lwd = 2,
         col = c("lightgreen", "salmon", "lightblue", "yellow"))
  }
}

# make QQ plots for paper:
for (path in c("plots/", 
               "../../bee_manuscript/figures/")){
  png(paste0(path, "QQ_plots_MVN_vs_data.png"),
      height = 8, width = 8, res = 300, units = "in")
  png("plots/QQ_plots_MVN_vs_data.png",
      height = 8, width = 8, res = 300, units = "in")
  # plot poisson binomial
  par(mfrow=c(2,2))
  # plot effect of truncation
  par(oma = c(4, 1, 1, 1))
  qqplot(meanA_MVNsim_zero, meanA_MVNsim_zero_bounded,
         main = "QQ plot - A Ancestry effect of truncation",
         xlab = "MVN before truncation [0, 1]",
         ylab = "MVN after truncation [0, 1]",
         col = "grey")
  make_qqplot_lines(meanA_MVNsim_zero_bounded, legend = F)
  qqplot(meanA_MVNsim_zero_bounded, meanA,
         main = "QQ plot fit - A Ancestry Combined Sample",
         xlab = "MVN simulation",
         ylab = "Observed",
         col = "grey")
  make_qqplot_lines(meanA_MVNsim_zero_bounded, legend = F)
  qqplot(meanA_MVNsim_CA_zero_bounded, meanA_CA,
         main = "QQ plot fit - A Ancestry N. America",
         xlab = "MVN simulation",
         ylab = "Observed",
         col = "grey")
  make_qqplot_lines(meanA_MVNsim_CA_zero_bounded, legend = F)
  qqplot(meanA_MVNsim_AR_zero_bounded, meanA_AR,
         main = "QQ plot fit - A Ancestry S. America",
         xlab = "MVN simulation",
         ylab = "Observed",
         col = "grey")
  make_qqplot_lines(meanA_MVNsim_AR_zero_bounded, legend = F)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",
         legend = c(0.99, 0.95, 0.75, 0.5), 
         title = "Quantile",
         lty = 1,
         lwd = 2,
         xpd = TRUE, horiz = TRUE, 
         inset = c(0, 0),
         bty = "n",
         col = c("lightgreen", "salmon", "lightblue", "yellow"))
  dev.off()
  
}
png("plots/mean_A_ancestry_CA_vs_AR_grey.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = alpha("grey", .1),
     main = "comparison of A ancestry in California vs. Argentina zone")
dev.off()

# for how long do my quantiles match up?
qs <- seq(0, 1, by = .01)
qs_meanA_MVNsim_bounded <- quantile(meanA_MVNsim_bounded, probs = qs) 
qs_meanA <- quantile(meanA, probs = qs) 
qs_meanA_AR <- quantile(meanA_AR, probs = qs) 
qs_meanA_CA <- quantile(meanA_CA, probs = qs) 
qs_meanA_MVNsim_CA_bounded <- quantile(meanA_MVNsim_CA_bounded, probs = qs)
qs_meanA_MVNsim_AR_bounded <- quantile(meanA_MVNsim_AR_bounded, probs = qs)
qqplot(meanA_MVNsim_bounded, meanA)
points(qs_meanA_MVNsim_bounded, qs_meanA, col = "blue")
tail(qs_meanA_MVNsim_bounded)
tail(qs_meanA)
head(qs_meanA_MVNsim_bounded)
head(qs_meanA)
lapply(list(qs_meanA_AR, qs_meanA_MVNsim_AR_bounded, 
            qs_meanA_CA, qs_meanA_MVNsim_CA_bounded, 
            qs_meanA, qs_meanA_MVNsim_bounded), tail)
lapply(list(qs_meanA_AR, qs_meanA_MVNsim_AR_bounded, 
            qs_meanA_CA, qs_meanA_MVNsim_CA_bounded, 
            qs_meanA, qs_meanA_MVNsim_bounded), head)
# probably just removing top 1% is fine, but I could be safe and remove top 3% outliers
high = "97%"
low = "3%"
is_outlier <- (meanA >= qs_meanA[high] | meanA_AR >=  qs_meanA_AR[high] | meanA_CA >= qs_meanA_CA[high] |
                 meanA <= qs_meanA[low] | meanA_AR <=  qs_meanA_AR[low] | meanA_CA <= qs_meanA_CA[low])
# what % are outliers high or low? About 12% of my data. I'll exclude this set for demographic inference
table(is_outlier)/length(is_outlier)


# plot K matrix with = numbers per pop (CA and AR 2018)
# and excluding outlier SNPs


# make K matrix
K_all8 <- make_K_calcs(t(A_all8))
plot(diag(K_all8$K)/K_all8$alpha, bees_all8$lat)
plot(diag(K_all8$K), K_all8$alpha)
plot(diag(K_all8$K)/(K_all8$alpha*(1-K_all8$alpha)), abs(bees_all8$lat),
     col = ifelse(bees_all8$geographic_location == "California", 
                  "blue", "red"))
melt(K_all8$K) %>%
  filter(Var1 != Var2) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("K covariance - 8 bees per pop")
ggsave("plots/k_matrix_all8_inds.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
# take out outliers
K_all8_noOut <- make_K_calcs(t(A_all8[!is_outlier,]))
melt(K_all8_noOut$K) %>%
  filter(Var1 != Var2) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("K covariance - 8 bees per pop. No outlier loci.")
ggsave("plots/k_matrix_all8_inds_noOutliers.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
# mean covariance by group:
melt(K_all8_noOut$K) %>%
  filter(Var1 != Var2) %>%
  left_join(., bees_all8[ , c("Bee_ID", "population", "subgroup")], by = c("Var1"="Bee_ID")) %>%
  left_join(., bees_all8[ , c("Bee_ID", "population", "subgroup")], by = c("Var2"="Bee_ID")) %>%
  mutate(comparison = paste(subgroup.x, subgroup.y, sep = ".")) %>%
  mutate(same_pop = population.x == population.y) %>%
  ggplot(data = ., aes(x = comparison, y = value, color = same_pop)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("K covariance - 8 bees per pop. No outlier loci.")
melt(K_all8_noOut$K) %>%
  filter(Var1 != Var2) %>%
  left_join(., bees_all8[ , c("Bee_ID", "population", "subgroup")], by = c("Var1"="Bee_ID")) %>%
  left_join(., bees_all8[ , c("Bee_ID", "population", "subgroup")], by = c("Var2"="Bee_ID")) %>%
  mutate(comparison = paste(subgroup.x, subgroup.y, sep = ".")) %>%
  group_by(comparison) %>%
  summarise(mean_cov = mean(value)) %>%
  ggplot(., aes(x = comparison, y = mean_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# plot 50%, 95%, 99% & .999% credible intervals for my empirical observations
LaplacesDemon::joint.pr.plot(meanA_CA, 
                             meanA_AR, # top 1% may not be strict enough to ID outliers
                             quantiles=c(0.5, 0.95, 0.99, .999))


png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_outliers_blue.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(apply(MVNsim_CA_bounded, 1, mean), .9999) | 
                                                   meanA_AR > quantile(apply(MVNsim_AR_bounded, 1, mean), .9999)) |
                                                  (meanA_CA < quantile(apply(MVNsim_CA_bounded, 1, mean), .0001) | 
                                                     meanA_AR < quantile(apply(MVNsim_AR_bounded, 1, mean), .0001)), 
                                                alpha("blue", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, blue = top 0.01% sim")
dev.off()


png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_results_orange.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = alpha("grey", .1),
     main = "comparison of A ancestry in California vs. Argentina zone, orange = 100k MVN simulations")
points(apply(MVNsim_CA_bounded, 1, mean), apply(MVNsim_AR_bounded, 1, mean), pch = 20, col = alpha("orange", .02))
dev.off()
#points(apply(MVNsim_CA, 1, mean), apply(MVNsim_AR, 1, mean), pch = 20, col = alpha("green", .1))

# plot with both outliers in blue and MVN sims in orange
png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_results_orange_outliers_blue.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(apply(MVNsim_CA_bounded, 1, mean), .9999) | 
                                                   meanA_AR > quantile(apply(MVNsim_AR_bounded, 1, mean), .9999)) |
                                                  (meanA_CA < quantile(apply(MVNsim_CA_bounded, 1, mean), .0001) | 
                                                     meanA_AR < quantile(apply(MVNsim_AR_bounded, 1, mean), .0001)), 
                                                alpha("blue", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, blue = top 0.01% sim")
points(apply(MVNsim_CA_bounded, 1, mean), apply(MVNsim_AR_bounded, 1, mean), pch = 20, col = alpha("orange", .02))
dev.off()

# how likely are these shared outliers in my simulation?
MVNsim_AR_bounded_quantile <- quantile(meanA_MVNsim_AR_bounded, .99)
MVNsim_CA_bounded_quantile <- quantile(meanA_MVNsim_CA_bounded, .99)
table(meanA_MVNsim_AR_bounded > MVNsim_AR_bounded_quantile & 
        meanA_MVNsim_CA_bounded > MVNsim_CA_bounded_quantile)/sim_n
.01^2 # it's about six times as likely under the MVN to get a double outlier at top .01% than it would be under true independence between N and S America
# but still same order of magnitude. Also I maybe need to simulate a larger # of trials to get accuracy in the tails (only expect 1000 outliers each out of 100k, very small # overlap)
# when I simulated 1 million MVN's I got .000521 as the probability of outliers in both CA and AR at .01%
MVNsim_AR_bounded_quantile_low <- quantile(meanA_MVNsim_AR_bounded, .01)
MVNsim_CA_bounded_quantile_low <- quantile(meanA_MVNsim_CA_bounded, .01)
table(meanA_MVNsim_AR_bounded < MVNsim_AR_bounded_quantile_low & 
        meanA_MVNsim_CA_bounded < MVNsim_CA_bounded_quantile_low)/sim_n
# on the low side it's similar. some variance from low # expected hits in my simulations


# these outliers are pretty surprising under a MVN framework
# how about under a poisson binomial framework? even more surprising
# first read in individual alpha estimates for mean A ancestry
CA_indAlpha <- do.call(rbind,
                       lapply(CA_pops_included$population, function(p) 
  read.table(paste0("results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/", p, ".alpha.anc"),
                                            stringsAsFactors = F)))
AR_indAlpha <- do.call(rbind,
                       lapply(AR_pops_included$population, function(p) 
                         read.table(paste0("results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/", p, ".alpha.anc"),
                                    stringsAsFactors = F)))

# compare distribution between binomial and MVN for 10^5 simulations
set.seed(101) # use same seed
PoiBinsim_CA <- apply(do.call(rbind,
                        lapply(CA_indAlpha$A, 
                       function(alpha)
                         rbinom(n = sim_n, size = 2, prob = alpha)))/2, 2, mean)
PoiBinsim_AR <- apply(do.call(rbind,
                                 lapply(AR_indAlpha$A, 
                                        function(alpha)
                                          rbinom(n = sim_n, size = 2, prob = alpha)))/2, 2, mean)
PoiBinsim_combined <- (PoiBinsim_CA*sum(CA_pops_included$n_bees) + PoiBinsim_AR*sum(AR_pops_included$n_bees))/sum(c(CA_pops_included$n_bees, AR_pops_included$n_bees))
                                                                                                                    
# MVN simulation bounded with no covariances
set.seed(101)
MVNsim_no_cov <- mvrnorm(n = sim_n, 
                       mu = AR_CA_K$alpha, 
                       Sigma = diag(diag(AR_CA_K$K), length(AR_CA_K$alpha), length(AR_CA_K$alpha)), 
                       tol = 1e-6, 
                       empirical = FALSE, 
                       EISPACK = FALSE)
MVNsim_no_cov_bounded <- MVNsim_no_cov
MVNsim_no_cov_bounded[MVNsim_no_cov < 0] <- 0
MVNsim_no_cov_bounded[MVNsim_no_cov > 1] <- 1
MVNsim_AR_no_cov_bounded <- MVNsim_no_cov_bounded[ , AR_pops_included$population]
MVNsim_CA_no_cov_bounded <- MVNsim_no_cov_bounded[ , CA_pops_included$population]
# mean across individuals
meanA_MVNsim_AR_no_cov_bounded <- apply(MVNsim_no_cov_bounded[ , AR_pops_included$population], 1, 
                                      function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_no_cov_bounded <- apply(MVNsim_no_cov_bounded[ , CA_pops_included$population], 1, 
                                      function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
# combined mean across individuals
meanA_MVNsim_no_cov_bounded <- apply(MVNsim_no_cov_bounded[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                                   function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))


# plot comparison
sim_compare <- 
  data.frame(poisson_binomial = PoiBinsim_combined,
           MVN_no_covariance = meanA_MVNsim_no_cov_bounded,
           #MVN_no_bounds = meanA_MVNsim,
           MVN_with_covariance = meanA_MVNsim_zero_bounded,
           observed_data = sample(meanA, sim_n, replace = F)) %>% # downsample data to match length of simulations
  tidyr::gather(., "distribution", "A")
# plot
p_sim_compare <- sim_compare %>%
  ggplot(., aes(x = A, color = distribution)) +
  geom_density(lwd = 1) +
  theme_bw() +
  xlab("Combined sample mean African ancestry") +
  scale_color_discrete(name = "Distribution", labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
p_sim_compare
ggsave("plots/distribution_data_vs_poibin_vs_MVN_sim.png",
       plot = p_sim_compare,
       device = "png",
       width = 8, height = 4, units = "in")

for (path in c("plots/", 
               "../../bee_manuscript/figures/")){
png(paste0(path, "QQ_plot_data_against_PoiBin.png"),
    height = 8, width = 8, res = 300, units = "in")
qqplot(PoiBinsim_combined, meanA,
       main = "QQ plot fit - A Ancestry Combined Sample",
       xlab = "Poisson Binomial Simulation",
       ylab = "Observed",
       col = "grey")
abline(0, 1, col = "black")
dev.off()
}
#make_qqplot_lines(meanA_MVNsim_zero_bounded, legend = T)

# plot QQ-plots:
#QQ-plots:
png("plots/QQ_plots_against_PoiBin.png",
    height = 10, width = 15, res = 300, units = "in")
par(mfrow=c(2,2))
qqplot(PoiBinsim_combined, meanA_MVNsim_bounded,
       main = "Poisson Binomial vs. MVN distribution null")
abline(0, 1, col = "blue")
qqplot(PoiBinsim_combined, meanA,
       main = "QQ-plot: Poisson Binomial dist. vs. data, combined mean")
abline(0, 1, col = "blue")
qqplot(PoiBinsim_CA, meanA_CA,
       main = "QQ-plot: Poisson Binomial dist. vs. data, CA mean")
abline(0, 1, col = "blue")
qqplot(PoiBinsim_AR, meanA_AR,
       main = "QQ-plot: Poisson Binomial dist. vs. data, AR mean")
abline(0, 1, col = "blue")
dev.off()


# What are my false discovery rates?
# get a 1% and 5% FDR for jointly shared outliers under the MVN (based on simulations):
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

# set FDR thresholds for A ancestry in Africanized bees
test_fdr_shared_high <- sapply(1:length(meanA), function(x) # takes a while
  fdr_shared_high(a1 = meanA_AR[x], a2 = meanA_CA[x], pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
test_fdr_shared_low <- sapply(1:length(meanA), function(x) # takes a while
  fdr_shared_low(a1 = meanA_AR[x], a2 = meanA_CA[x], pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))

# putting FDRs together with data
meanA_both <- data.frame(AR = meanA_AR, CA = meanA_CA, both = meanA,
                         FDR_shared_high = test_fdr_shared_high,
                         FDR_shared_low = test_fdr_shared_low)
table(meanA_both$FDR_shared_high < .05)/nrow(meanA_both) # 0.45% of sites appear positively selected in both zones at 5% FDR  
summary(meanA_both$both[meanA_both$FDR_shared_high < .05 | is.na(meanA_both$FDR_shared_high)])
hist(meanA_both$both[meanA_both$FDR_shared_high < .05 | is.na(meanA_both$FDR_shared_high)])
hist(meanA_both$CA[meanA_both$FDR_shared_high < .05 | is.na(meanA_both$FDR_shared_high)])
hist(meanA_both$AR[meanA_both$FDR_shared_high < .05 | is.na(meanA_both$FDR_shared_high)])

meanA_both %>%  
  filter(FDR_shared_high < .1 | is.na(FDR_shared_high)) %>%
  ggplot(aes(x = CA, y = AR, color = FDR_shared_high)) +
  geom_point() # doesn't look quite right -- how do you do a joint probability false-discovery-rate?
summary(meanA_both$both)
summary(meanA_both$both[is.na(meanA_both$FDR_shared_high)])
summary(meanA_both$CA)
# what is the FDR for top 1% overlap between both zones? About 8% FDR
test_fdr_shared_high_quantiles <- sapply(c(0.9, 0.95, 0.99, 0.991, 0.992, 0.993, 0.994, 0.995, 0.999, 0.9999), function(x) # takes a while
  fdr_shared_high(a1 = quantile(meanA_AR, x), a2 = quantile(meanA_CA, x), pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
test_fdr_shared_high_quantiles
# alternatively I could walk along individual FDRs, rather than quantiles. 
# Or walk along sd's above mean ancestry for each zone
sd_CA <- sd(meanA_CA)
mu_CA <- mean(meanA_CA)
sd_AR <- sd(meanA_AR)
mu_AR <- mean(meanA_AR)
sd_range <- seq(0, 5, by = .01) # I can make this more precise later if I want

# get standard deviation cutoffs for FDR shared high loci
FDR_values = c(.1, .05, .01)
test_fdr_shared_high_sds <- sapply(sd_range, function(x)
  fdr_shared_high(a1 = mu_AR+x*sd_AR, a2 = mu_CA+x*sd_CA, pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
FDRs_shared_high_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_shared_high_sds<p], na.rm = T))

# shared low # standard deviations below mean
test_fdr_shared_low_sds <- sapply(sd_range, function(x)
  fdr_shared_low(a1 = mu_AR-x*sd_AR, a2 = mu_CA-x*sd_CA, pop1 = meanA_AR, pop2 = meanA_CA, 
                  sims1 = meanA_MVNsim_AR_bounded, sims2 = meanA_MVNsim_CA_bounded))
FDRs_shared_low_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_shared_low_sds<p], na.rm = T))

# high CA and high AR separately (less powerful):
test_fdr_CA_high_sds <- sapply(sd_range, function(x) 
  fdr_1pop_high(a = mu_CA+x*sd_CA, pop = meanA_CA, 
                  sims = meanA_MVNsim_CA_bounded))
FDRs_CA_high_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_CA_high_sds<p], na.rm = T))
test_fdr_AR_high_sds <- sapply(sd_range, function(x)
  fdr_1pop_high(a = mu_AR+x*sd_AR, pop = meanA_AR, 
                sims = meanA_MVNsim_AR_bounded))
FDRs_AR_high_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_AR_high_sds<p], na.rm = T))

# low CA and high AR separately (less powerful):
test_fdr_CA_low_sds <- sapply(sd_range, function(x) 
  fdr_1pop_low(a = mu_CA-x*sd_CA, pop = meanA_CA, 
                sims = meanA_MVNsim_CA_bounded))
FDRs_CA_low_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_CA_low_sds<p], na.rm = T))
# we do not find any low outliers in CA based on a .1, .05 or .01 FDR (underpowered)
test_fdr_AR_low_sds <- sapply(sd_range, function(x)
  fdr_1pop_low(a = mu_AR-x*sd_AR, pop = meanA_AR, 
                sims = meanA_MVNsim_AR_bounded))
FDRs_AR_low_sds <- sapply(FDR_values, function(p) min(sd_range[test_fdr_AR_low_sds<p], na.rm = T))
# get standard deviations for false-discovery rates
FDRs_sds = data.frame(FDR_values = FDR_values,
                  shared_high = FDRs_shared_high_sds,
                  shared_low = FDRs_shared_low_sds,
                  CA_high = FDRs_CA_high_sds,
                  CA_low = FDRs_CA_low_sds,
                  AR_high = FDRs_AR_high_sds,
                  AR_low = FDRs_AR_low_sds,
                  stringsAsFactors = F)
write.table(FDRs_sds, "results/FDRs_MVN_01_high_low_A_SDs.txt", quote = F, col.names = T, row.names = F, sep = "\t")
# translate SD cutoffs to % A ancestry cutoffs:
FDRs = data.frame(FDR_values = FDR_values,
                  shared_high_CA = mu_CA + FDRs_sds$shared_high*sd_CA,
                  shared_high_AR = mu_AR + FDRs_sds$shared_high*sd_AR,
                  shared_low_CA = mu_CA - FDRs_sds$shared_low*sd_CA,
                  shared_low_AR = mu_AR - FDRs_sds$shared_low*sd_AR,
                  CA_high = mu_CA + FDRs_sds$CA_high*sd_CA,
                  CA_low = mu_CA - FDRs_sds$CA_low*sd_CA,
                  AR_high = mu_AR + FDRs_sds$AR_high*sd_AR,
                  AR_low = mu_AR - FDRs_sds$AR_low*sd_AR,
                  stringsAsFactors = F)
write.table(FDRs, "results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", quote = F, col.names = T, row.names = F, sep = "\t")
FDRs <- read.table("results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", 
                   stringsAsFactors = F, header = T, sep = "\t")

# what percent of the genome matches these cutoffs?
table(meanA_CA > FDRs$CA_high[FDRs$FDR_values==.05])/length(meanA_CA) # 0.28% of loci (very few!)
table(sites$scaffold[meanA_CA > FDRs$CA_high[FDRs$FDR_values==.05]]) # where are they?

# about .3% of the genome is found in high shared sites at 5% FDR
table(meanA_CA > FDRs$shared_high_CA[FDRs$FDR_values==.05] & meanA_AR > FDRs$shared_high_AR[FDRs$FDR_values == .05])/length(meanA_CA) # 0.28% of loci (very few!)
table(sites$scaffold[meanA_CA > FDRs$shared_high_CA[FDRs$FDR_values==.1] & meanA_AR > FDRs$shared_high_AR[FDRs$FDR_values == .1]])


# make a grid where you calculate quantile for joint distribution based off of MVN simulation
grid_n <- 10 # smoother if you use more like 100
range_CA <- seq(from = range(meanA_CA)[1], to = range(meanA_CA)[2], length.out = grid_n)
range_AR <- seq(from = range(meanA_AR)[1], to = range(meanA_AR)[2], length.out = grid_n)
quantile_MVN <- matrix(0, grid_n, grid_n)
for (i in 1:grid_n){
  for (j in 1:grid_n){
    quantile_MVN[i, j] <- mean(abs(meanA_MVNsim_CA_zero_bounded - mu_CA) <= abs(range_CA[i] - mu_CA) & 
                                 abs(meanA_MVNsim_AR_zero_bounded - mu_AR) <= abs(range_AR[j] - mu_AR)) 
  }
}
png("plots/contour_plot_joint_outliers_CA_AR_zero_bounded_MVN.png",
    height = 6, width = 8, units = "in", res = 300)
contour(x = range_CA, y = range_AR, z = quantile_MVN, 
        levels = c(0.25, 0.5, .95, .99, 1),
        xlab = "mean A ancestry in CA hybrid zone",
        ylab = "mean A ancestry in AR hybrid zone",       
        main = "contour plot for joint outliers, MVN zero bounded simulation")
points(meanA_CA, meanA_AR, col = alpha("blue", .05), pch = 20)
dev.off()

# actually maybe what I want is simpler -- just a density plot for the MVN sim underneath my data points
ggplot() +
  geom_point(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
             aes(x = meanA_CA, y = meanA_AR)) +
  geom_density_2d(data = data.frame(MVN_CA = meanA_MVNsim_CA_zero_bounded,
                                    MVN_AR = meanA_MVNsim_AR_zero_bounded), 
                  aes(x=MVN_CA, y=MVN_AR))
# I can use stat_contour if I have a z gridded for what I want to display
png("plots/density_MVN_overlay_on_scatterplot.png", 
    height = 8, width = 8, units = "in", res = 300)
plot(meanA_CA, meanA_AR, col = alpha("grey", alpha = .3), pch = 20,
     main = "Ancestry frequencies in both hybrid zones across SNPs",
     xlab = "Mean African ancestry in California bees",
     ylab = "Mean African ancestry in Argentina bees")
#plot(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded, col = "grey", pch = 20)
emdbook::HPDregionplot(cbind(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded), 
                       prob = c(0.99, 0.95, 0.75, 0.5), 
                       n = 100, # number of grid points
                       #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
                       col=c("lightgreen", "salmon", "lightblue", "yellow"), 
                       lwd = 3, 
                       add = TRUE)
legend("bottomright", legend = c(0.99, 0.95, 0.75, 0.5), 
       title = "Neutral simulation density",
       lty = 1,
       lwd = 4,
       col=c("lightgreen", "salmon", "lightblue", "yellow"))
dev.off()

png("plots/density_MVN_overlay_on_scatterplot_w_ind_zone_quantiles.png", 
    height = 8, width = 8, units = "in", res = 300)
plot(meanA_CA, meanA_AR, col = alpha("grey", alpha = .3), pch = 20,
     main = "Ancestry frequencies in both hybrid zones across SNPs",
     xlab = "Mean African ancestry in California bees",
     ylab = "Mean African ancestry in Argentina bees")
#plot(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded, col = "grey", pch = 20)
emdbook::HPDregionplot(cbind(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded), 
                       prob = c(0.99, 0.95, 0.75, 0.5), 
                       n = 100, # number of grid points
                       #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
                       col=c("lightgreen", "salmon", "lightblue", "yellow"), 
                       lwd = 3, 
                       add = TRUE)

abline(v = quantile(meanA_MVNsim_CA_zero_bounded, c(0.01, 0.99)), col = "lightgreen")
abline(h = quantile(meanA_MVNsim_AR_zero_bounded, c(0.01, 0.99)), col = "lightgreen")
legend("bottomright", legend = c(0.99, 0.95, 0.75, 0.5), 
       title = "Neutral simulation (MVN) density",
       lty = 1,
       lwd = 4,
       col=c("lightgreen", "salmon", "lightblue", "yellow"))
dev.off()


# any genes in outlier regions?

# write a bed file with outlier regions (to compare with honeybee genes):
# ancestry calls are extended to halfway between any two calls and end at the position of the last/first call for a scaffold
# note: this does not extend all the way to the ends of the chromosomes (so some genes may not have ancestry calls)
A_plus_sites <- bind_cols(sites, data.frame(CA = meanA_CA, AR = meanA_AR, stringsAsFactors = F)) 
A_plus_sites$half_right <- c(diff(A_plus_sites$pos, 1), 1)/2 # halfway between that and next position

A_plus_sites$half_right[c(diff(A_plus_sites$scaffold_n, 1), 0) != 0] <- 0.5 # new scaffold
A_plus_sites$half_left <- c(0.5, A_plus_sites$half_right[1:(nrow(A_plus_sites) - 1)]) # halfway between focal locus and previous position
# make bed start and end points
A_plus_sites$end <- floor(A_plus_sites$pos + A_plus_sites$half_right)
A_plus_sites$start <- floor(A_plus_sites$pos - A_plus_sites$half_left)

# write bed file with mean ancestry for Argentina and California included bees
# also include whether a site meets a FDR threshold for selection
A_plus_sites %>%
  mutate(FDR_shared_high = ifelse(CA >= FDRs[FDRs$FDR_values == .01, "shared_high_CA"] & AR >= FDRs[FDRs$FDR_values == .01, "shared_high_AR"],
                                  .01, ifelse(CA >= FDRs[FDRs$FDR_values == .05, "shared_high_CA"] & AR >= FDRs[FDRs$FDR_values == .05, "shared_high_AR"],
                                              .05, ifelse(CA >= FDRs[FDRs$FDR_values == .1, "shared_high_CA"] & AR >= FDRs[FDRs$FDR_values == .1, "shared_high_AR"],
                                                          .1, NA)))) %>%
  mutate(FDR_CA_high = ifelse(CA >= FDRs[FDRs$FDR_values == .01, "CA_high"],
                              .01, ifelse(CA >= FDRs[FDRs$FDR_values == .05, "CA_high"],
                                          .05, ifelse(CA >= FDRs[FDRs$FDR_values == .1, "CA_high"],
                                                      .1, NA)))) %>%
  mutate(FDR_AR_high = ifelse(AR >= FDRs[FDRs$FDR_values == .01, "AR_high"],
                              .01, ifelse(AR >= FDRs[FDRs$FDR_values == .05, "AR_high"],
                                          .05, ifelse(AR >= FDRs[FDRs$FDR_values == .1, "AR_high"],
                                                      .1, NA)))) %>%
  mutate(FDR_shared_low = ifelse(CA <= FDRs[FDRs$FDR_values == .01, "shared_low_CA"] & AR <= FDRs[FDRs$FDR_values == .01, "shared_low_AR"],
                                 .01, ifelse(CA <= FDRs[FDRs$FDR_values == .05, "shared_low_CA"] & AR <= FDRs[FDRs$FDR_values == .05, "shared_low_AR"],
                                             .05, ifelse(CA <= FDRs[FDRs$FDR_values == .1, "shared_low_CA"] & AR <= FDRs[FDRs$FDR_values == .1, "shared_low_AR"],
                                                         .1, NA)))) %>%
  mutate(FDR_CA_low = ifelse(CA <= FDRs[FDRs$FDR_values == .01, "CA_low"],
                             .01, ifelse(CA <= FDRs[FDRs$FDR_values == .05, "CA_low"],
                                         .05, ifelse(CA <= FDRs[FDRs$FDR_values == .1, "CA_low"],
                                                     .1, NA)))) %>%
  mutate(FDR_AR_low = ifelse(AR <= FDRs[FDRs$FDR_values == .01, "AR_low"],
                             .01, ifelse(AR <= FDRs[FDRs$FDR_values == .05, "AR_low"],
                                         .05, ifelse(AR <= FDRs[FDRs$FDR_values == .1, "AR_low"],
                                                     .1, NA)))) %>%
  dplyr::select(scaffold, start, end, snp_id, AR, CA, FDR_shared_high, FDR_AR_high, FDR_CA_high, 
                FDR_shared_low, FDR_AR_low, FDR_CA_low) %>%
  write.table(., 
              "results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/mean_ancestry_AR_CA_included.bed", 
              sep = "\t", quote = F, col.names = F, row.names = F)


# I think I'm making some independence assumption here with the FDR rates
# It might be better to bin ancestry into broader windows, before calculating covariances.

# Also Harpur's study blasted the markers against the genome + 5kb to get coordinates for prior varroa-hygeine QTLs 

# plotting K covariance in ancestry matrix
melt(AR_CA_K$K) %>%
  #filter(!(Var1 == "Avalon_2014" | Var2 == "Avalon_2014")) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees AR & CA")
ggsave("plots/k_matrix_CA_AR_pops.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
# plot as correlation:
melt(cov2cor(AR_CA_K$K)) %>%
  #filter(!(Var1 == "Avalon_2014" | Var2 == "Avalon_2014")) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix correlation - bees AR & CA")
ggsave("plots/k_matrix_correlation_CA_AR_pops.png", 
       height = 6, width = 8, 
       units = "in", device = "png")

melt(AR_CA_K$K) %>% # filter out Avalon to get better color distinction between the others
  filter(!(Var1 == "Avalon_2014" | Var2 == "Avalon_2014")) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees AR & CA no Avalon")
ggsave("plots/k_matrix_CA_AR_pops_no_Avalon.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
# get rid of diagonal to better see negative covariances:
melt(AR_CA_K$K) %>%
  mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees AR & CA no within-pop var")
ggsave("plots/k_matrix_CA_AR_pops_no_within_pop.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
# no diagonal look at correlations:
melt(cov2cor(AR_CA_K$K)) %>%
  mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix correlations - bees AR & CA no within-pop var")
ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_correlations.png", 
       height = 6, width = 8, 
       units = "in", device = "png")



# take away diagonal contribution of binomial sampling variance
# the diagonal elements have variance due to drift + binomial sampling variance
# = drift var + (n/n-1)*n*p*q/(n^2) where n/n-1 is the effect of small sample size n and npq is the standard binomial variance
# p is the African ancestry allele frequency, q = 1-p and n = number of haplotypes = 2*n_bees in a population
# but then because I'm using variance in frequencies, not counts, I divide by n^2
# except this isn't quite right because it should be observed allele freq at a locus, not genomewide mean
n_bees_per_pop <- sapply(names(AR_CA_K$alpha), function(p) meta.pop[meta.pop$population == p, "n_bees"])
sampling_var_diag <- ((2*n_bees_per_pop/(2*n_bees_per_pop - 1))*n_bees_per_pop*2*AR_CA_K$alpha*(1 - AR_CA_K$alpha))/(2*n_bees_per_pop)^2 
AR_CA_K_minus_sampling <- AR_CA_K$K
for (i in 1:length(sampling_var_diag)){
  AR_CA_K_minus_sampling[i,i] <- AR_CA_K_minus_sampling[i,i] - sampling_var_diag[i]
}
AR_CA_K_minus_sampling %>% # not quite
  melt() %>%
  dplyr::filter(!(Var1 == "Avalon_2014" | Var2 == "Avalon_2014")) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees AR & CA")
# simulate the variance:
a = .5
n_binom = 10
n_sim = 10^5
# simulate 100 loci, just binomial variance
sim_pop1 <- rbinom(size = n_binom, n = n_sim, p = a)/n_binom  
var(sim_pop1)
mean(sim_pop1)
cov(sim_pop1, sim_pop1)
mean((sim_pop1 - mean(sim_pop1))^2) # observed variance
mean(sim_pop1)*(1-mean(sim_pop1))/n_binom # theoretical no small sample size correction - looks correct
mean(sim_pop1)*(1-mean(sim_pop1))/(n_binom-1) # theoretical w/ small sample size correction - incorrect
# now add normal noise around allele freq due to drift before binomial sampling:
sim_drift2 <- rnorm(n = n_sim, mean = a, sd = .1)
sim_pop2 <- rbinom(size = n_binom, n = sim_drift2, p = a)/n_binom
var(sim_drift2)
mean(sim_pop2)
var(sim_pop2) # observed variance
mean((sim_pop2 - mean(sim_pop2))^2)
mean(sim_pop2)*(1-mean(sim_pop2))/n_binom # mean theoretical


log(1/(1 - diag(AR_CA_K$K)/(AR_CA_K$alpha*(1 - AR_CA_K$alpha))))
log(diag(AR_CA_K$K)/(AR_CA_K$alpha*(1 - AR_CA_K$alpha)))

# goal: set cutoffs for when distribution doesn't match MVN anymore
# remove outliers SNPs beyond these cutoffs
# recalculate K matrix
# plot new K matrix without outliers
# group individuals into pops & calculate cov/(a0*(1-a0)*(aCA)*(aAR)) where
# a0 is the initial pulse into the shared population before the split (vary it 0.99-0.82)
# and (1 - aCA) is the additional European ancestry added to CA after the split
# (so observed African ancestry in CA today is alphaCA = a0*aCA)
# Likewise aAR is the additional European ancestry added to AR after the split
# with this covariance, calculate N for a 1 generation bottleneck. cov = t/(2N) 
# .. then I'll have to translate 2N into a bee population size because of haplo-diploidy

# It would be good to test robustness of any results with bootstraps for covariances
# as in, bootstrap over regions of the genome.

# simulate W-F evolution (drift):
# 10 generations, population size N = 100, starting freq = .5
start_freq = .5
n_gen = 10
N = 100
sim_WF <- function(n_gen, N, start_freq){
  freqs <- numeric(n_gen + 1)
  freqs[1] <- start_freq
  for (i in 2:(n_gen + 1)){ # each generation of drift
    freqs[i] <- rbinom(1, 2*N, p = freqs[i - 1])/(2*N)
  }
  return(freqs[n_gen+1])
}
wf1 <- sapply(1:10000, function(x) sim_WF(n_gen = 10, N = 100, start_freq = .5))
mean((wf1 - .5)^2)/(.5^2)
10/(2*100)
wf2 <- sapply(1:10000, function(x) sim_WF(n_gen = 20, N = 100, start_freq = .5))
mean((wf2 - .5)^2)/(.5^2)
20/(2*100)
wf3 <- sapply(1:10000, function(x) sim_WF(n_gen = 20, N = 1000, start_freq = .5))
mean((wf3 - .5)^2)/(.5^2)
20/(2*1000)
wf4 <- sapply(1:100000, function(x) sim_WF(n_gen = 5, N = 10000, start_freq = .2))
mean((wf4 - .2)^2)/(.2*.8)
5/(2*10000)
# add in sample variance for sample size 16 haplotypes
wf2_16hap <- sapply(wf2, function(x) sim_WF(n_gen = 1, N = 8, start_freq = x))
# and 10 haplotypes
wf2_10hap <- sapply(wf2, function(x) sim_WF(n_gen = 1, N = 5, start_freq = x))
cov(wf2_16hap, wf2_10hap)
var(wf2_10hap)
var(wf2)
AR_CA_K$K
cov(meanA_CA, meanA_AR)/(.85*.15)*200
cov(meanA_CA, meanA_AR)/(.5*.5)*200
cov_CA_AR <- mean((meanA_CA-mean(meanA_CA))*(meanA_AR - mean(meanA_AR)))
mean(meanA)
# t/2n = cov
# n = t/(cov * 2)
# so if I assume 1 generation of bottleneck, n is about 200:
# note: later I can correct for haploid/diploid
cov_CA_AR
1/(cov_CA_AR/(.85*.15)*2)
1/(cov_CA_AR/(mean(meanA)*(1-mean(meanA)))*2)
AR_CA_K$K
1/(.001/(mean(meanA)*(1-mean(meanA)))*2)
# how do I do the w/in pop small pop correction??
(mean((wf2_16hap - .5)^2))/(.5*.5)
(mean((wf2_16hap - .5)^2)*(15/16))/(.5*.5)
(mean((wf2_16hap - .5)^2)*(1/16))/(.5*.5)
(mean((wf2_16hap - .5)^2))/(.5*.5) - 1/16
(mean((wf2_16hap - .5)^2))/(.5*.5) - 1/16
(mean((wf2_16hap - .5)^2) - 1/16)/(.5*.5)
20/(2*100) + (.5^2)/16
hap16 <- sapply(rep(.5, 10000), function(x) sim_WF(n_gen = 1, N = 8, start_freq = x))
mean((hap16-.5)^2)
(mean((wf2_16hap - .2)^2))/(.5*.5) - mean((hap16-.5)^2)/(.5*.5)
((.5^2)/16)*1/16
var(hap16)
(.5^2)/16

# simulate quick pop split to confirm.


r <- read.table("results/SNPs/thin1kb_common3/included_pos_on_Wallberg_Amel4.5_chr.r.bed",
                sep = "\t", stringsAsFactors = F) %>%
  mutate(snp_id = read.table("results/SNPs/thin1kb_common3/included_pos_on_Wallberg_Amel4.5_chr.map",
                             stringsAsFactors = F, header = F, sep = "\t")$V2) %>%
  rename(chr = V1) %>%
  rename(cM_Mb = V4) %>%
  mutate(pos_chr = V2 + 1) %>% # because .bed files are 0 indexed
  dplyr::select(c("snp_id", "chr", "pos_chr", "cM_Mb"))

# **** check on the recombination rate files!

d <- bind_cols(sites, A) %>%
  left_join(., r, by = c("chr", "snp_id"))
d$r_bin10 = with(d, cut(cM_Mb,  # note need to load map from scaffolds_to_chr.R
                             breaks = unique(quantile(c(0, r$cM_Mb, 100), # I extend bounds so that every value gets in a bin
                                                      p = seq(0, 1, by = .1))),
                             right = T,
                             include.lowest = T))
table(d$r_bin10) # fewer SNPs represented in lowest recomb. bins
head(d[is.na(d$r_bin10),]) # good, should be none

d$r_bin5 = with(d, cut(cM_Mb,  # note need to load map from scaffolds_to_chr.R
                        breaks = unique(quantile(c(0, r$cM_Mb, 100), # I extend bounds so that every value gets in a bin
                                                 p = seq(0, 1, by = .2))),
                        right = T,
                        include.lowest = T))
table(d$r_bin5) # fewer SNPs represented in lowest recomb. bins
table(d$r_bin5, is_outlier)
head(d[is.na(d$r_bin5),]) # good, should be none

# make K matrix again with only high recombination rate region of genome

#a <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "(38.8,100]", ]))
#a2 <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 %in% c("[0,8.49]","(8.39,15.6]"),]))
a <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "(39.8,100]", ]))
a2 <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "[0,14]", ]))

#a <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "(39.8,100]" & !is_outlier, ]))
#a2 <- make_K_calcs(t(cbind(AR_A, CA_A)[d$r_bin10 == "[0,14]" & !is_outlier, ]))

# note: I need some way of summarising these K matrices statistically 
# to test hypotheses about low vs. high r regions.

melt(a$K) %>%
  mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  ggtitle("K matrix covariances r > 38cM/Mb - bees AR & CA no within-pop var")
ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_highr.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
melt(cov2cor(a$K)) %>%
  mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  ggtitle("K matrix correlations r > 38cM/Mb - bees AR & CA no within-pop var")
ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_highr_correlations.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
melt(a2$K) %>%
  mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  ggtitle("K matrix covariances r < 14cM/Mb - bees AR & CA no within-pop var")
ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_lowr.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
melt(cov2cor(a2$K)) %>%
  mutate(value = ifelse(Var1 == Var2, 0, value)) %>% # arbitarily set diagonal to zero
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +
  ggtitle("K matrix correlations r < 14cM/Mb - bees AR & CA no within-pop var")
ggsave("plots/k_matrix_CA_AR_pops_no_within_pop_lowr_correlations.png", 
       height = 6, width = 8, 
       units = "in", device = "png")


# a better null: maybe I should take the full covariance matrix for ind. ancestries, 
# subtract off the diagonal for what I know the binomial sampling variance to be,
# run a MVN, and then do binomial sampling at the end. 
# Graham *might* say this is excessive since we don't know if some variance comes from the ancestry caller
# but I think with the posteriors it's actually more conservative than a binomial
# which is more likely to get to the bounds 0/1 than a posterior. 
# (and I could use the ancestry call from ancestry_hmm if I wanted)
# also if I have enough bees, this might just give me the same answer as the MVN

AbyPopr10 <- d %>%
  gather("pop", "A", pops) %>%
  group_by(., pop, r_bin10) %>%
  summarise(meanA = mean(A))
AbyPopr10 %>%
  ggplot(aes(x = pop, y = meanA)) +
  geom_point() +
  facet_wrap(~r_bin10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

AbyPopr5 <- d %>%
  gather("pop", "A", pops) %>%
  group_by(., pop, r_bin5) %>%
  summarise(meanA = mean(A))
AbyPopr5 %>%
  ggplot(aes(x = pop, y = meanA, color = r_bin5)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# it would be much better to use NGSAdmix for different SNP sets
# to get global ancestry estimates

# plot ancestry across all populations by mean latitude in that pop, separate by hybrid zone
# I used Riverside 2014 coordinates for Riverside 1999. 
# Note there are two Riverside collection sites, close in latitude, but different in temp etc. due to elevation
a <- cbind(A_AR_CA, A) %>%
  left_join(., dplyr::select(sites, c("chr", "pos", "chr_n", "snp_id")), by = "snp_id") %>%
  tidyr::gather(., "population", "A_ancestry", colnames(A)) %>%
  left_join(., meta.pop, by = "population")
a_mean <- a %>%
  group_by(population) %>%
  summarise(A_ancestry = mean(A_ancestry)) %>%
  left_join(meta.pop, by = "population") %>%
  mutate(abs_lat = abs(lat)) %>%
  mutate(abs_lat_c = abs_lat - mean(abs_lat))

# plot some ancestry clines across latitude for top SNPs:
# start with high shared outliers (there's only 15 regions):
high.shared.outliers
i = 1
a_shared_high_outlier_.05sig <- a %>%
  mutate(combined = (AR + CA)/2) %>%
  dplyr::select(combined, scaffold, pos, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_shared_high != "NA") %>%
  unique() %>%
  filter(scaffold == high.shared.outliers[i, "chr"] &
           pos >= high.shared.outliers[i, "start"] &
           pos <= high.shared.outliers[i, "end"]) %>%
  dplyr::arrange(combined) %>%
  .[1, "snp_id"]

# ID SNPs from very top outliers:
shared_high_top_snp <- a %>%
  mutate(combined = (AR + CA)/2) %>%
  dplyr::select(combined, scaffold, pos, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_shared_high != "NA") %>%
  unique() %>%
  dplyr::arrange(desc(combined)) %>%
  .[1, "snp_id"]
shared_low_top_snp <- a %>%
  mutate(combined = (AR + CA)/2) %>%
  dplyr::select(combined, scaffold, pos, FDR_shared_low, snp_id) %>%
  dplyr::filter(FDR_shared_low != "NA") %>%
  unique() %>%
  dplyr::arrange(combined) %>%
  .[1, "snp_id"]

# just AR
AR_high_only_top_snp <- a %>%
  dplyr::select(scaffold, pos, AR, CA, FDR_AR_high, FDR_CA_high, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_AR_high != "NA" & FDR_shared_high == "NA") %>%
  unique() %>%
  dplyr::arrange(desc(AR)) %>%
  .[1, "snp_id"]
AR_low_only_top_snp <- a %>%
  dplyr::select(scaffold, pos, AR, CA, FDR_AR_low, FDR_CA_low, FDR_shared_low, snp_id) %>%
  dplyr::filter(FDR_AR_low != "NA" & FDR_shared_low == "NA") %>%
  unique() %>%
  dplyr::arrange(AR) %>%
  .[1, "snp_id"]
# just CA
CA_high_only_top_snp <- a %>%
  dplyr::select(scaffold, pos, AR, CA, FDR_AR_high, FDR_CA_high, FDR_shared_high, snp_id) %>%
  dplyr::filter(FDR_CA_high != "NA" & FDR_shared_high == "NA") %>%
  unique() %>%
  dplyr::arrange(desc(CA)) %>%
  .[1, "snp_id"]
# put together into data frame:
some_snp_outliers <- data.frame(
  snp_id = c(a_shared_high_outlier_.05sig,
             shared_high_top_snp,
             shared_low_top_snp,
             AR_high_only_top_snp,
             AR_low_only_top_snp,
             CA_high_only_top_snp),
  label = c("a 5% FDR outlier shared high A ancestry",
            "the top outlier shared high A ancestry",
            "the top outlier shared low A ancestry",
            "the top outlier Argentina-only high A ancestry",
            "the top outlier Argentina-only low A ancestry",
            "the top outlier California-only high A ancestry"),
  name = c("a_shared_high_outlier_.05sig",
           "shared_high_top_snp",
           "shared_low_top_snp",
           "AR_high_only_top_snp",
           "AR_low_only_top_snp",
           "CA_high_only_top_snp"))
  
# make plots:
for (i in 1:nrow(some_snp_outliers)){
my_plot <- a %>%
    filter(snp_id == some_snp_outliers[i, "snp_id"]) %>%
    ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
    geom_point() +
    geom_point(data = a_mean, shape = 2) +
    ggtitle(paste0("cline across latitude for ", some_snp_outliers[i, "label"], ": ", some_snp_outliers[i, "snp_id"])) +
    xlab("absolute latitude") +
    ylab("population mean African ancestry") +
    facet_grid(zone ~ .)
ggsave(paste0("plots/outlier_clines_", some_snp_outliers[i, "name"], ".png"),
       plot = my_plot,
       height = 5, 
       width = 10, 
       units = "in", 
       device = "png")
}

# new approach: get 1 SNP per outlier. So find top SNP in each outlier region:





# Is this consistent with no spatial pattern to selection? Either selected only in Brazil.
# or selected at low selection coefficient across the whole S. American hybrid zone
# ok so in Brazil bees are 84% A ancestry and at this locus they are 54% A ancestry.
# so if all the selection happened in Brazil, then spread South neutrally,
# we'd expect 64% (54/84) x (mean A ancestry genomewide) at this locus for S. American bees
# This is a conservative test of early selection in Brazil (assuming it all happened early), 
# then neutral spread, and it doesn't fit well in Argentina (too high A) or CA (too low A). 
a_mean_64percent <- a_mean %>%
  #filter(zone == "S. America") %>%
  mutate(A_ancestry = .64*A_ancestry)
a %>%
  filter(snp_id == AR_low_only_top_snp) %>%
  #filter(zone == "S. America") %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
  geom_point() +
  geom_point(data = a_mean_64percent, shape = 2) +
  ggtitle(paste0("cline across latitude for ", "low region in AR compared to 64% mean A ancestry", ": ", AR_low_only_top_snp)) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  facet_grid(zone ~ .)
ggsave(paste0("plots/outlier_clines_low_AR_against_64perc_dilution_hypothesis.png"),
       height = 5, 
       width = 10, 
       units = "in", 
       device = "png")

# fit mean ancestry beta regression model
m_lat.hmm <- map( # quadratic approximation of the posterior MAP
  alist(
    A_ancestry ~ dbeta2(prob = p, theta = theta),
    logit(p) <- mu + b_lat*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 30),
    b_lat ~ dnorm(0, 5)
  ),
  data = data.frame(a_mean, stringsAsFactors = F),
  start = list(mu = -1, theta = 5, b_lat = 0.5))
precis(m_lat.hmm)
pairs(m_lat.hmm)
m_lat.hmm.betareg <- betareg(A_ancestry ~ abs_lat_c,
                         link = "logit",
                         data = a_mean)

# m_lat best single-predictor model:
full_range_data_hmm = list(abs_lat_c = seq(from = range(a_mean$abs_lat_c)[1], 
                                       to = range(a_mean$abs_lat_c)[2], 
                                       length.out = 20),
                       abs_lat = seq(from = range(a_mean$abs_lat_c)[1], 
                                     to = range(a_mean$abs_lat_c)[2], 
                                     length.out = 20) + mean(a_mean$abs_lat))
m_lat.hmm.post <- sim(m_lat.hmm, full_range_data_hmm, n = 10000)
m_lat.hmm.post.mu <- apply(m_lat.hmm.post, 2, mean)
m_lat.hmm.post.HDPI <- apply(m_lat.hmm.post, 2, HPDI, prob=0.95)



# plot model prediction for m_lat:
png("plots/m_lat_model_prediction_hmm_with_top_shared_outlier.png", height = 6, width = 8, units = "in", res = 300)
plot(A_ancestry ~ abs_lat, data = a_mean, col = ifelse(a_mean$zone == "S. America", rangi2, "skyblue"),
     ylim = c(0, 1), main = "African ancestry predicted by latitude", xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion (ancestry_hmm)")

# plot MAP line from model coefficients:
curve(logistic(coef(m_lat.hmm)["mu"] + (x - mean(a_mean$abs_lat))*coef(m_lat.hmm)["b_lat"]), 
               from = range(a_mean$abs_lat)[1], to = range(a_mean$abs_lat)[2], n = 1000, add = T,
      col = "darkgrey")
# plot a shaded region for 95% HPDI
shade(m_lat.hmm.post.HDPI, full_range_data_hmm$abs_lat)
# plot snp of interest:
a %>%
  filter(snp_id == shared_high_top_snp) %>%
  with(., points(x = abs(lat), y = A_ancestry, pch = 17,
                 col = ifelse(a_mean$zone == "S. America", rangi2, "skyblue")))
legend("topright", c("S. America", "N. America", "model prediction (MAP)", "95% confidence (HPDI)"),
       pch = c(1, 1, NA, 15), lty = c(NA, NA, 1, NA), col = c(rangi2, "skyblue", "black", "grey"))
dev.off()

# get logistic curve for outlier SNP:
m_lat.shared.high.outlier.betareg <- a %>%
  filter(snp_id == shared_high_top_snp) %>%
  mutate(abs_lat_c = abs(lat) - mean(abs(lat))) %>%
  betareg(A_ancestry ~ abs_lat_c,
          link = "logit",
          data = .)
# plot curve
# curve:
curve(logistic(coef(m_lat.shared.high.outlier.betareg)["(Intercept)"] + 
                 (x - mean(a_mean$abs_lat))*coef(m_lat.shared.high.outlier.betareg)["abs_lat_c"]), 
      from = range(a_mean$abs_lat)[1], to = range(a_mean$abs_lat)[2], n = 1000, add = F,
      main = "estimated logistic cline at top shared high A outlier SNP",
      xlab = "degrees latitude from the equator",
      ylab = "mean population A ancestry")
a %>%
  filter(snp_id == shared_high_top_snp) %>%
  with(., points(x = abs(lat), y = A_ancestry, pch = 17,
                 col = ifelse(a_mean$zone == "S. America", rangi2, "skyblue")))

# it's hard to fit a logistic because I'm only observing a small portion of the cline for this SNP
# basically it's not bounded high or low, just consistently high A across the zone




# TO DO: I want to make a plot for top SNP within each peak region (only a few regions really, I could include them all):
# also add the logit fit for mean african ancestry (from ancestry_hmm calls) to the plot. 
# For the supplement, I need to compare the ancestry_hmm mean with the NGSadmix mean

# For each peak, ID the top SNP. Then plot it's frequency across latitude with logit in background.

# For low ancestry peaks, check if they are M vs. C

# what are my tract lengths like? Start by looking at high confidence tracts > 0.8
some_inds <- c("AR0101", "AR0506", "AR1112", "AR1910", "AR2701", "CA0529", "CA1109", "CA1419", "SRCD6C")
some_tracts <- lapply(some_inds, function(ind) 
  read.table(paste0("results/tracts/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/AA/", ind, ".bed"),
                          stringsAsFactors = F, header = F, sep = "\t") %>%
  data.table::setnames(c("chr", "start", "end")) %>%
    mutate(length = end - start))
par(mfrow=c(3,3))
lapply(1:9, function(x) hist(some_tracts[[x]]$length, freq = F, xlim = c(0, 500000), main = some_inds[x]))
par(mfrow=c(1,1))
# I'm not seeing a bunch (a higher %) of super short tracts for low A individuals
# (the high confidence filter probably takes care of this)
lapply(1:9, function(x) table(some_tracts[[x]]$length <= 10000))
lapply(1:9, function(x) table(some_tracts[[x]]$length <= 5000))
lapply(1:9, function(x) table(some_tracts[[x]]$length <= 1000))
lapply(1:9, function(x) table(some_tracts[[x]]$length <= 2000))





# compare the priors with the mean inferred ancestry from ancestry_hmm
# make individual ancestry plots too:
compare_anc_ind <- rbind(CA_indAlpha, AR_indAlpha) %>% # ancestry_hmm genomewide estimates for individuals
  tidyr::gather(., "ancestry", "ancestry_hmm", c("A", "C", "M")) %>%
  left_join(., tidyr::gather(d_ACM[ , c("Bee_ID", "A", "C", "M")], "ancestry", "NGSAdmix", c("A", "C", "M")), # NGSAdmix genomewide estimates for individuals
            by = c("ancestry"="ancestry", "ID"="Bee_ID"))
compare_anc_pop <- admix_proportions %>% # ancestry_hmm genomewide estimates for individuals
  tidyr::gather(., "ancestry", "ancestry_hmm", c("A", "C", "M")) %>%
  left_join(., tidyr::gather(d_pop_ACM[ , c("population", "A", "C", "M")], "ancestry", "NGSAdmix", c("A", "C", "M")), # NGSAdmix genomewide estimates for individuals
            by = c("ancestry"="ancestry", "population"="population"))
p_ind <- compare_anc_ind %>%
  ggplot(., aes(NGSAdmix, ancestry_hmm, color = ancestry)) +
  geom_point(alpha = .5) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
  xlab("NGSAdmix") +
  ylab("ancestry_hmm") +
  ggtitle("Individual ancestry estimates")
p_pop <- compare_anc_pop %>%
  ggplot(., aes(NGSAdmix, ancestry_hmm, color = ancestry)) +
  geom_point(pch = 1) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
  xlab("NGSAdmix") +
  ylab("ancestry_hmm") +
  ggtitle("Population ancestry estimates")
p_pop_ind <- grid.arrange(p_pop + theme(legend.position = "none"), p_ind, nrow = 1, ncol = 2, right = 2)
p_pop_ind
# save plot locally and in bee_manuscript figures folder.
ggsave("plots/mean_ancestry_prior_posterior_ancestry_hmm.png",
       plot = p_pop_ind,
       device = "png",
       width = 8, height = 5, units = "in")
ggsave("../../bee_manuscript/figures/mean_ancestry_prior_posterior_ancestry_hmm.png",
       plot = p_pop_ind,
       device = "png",
       width = 8, height = 5, units = "in")

# I can use the older version (below) if I just want a population-level plot

colors_3 <- scales::hue_pal()(3)

for (path in c("plots/mean_pop_ancestry_prior_posterior_ancestry_hmm.png", 
               "../../bee_manuscript/figures/mean_pop_ancestry_prior_posterior_ancestry_hmm.png")){
  png(path,
      height = 5, width = 6, units = "in", res= 300)
  left_join(admix_proportions, filter(admix_times, ancestry == "A"), by = "population") %>%
    with(., plot(proportion, A, xlim = 0:1, ylim = 0:1, col = colors_3[1], 
                 xlab = "mean ancestry prior (NGSAdmix)",
                 ylab = "mean ancestry posterior (ancestry_hmm)", 
                 main = "Effect of HMM on inferred mean population ancestry"))
  left_join(admix_proportions, filter(admix_times, ancestry == "M"), by = "population") %>%
    with(., points(proportion, M, xlim = 0:1, ylim = 0:1, col = colors_3[3]))
  left_join(admix_proportions, filter(admix_times, ancestry == "C"), by = "population") %>%
    with(., points(proportion, C, xlim = 0:1, ylim = 0:1, col = colors_3[2]))
  abline(0, 1)
  legend("bottomright", legend = c("A", "C", "M"),
         col = colors_3,
         pch = 1,
         title = "Ancestry")
  dev.off()
}

