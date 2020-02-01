# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(reshape2)
library(viridis)
library(reshape2) # melt function
library(LaplacesDemon)
library(emdbook)
library(betareg)
library(gridExtra)
library(MASS) # for mvrnorm
library(ggpointdensity)
library(coda) # for hpdi calc
#library(hex) # for hex plot
#library(ggridges) # to get density plot without bottom line on x axis
source("../colors.R") # for color palette
source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions
#source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions
source("calc_FDRs.R") # scripts to calculate false discovery rates

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
#pops <- read.table("../bee_samples_listed/byPop/pops_included.list", stringsAsFactors = F)$V1
pops <- read.table("../bee_samples_listed/byPop/combined_sept19_pops.list", stringsAsFactors = F)$V1
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
  left_join(data.frame(population = pops, stringsAsFactors = F), # limits to included populations
            ., by = "population") %>%
  # not relevant anymore, not using Riverside 1999 historical samples:
  mutate(lat = ifelse(population == "Riverside_1999", .[.$population == "Riverside_2014", "lat"], lat)) %>%
  mutate(zone = ifelse(group == "AR_2018", "S. America", "N. America")) %>%
  arrange(lat)

# included populations by latitude:
pops_by_lat <- meta.pop$population[order(meta.pop$lat)]
meta.AR.order.by.lat <- data.frame(population = pops_by_lat, stringsAsFactors = F) %>%
  left_join(., meta.pop, by = "population") %>%
  filter(zone == "S. America") %>%
  mutate(abs_lat_SA_c = abs(lat) - mean(abs(lat))) # absolute latitude centered for SA
save(file = "results/pops_by_lat.RData", list = c("pops_by_lat", "meta.pop", "meta.AR.order.by.lat"))
save(file = "results/meta.RData", list = c("meta.ind", "meta.pop", "pops_by_lat", "meta.AR.order.by.lat"))

# get SNP sites where ancestry was called
sites0 <- read.table("results/SNPs/combined_sept19/chr.var.sites", stringsAsFactors = F,
                     sep = "\t", header = F)[ , 1:2]
colnames(sites0) <- c("scaffold", "pos")
chr_lengths <- cbind(read.table("../data/honeybee_genome/chr.names", stringsAsFactors = F),
                     read.table("../data/honeybee_genome/chr.lengths", stringsAsFactors = F)) %>%
  data.table::setnames(c("chr", "scaffold", "chr_length")) %>%
  mutate(chr_n = 1:16) %>%
  mutate(chr_end = cumsum(chr_length)) %>%
  mutate(chr_start = chr_end - chr_length) %>%
  mutate(chr_mid = (chr_start + chr_end)/2)

sites <- left_join(sites0, chr_lengths[ , c("chr", "scaffold", "chr_n", "chr_start")], by = "scaffold") %>%
  mutate(cum_pos = pos + chr_start)


# get ancestry frequencies for each population across the genome
dir_results <- "results/ancestry_hmm/combined_sept19/posterior"

popA <- lapply(pops_by_lat, function(p) read.table(paste0(dir_results, "/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A <- do.call(cbind, popA)
rm(popA)
colnames(A) <- pops_by_lat

popM <- lapply(pops_by_lat, function(p) read.table(paste0(dir_results, "/anc/", p, ".M.anc"),
                                            stringsAsFactors = F))
M <- do.call(cbind, popM)
colnames(M) <- pops_by_lat
rm(popM)
popC <- lapply(pops_by_lat, function(p) read.table(paste0(dir_results, "/anc/", p, ".C.anc"),
                                            stringsAsFactors = F))
C <- do.call(cbind, popC)
colnames(C) <- pops_by_lat
rm(popC)

# mean ancestry across populations
meanA0 <- apply(A, 1, mean)
meanC0 <- apply(C, 1, mean)
meanM0 <- apply(M, 1, mean)

# mean ancestry across the genome for each pop
admix_proportions <- data.frame(population = pops_by_lat,
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
# admix.pops and admix.ind are loaded from plot_NGSadmix.R script
left_join(admix.pops, filter(admix_times, ancestry == "A"), by = "population") %>%
  with(., plot(A, proportion, xlim = 0:1, ylim = 0:1, col = "blue", main = "should be perfect match -- showing prior"))
abline(0, 1)


# plot ancestry times vs. mean ancestry proportion:
admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  filter(., ancestry != "C") %>% # no time, first ancestry
  ggplot(., aes(x = abs(lat), y = time, color = zone, shape = factor(year))) +
  geom_point(alpha = .75) +
  facet_grid(. ~ ancestry) +
  scale_color_manual(values = col_NA_SA_both, name = "Hybrid zone") +
  #ggtitle("Inferred time of admixture pulses from HMM") +
  ylab("Time (generations)") +
  xlab("Degrees latitude from the equator") +
  labs(shape = "Year") +
  theme_classic() +
  scale_shape_manual(values = c(17, 19))
ggsave("plots/time_of_admixture_vs_latitude.png",
       height = 3, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures/time_of_admixture_vs_latitude.png",
       height = 3, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/time_of_admixture_vs_latitude.tiff",
       height = 3, width = 5.2, units = "in", dpi = 600)

# compare California and Argentina: 
# For the same ancestry proportions, do they have similar inferred times of admixture?
# This would be evidence for different amounts of ongoing gene flow in the two zones.
admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  #filter(., ancestry == "A") %>%
  filter(., ancestry != "C") %>%
  mutate(ancestry_labels = paste(ancestry, "ancestry admixture pulse")) %>%
  ggplot(., aes(x = proportion, y = time, color = zone, shape = factor(year))) +
  geom_point(alpha = .75) +
  #ggtitle("Time since admixture estimated from ancestry block lengths") +
  facet_grid(. ~ ancestry) +
  scale_color_manual(values = col_NA_SA_both, name = "Hybrid zone") +
  ylab("Time (generations in the past)") +
  xlab("Admixture proportion") +
  labs(color = "Hybrid Zone", shape = "Year") +
  theme_classic() +
  scale_shape_manual(values = c(17, 19))
ggsave("plots/California_has_shorter_ancestry_blocks.png",
       height = 3, width = 6, units = "in")
ggsave("../../bee_manuscript/figures/California_has_shorter_ancestry_blocks.png",
       height = 3, width = 5.2, units = "in", dpi = 600, device = "pdf")
ggsave("../../bee_manuscript/figures_supp/California_has_shorter_ancestry_blocks.tiff",
       height = 3, width = 5.2, units = "in", dpi = 600)


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
A_bees_all8 <- lapply(bees_all8$Bee_ID, function(p) read.table(paste0("Amel4.5_results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A_all8 <- do.call(cbind, A_bees_all8)
colnames(A_all8) <- bees_all8$Bee_ID


#CA_pops_included <- meta.pop[meta.pop$group %in% c("N_CA", "S_CA", "CA_2018") & meta.pop$year >= 2014, ]
CA_pops_included <- meta.pop[meta.pop$zone == "N. America", ]
CA_A <- A[, CA_pops_included$population]
CA_A_earlier <- A[, meta.pop$population[meta.pop$group %in% c("N_CA", "S_CA") & meta.pop$year < 2014]]
CA_A_2014 <- A[, meta.pop$population[meta.pop$group %in% c("N_CA", "S_CA") & meta.pop$year == 2014]]
CA_A_2018 <- A[, meta.pop$population[meta.pop$group %in% c("CA_2018")]]
plot(apply(CA_A_2014, 1, mean), apply(CA_A_2018, 1, mean))
plot(apply(CA_A_earlier, 1, mean), apply(CA_A_2018, 1, mean))
AR_pops_included <- meta.pop[meta.pop$zone == "S. America", ]
AR_A <- A[, AR_pops_included$population]

# take mean across individuals, not across populations:
meanA_CA0 <- apply(CA_A, 1, mean)
meanA_CA <- apply(CA_A, 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanA_AR0 <- apply(AR_A, 1, mean)
meanA_AR <- apply(AR_A, 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))

# with only included pops, mean of all individuals:
meanA <- apply(A, 1, function(x) sum(x*meta.pop$n_bees)/sum(meta.pop$n_bees))

# C and M ancestries:
meanC_CA <- apply(C[ , CA_pops_included$population], 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanC_AR <- apply(C[ , AR_pops_included$population], 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanC <- apply(C, 1, function(x) sum(x*meta.pop$n_bees)/sum(meta.pop$n_bees))
meanM_CA <- apply(M[ , CA_pops_included$population], 1, function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))
meanM_AR <- apply(M[ , AR_pops_included$population], 1, function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanM <- apply(M, 1, function(x) sum(x*meta.pop$n_bees)/sum(meta.pop$n_bees))

# save ancestry data for later access:
save(file = "results/A.RData", list = c("A", "sites", "meanA", "meanA_CA", "meanA_AR"))
save(file = "results/C.RData", list = c("C", "sites", "meanC", "meanC_CA", "meanC_AR"))
save(file = "results/M.RData", list = c("M", "sites", "meanM", "meanM_CA", "meanM_AR"))

# make plots - all loci, mean ancestry across all pops
png("plots/histogram_A_ancestry_all_loci.png", height = 6, width = 8, units = "in", res = 300)
hist(meanA, main = "all loci: mean ancestry across populations")
abline(v = max(meanA), col = "blue")
abline(v = min(meanA), col = "blue")
abline(v = quantile(meanA, .99), col = "orange")
abline(v = quantile(meanA, .01), col = "orange")
dev.off()

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
  #filter(chr == "Group1" & scaffold_n %in% 20:23) %>%
  filter(chr == "Group1") %>%
  gather("pop", "Afreq", c("CA", "AR")) %>%
  ggplot(aes(x = pos, y = Afreq, color = pop)) +
  geom_point(size = .1) +
  #facet_wrap(~scaffold_n) + 
  ggtitle("African ancestry frequencies in CA and AR")
ggsave("plots/mean_A_ancestry_CA_and_AR_Group1_scaffolds_20-23.png", 
       device = "png", height = 10, width = 12, units = "in")


# plot K matrix (! order by latitude!)
zAnc_bees = make_K_calcs(t(A[ , pops_by_lat]))

save(file = "results/zAnc.RData", list = "zAnc_bees")
# what are mean correlations?
# within CA
# within AR S
# within AR N
# CA -> AR S
# CA -> AR N
# AR S -> AR N
AR_pops_S <- names(zAnc_bees$alpha[zAnc_bees$alpha < .5 & names(zAnc_bees$alpha) %in% AR_pops_included$population])
AR_pops_N <- names(zAnc_bees$alpha[zAnc_bees$alpha >= .5 & names(zAnc_bees$alpha) %in% AR_pops_included$population])
CA_pops <- CA_pops_included$population

get_mean_from_K <- function(K){ # for covariance use K, for correlation set K = cov2cor(K)
  lower_tri <- K
  lower_tri[lower.tri(lower_tri, diag = T)] <- NA # ignore variances and repeats
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
mean_corr_k <- get_mean_from_K(cov2cor(zAnc_bees$K))
mean_corr_k
mean_corr_k %>%
  write.table(., "results/mean_anc_corr_grouped.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
mean_cov_k <- get_mean_from_K(zAnc_bees$K)
mean_cov_k
mean_cov_k %>%
  write.table(., "results/mean_anc_cov_grouped.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
# plot the lower triangle of K
lower_tri <- zAnc_bees$K
lower_tri[lower.tri(lower_tri, diag = T)] <- NA
melt(lower_tri) %>% 
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()




melt(zAnc_bees$K) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_viridis() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("K covariance - bees")
ggsave("plots/k_matrix_all_pops.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
melt(zAnc_bees$K) %>%
  filter(Var1 != Var2) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_viridis() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("K covariance - bees")
ggsave("plots/k_matrix_all_pops_no_var.png", 
       height = 6, width = 8, 
       units = "in", device = "png")
k_plot_all <- melt(cov2cor(zAnc_bees$K)) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     #limits = c(-.35, .55)) +
                     limits = c(-.12, .32),
                     name = "") +
  xlab("") +
  ylab("") +
  ggtitle("K correlation matrix - bees")
k_plot_all
ggsave("plots/k_correlation_matrix_all_pops.png", 
       plot = k_plot_all,
       height = 6, width = 8, 
       units = "in", device = "png")
ggsave("../../bee_manuscript/figures/k_correlation_matrix_all_pops.pdf", 
       plot = k_plot_all + theme(axis.text.x = NULL, axis.text.y = NULL) + ggtitle(""),
       height = 7, width = 8, 
       units = "in", device = "pdf")

# make a new K matrix but omit outlier points:
ZAnc_bees_noOutliers = make_K_calcs(t(A[!(meanA > quantile(meanA, .9) | 
                                            meanA < quantile(meanA, .1)), 
                                        pops_by_lat]))
melt(cov2cor(ZAnc_bees_noOutliers$K)) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     #limits = c(-.35, .55)) +
                     limits = c(-.12, .32)) +
  xlab("") +
  ylab("") +
  ggtitle("K correlation matrix - bees")
ggsave("plots/k_NO_OUTLIERS_correlation_matrix_all_pops.png", 
       height = 7, width = 8, 
       units = "in", device = "png")
ggsave("../../bee_manuscript/figures/k_NO_OUTLIERS_correlation_matrix_all_pops.pdf", 
       height = 7, width = 8, 
       units = "in", device = "pdf")


summary(meanA)
sd(meanA)

# simulate from MVN distribution:
#n_sim <- 10^6 # slow
n_sim <- 10^5
#n_sim = length(meanA) # simulate data of same length
set.seed(101)
MVNsim <- mvrnorm(n = n_sim, 
                  mu = zAnc_bees$alpha, 
                  Sigma = zAnc_bees$K, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)



# sets bounds at 0 and 1
MVNsim_bounded <- data.frame(MVNsim, stringsAsFactors = F) 
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
# before bounding:
# mean across individuals
meanA_MVNsim_AR_unbounded <- apply(MVNsim_AR[ , AR_pops_included$population], 1, 
                                 function(x) sum(x*AR_pops_included$n_bees)/sum(AR_pops_included$n_bees))
meanA_MVNsim_CA_unbounded <- apply(MVNsim_CA[ , CA_pops_included$population], 1, 
                                 function(x) sum(x*CA_pops_included$n_bees)/sum(CA_pops_included$n_bees))


# combined mean across individuals
meanA_MVNsim_bounded <- apply(cbind(MVNsim_AR_bounded, MVNsim_CA_bounded)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))
meanA_MVNsim <- apply(cbind(MVNsim_AR, MVNsim_CA)[ , c(AR_pops_included$population, CA_pops_included$population)], 1, 
                               function(x) sum(x*c(AR_pops_included$n_bees, CA_pops_included$n_bees))/sum(c(AR_pops_included$n_bees, CA_pops_included$n_bees)))

# save MVN simulation as data objects
save(MVNsim_bounded, meanA_MVNsim_bounded, meanA_MVNsim_AR_bounded, meanA_MVNsim_CA_bounded,
     file = "results/MVNsim_bounded.RData")

# Add in constraint that all covariances must be 0 or positive
# (effectively I zero out negative covariances):
summary(zAnc_bees$K[zAnc_bees$K < 0])
sum(zAnc_bees$K < 0)/prod(dim(zAnc_bees$K)) # ~13% of covariances < 0
summary(zAnc_bees$K[zAnc_bees$K > 0])
# conclusion: negative covariances are on the same order as positive covariances
summary(zAnc_bees$K[T])
summary(as.matrix(zAnc_bees$K)[lower.tri(as.matrix(zAnc_bees$K), diag = F)]) # don't include diagonals
# simulate from MVN distribution:
K_zero <- zAnc_bees$K
K_zero[zAnc_bees$K < 0] <- 0
set.seed(101)
MVNsim_zero <- mvrnorm(n = n_sim, 
                  mu = zAnc_bees$alpha, 
                  Sigma = K_zero, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)

# how many simulations exceed the bounds?
# overall few, but more than 25% for some populations
sum(MVNsim < 0)/sum(table(MVNsim<0)) # low outliers need to be set to bound
sum(MVNsim > 1)/sum(table(MVNsim>1))
sum(MVNsim_zero < 0)/sum(table(MVNsim_zero<0)) # high outliers need to be set to bound
sum(MVNsim_zero > 1)/sum(table(MVNsim_zero>1))
hist(MVNsim[MVNsim_zero<0])
table(MVNsim_bounded[MVNsim_zero < 1])
apply(MVNsim, 2, function(x) sum(x > 1))/(nrow(MVNsim_zero))
summary(apply(MVNsim, 2, function(x) sum(x <0))/(nrow(MVNsim_zero)))
perc_sim_freq_out_of_bounds = data.frame(population = colnames(MVNsim),
           Lower = apply(MVNsim, 2, function(x) sum(x < 0))/nrow(MVNsim),
           Upper = apply(MVNsim, 2, function(x) sum(x > 1))/nrow(MVNsim),
           stringsAsFactors = F) %>%
  pivot_longer(data = ., cols = c("Lower", "Upper"), names_to = "bound", values_to = "p") %>% 
  left_join(., admix_proportions, by = "population") %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = A, y = p, color = zone)) +
  geom_point(alpha = .75) +
  facet_wrap(~bound) +
  xlab("Mean A ancestry proportion (ancestry_hmm)") +
  ylab("Percent sims exceeding bound") +
  #ggtitle("Percent simulated population ancestry frequencies < 0 (MVN)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic()
perc_sim_freq_out_of_bounds
ggsave("plots/percents_MVN_sims_need_truncation_low.png",
       plot = perc_sim_freq_out_of_bounds,
       height = 3, width = 6, units = "in")



MVNsim_zero_bounded <- data.frame(MVNsim_zero, stringsAsFactors = F) # sets bounds at 0 and 1
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

# save data objects from this simulation:
save(MVNsim_zero_bounded, meanA_MVNsim_zero_bounded, meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_CA_zero_bounded,
     file = "results/MVNsim_zero_bounded.RData")

# other simulations without covariances:
# compare distribution between binomial and MVN for 10^5 simulations
# first read in individual alpha estimates for mean A ancestry
CA_indAlpha <- do.call(rbind,
                       lapply(CA_pops_included$population, function(p) 
                         read.table(paste0("results/ancestry_hmm/combined_sept19/posterior/anc/", p, ".alpha.anc"),
                                    stringsAsFactors = F, header = T)))
AR_indAlpha <- do.call(rbind,
                       lapply(AR_pops_included$population, function(p) 
                         read.table(paste0("results/ancestry_hmm/combined_sept19/posterior/anc/", p, ".alpha.anc"),
                                    stringsAsFactors = F, header = T)))

set.seed(101) # use same seed
PoiBinsim_CA <- apply(do.call(rbind,
                              lapply(CA_indAlpha$A, 
                                     function(alpha)
                                       rbinom(n = n_sim, size = 2, prob = alpha)))/2, 2, mean)
PoiBinsim_AR <- apply(do.call(rbind,
                              lapply(AR_indAlpha$A, 
                                     function(alpha)
                                       rbinom(n = n_sim, size = 2, prob = alpha)))/2, 2, mean)
PoiBinsim_combined <- (PoiBinsim_CA*sum(CA_pops_included$n_bees) + PoiBinsim_AR*sum(AR_pops_included$n_bees))/sum(c(CA_pops_included$n_bees, AR_pops_included$n_bees))

# MVN simulation bounded with no covariances
set.seed(101)
MVNsim_no_cov <- mvrnorm(n = n_sim, 
                         mu = zAnc_bees$alpha, 
                         Sigma = diag(diag(zAnc_bees$K), length(zAnc_bees$alpha), length(zAnc_bees$alpha)), 
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
summary(diag(zAnc_bees$K))
summary(zAnc_bees$K[lower.tri(zAnc_bees$K, diag = F)])

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



                  
plot_QQ <- function(NA1, NA2, SA1, SA2, axis1, axis2, ps_qq = seq(0, 1, length.out = 10000)){
  d_qq <- bind_rows(data.frame(p = ps_qq,
                                           zone = "N. America", stringsAsFactors = F) %>%
                                  mutate(q1 = quantile(NA1, p), 
                                         q2 = quantile(NA2, p)),
                                data.frame(p = ps_qq,
                                           zone = "S. America", stringsAsFactors = F) %>%
                                  mutate(q1 = quantile(SA1, p), 
                                         q2 = quantile(SA2, p)))
  ggplot(data = d_qq) + 
    geom_point(aes(x = q1, y = q2, color = p)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab(axis1) +
    ylab(axis2) +
    scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
    theme_classic() +
    facet_grid(.~zone)#, scales = "free_x")
}
# make QQ plots for paper:
# effect of truncating at [0,1] is small:
MVN_vs_01_truncated_MVN_qq = plot_QQ(NA1 = meanA_MVNsim_CA_bounded, NA2 = meanA_MVNsim_CA_unbounded, 
        SA1 = meanA_MVNsim_AR_bounded, SA2 = meanA_MVNsim_AR_unbounded, 
        axis1 = "MVN simulation", axis2 = "MVN simulation unbounded", 
        ps_qq = seq(0, 1, length.out = 10000))
ggsave("plots/qq_effect_truncation_on_MVN.png",
       plot = MVN_vs_01_truncated_MVN_qq,
       height = 3, width = 6, 
       units = "in", device = "png")
effect_truncation_on_MVN <- arrangeGrob(perc_sim_freq_out_of_bounds + ggtitle("A"),
                                        MVN_vs_01_truncated_MVN_qq + ggtitle("B"),
                                         nrow = 2,
                                         ncol = 1)
ggsave("../../bee_manuscript/figures/effect_truncation_on_MVN.pdf", 
       plot = effect_truncation_on_MVN,
       height = 6.5, width = 6, 
       units = "in", device = "pdf")
# overall good fit of data to MVN:
mvn_vs_data_qq = plot_QQ(NA1 = meanA_MVNsim_CA_bounded, NA2 = meanA_CA, 
        SA1 = meanA_MVNsim_AR_bounded, SA2 = meanA_AR, 
        axis1 = "MVN simulation", axis2 = "Observed A frequencies", 
        ps_qq = seq(0, 1, length.out = 10000))
mvn_vs_data_qq
# very high expected false-positive rate with poisson-binomial null:
poibin_vs_data_qq = plot_QQ(NA1 = PoiBinsim_CA, NA2 = meanA_CA, 
        SA1 = PoiBinsim_AR, SA2 = meanA_AR, 
        axis1 = "Poisson binomial simulation", axis2 = "Observed A frequencies", 
        ps_qq = seq(0, 1, length.out = 10000))
# and high expected false-positive rate with no-covariance MVN:
mvn_no_cov_vs_data_qq = plot_QQ(NA1 = meanA_MVNsim_CA_no_cov_bounded, NA2 = meanA_CA, 
        SA1 = meanA_MVNsim_AR_no_cov_bounded, SA2 = meanA_AR, 
        axis1 = "MVN simulation with zero covariances", axis2 = "Observed A frequencies", 
        ps_qq = seq(0, 1, length.out = 10000))
# make joint plot for these QQs:
qq_vs_data_plots_combined <- arrangeGrob(mvn_vs_data_qq + ggtitle("A"),
                                         mvn_no_cov_vs_data_qq + ggtitle("B"),
                                         poibin_vs_data_qq + ggtitle("C"),
                                     nrow = 3,
                                     ncol = 1)
# plot them all together as a facet wrap instead:
sim_NA = list(meanA_MVNsim_CA_bounded, meanA_MVNsim_CA_no_cov_bounded, PoiBinsim_CA)
sim_SA = list(meanA_MVNsim_AR_bounded, meanA_MVNsim_AR_no_cov_bounded, PoiBinsim_AR)
sim_names = c("MVN", "MVN zero covariances", "Poisson binomial")
ps_qq = seq(0, 1, length.out = 10000)

d_qq <- do.call(rbind, lapply(1:3, function(i)
  bind_rows(data.frame(p = ps_qq,
             zone = "N. America",
             sim = sim_names[i],
             stringsAsFactors = F) %>%
    mutate(q2 = quantile(meanA_CA, p),
      q1 = quantile(sim_NA[[i]], p)),
    data.frame(p = ps_qq,
             zone = "S. America",
             sim = sim_names[i],
             stringsAsFactors = F) %>%
      mutate(q2 = quantile(meanA_AR, p), 
             q1 = quantile(sim_SA[[i]], p)))))
ggplot(data = d_qq) + 
    geom_point(aes(x = q1, y = q2, color = p), size = 0.5) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Simulated A frequencies") +
    ylab("Observed A frequencies") +
    scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
    theme_classic() +
    facet_grid(zone~sim)#, scales = "free_x")
ggsave("../../bee_manuscript/figures/qq_vs_data_sim_comparison.pdf",
       device = "pdf",
       width = 6, 
       height = 5, units = "in")
ggsave("plots/qq_vs_data_sim_comparison.png",
       device = "png",
       width = 6, 
       height = 5, units = "in")

library(plotly)
library(MASS)

# Compute kde2d
kd_data <- with(data = data.frame(x = meanA_CA, y = meanA_AR), 
           MASS::kde2d(y = y, x = x, n = 100,
                       lims = c(c(0,.6),c(0,.6))))
kd_mvn <- with(data = data.frame(x = meanA_MVNsim_CA_unbounded, 
                                 y = meanA_MVNsim_AR_unbounded), 
               MASS::kde2d(y = y, x = x, n = 100,
                           lims = c(c(0,.6),c(0,.6))))
ggplot(data.frame(kd_mvn)) +
  geom_density(aes(x = x, y = y, z = z))
# Plot with plotly
plot_ly(x = kd_data$x, y = kd_data$y, z = kd_data$z) %>% add_surface()
plot_ly(x = kd_mvn$x, y = kd_mvn$y, z = kd_mvn$z) %>% add_surface()
plot_ly(x = kd_data$x, 
        y = kd_data$y, 
        z = kd_data$z - kd_mvn$z,
        xlab = "N. America",
        ylab = "S. America") %>% 
  add_surface() %>% 
  layout(
    title = "2D Density Observed A frequency minus MVN",
    scene = list(
      xaxis = list(title = "NA", range = c(1, 0)),
      yaxis = list(title = "SA", range = c(1, 0)),
      zaxis = list(title = "Diff.")
    ))


# confirm identical 2d grid points for each density
identical(kd_data$x, kd_mvn$x) == T
identical(kd_data$y, kd_mvn$y) == T

# calc difference between density estimates
kd_diff = kd_data 
kd_diff$z = kd_data$z - kd_mvn$z

image(kd_diff, col = viridis(30))

# melt z data to long format from 100 x 100 matrix
rownames(kd_diff$z) = kd_diff$x
colnames(kd_diff$z) = kd_diff$y
kd_diff$z %>% # melt
  melt(., id.var = rownames(kd_diff)) %>%
  data.table::setnames(c("A_NA", "A_SA", "A_DIFF")) %>%
  ggplot(., aes(x = A_NA, A_SA, z = A_DIFF, fill = A_DIFF)) +
  stat_density2d(aes(alpha = A_DIFF),
                 geom="polygon", bins = 100)
  
  geom_tile() +
  stat_contour(aes(colour = ..level..), binwidth = .01) +
  scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
  scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"), midpoint=0) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  guides(colour=FALSE)

plot_2d_density <- function(){
  contour(kd_mvn, col = viridis(10)[1], nlevels = 7, 
          xlim = c(0.12, 0.35),
          ylim = c(0.32, 0.5),
          xlab = "N. America",
          ylab = "S. America")
  contour(kd_data, add = T, col = viridis(10)[9], nlevels = 7)
  legend("bottomright", legend = c("Observed A frequency", "Simulated A frequency (MVN)"),
         lwd = 1, lty = 1, col = viridis(10)[c(9,1)], cex = 0.5)
}
  
pdf(file = "../../bee_manuscript/figures/comparison_2d_density_data_mvn.pdf",
    height = 4, width = 6)
plot_2d_density()
dev.off()
png(file = "plots/comparison_2d_density_data_mvn.png",
    height = 4, width = 6, units = "in", res = 300)
plot_2d_density()
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
        meanA_MVNsim_CA_bounded > MVNsim_CA_bounded_quantile)/n_sim
.01^2 # it's about 5-6 times as likely under the MVN to get a double outlier at top .01% than it would be under true independence between N and S America
# but still same order of magnitude. Also I maybe need to simulate a larger # of trials to get accuracy in the tails (only expect 1000 outliers each out of 100k, very small # overlap)
# when I simulated 1 million MVN's I got .000521 as the probability of outliers in both CA and AR at .01%
MVNsim_AR_bounded_quantile_low <- quantile(meanA_MVNsim_AR_bounded, .01)
MVNsim_CA_bounded_quantile_low <- quantile(meanA_MVNsim_CA_bounded, .01)
table(meanA_MVNsim_AR_bounded < MVNsim_AR_bounded_quantile_low & 
        meanA_MVNsim_CA_bounded < MVNsim_CA_bounded_quantile_low)/n_sim
# on the low side it's similar. some variance from low # expected hits in my simulations


# these outliers are pretty surprising under a MVN framework
# how about under a poisson binomial framework? even more surprising

# compare simulations
# plot comparison
sim_compare <- 
  data.frame(poisson_binomial = PoiBinsim_combined,
           MVN_no_covariance = meanA_MVNsim_no_cov_bounded,
           #MVN_no_bounds = meanA_MVNsim,
           MVN_with_covariance = meanA_MVNsim_zero_bounded,
           #MVN_with_covariance_and_negs = meanA_MVNsim_bounded,
           #observed_data = meanA) %>%
           observed_data = sample(meanA, n_sim, replace = F)) %>% # downsample data to match length of simulations
  tidyr::gather(., "distribution", "A")
# plot
p_sim_compare <- sim_compare %>%
  ggplot(., aes(x = A, color = distribution)) +
  #geom_density_line(lwd = 1) +
  geom_line(stat = "density", #alpha = .75, 
            aes(linetype = distribution)) +
  theme_classic() +
  xlab("Combined sample mean African ancestry") +
  ylab("Density") +
  scale_color_viridis_d(option = "viridis", name = NULL, 
                        limits = c("observed_data", "MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_linetype_manual(name = NULL, values = c(1,2,3,4),
                        limits = c("observed_data", "MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
#scale_color_discrete(name = "Distribution", labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
p_sim_compare
p_sim_compare2 <- sim_compare %>%
  ggplot(., aes(x = A, color = distribution)) +
  #geom_density_line(lwd = 1) +
  geom_line(stat = "density", alpha = .75, lwd = 1) +
  theme_classic() +
  xlab("Mean African ancestry") +
  ylab("Density") +
  scale_color_viridis_d(option = "viridis", name = NULL, 
                        limits = c("observed_data", "MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("observed_data"="Observed data", "poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN"))
# use alpha instead of different line types
p_sim_compare2

ggsave("plots/distribution_data_vs_poibin_vs_MVN_sim.png",
       plot = p_sim_compare2,
       device = "png",
       width = 6, height = 4, units = "in")
ggsave("../../bee_manuscript/figures/distribution_data_vs_poibin_vs_MVN_sim.pdf",
       plot = p_sim_compare2,
       device = "pdf",
       width = 6, height = 4, units = "in")

# make joint plot of kinship matrix and these distribution comparisons:
dist_k_plots_combined <- arrangeGrob(k_plot_all + theme(axis.text = element_blank(), 
                                                         axis.ticks = element_blank()) + 
                                        ggtitle("A") +
                                        xlab("Ancestry correlation matrix"),
                                      p_sim_compare2 + ggtitle("B"),
                                      ncol = 2,
                                      widths = c(3, 5))
ggsave("../../bee_manuscript/figures/k_matrix_and_poi_bin_mvn_dist_comparison.pdf",
       plot = dist_k_plots_combined,
       device = "pdf",
       width = 8, 
       height = 3, units = "in")

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

# walk along sd's above mean ancestry for each zone
sd_CA <- sd(meanA_CA)
mu_CA <- mean(meanA_CA)
sd_AR <- sd(meanA_AR)
mu_AR <- mean(meanA_AR)
sd_range <- seq(0, 5, by = .005) # I can make this more precise later if I want

# get standard deviation cutoffs for FDR shared high loci

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
data.frame(mu = c(mu_CA, mu_AR), sd = c(sd_CA, sd_AR), 
           zone = c("N. America", "S. America"), short_name = c("CA", "AR"),
           stringsAsFactors = F) %>% 
  write.table(., "results/mu_sd_CA_AR.txt", quote = F, col.names = T, row.names = F, sep = "\t")

# what percent of the genome matches these cutoffs?
table(meanA_CA > FDRs$CA_high[FDRs$FDR_values==.05])/length(meanA_CA) # 0.26% of loci
table(sites$scaffold[meanA_CA > FDRs$CA_high[FDRs$FDR_values==.05]]) # where are they?

# high AR
table(meanA_AR > FDRs$AR_high[FDRs$FDR_values==.05])/length(meanA_AR) # 0.06% of loci (very few!)
table(sites$scaffold[meanA_AR > FDRs$AR_high[FDRs$FDR_values==.05]]) # where are they?

# about .3% of the genome is found in high shared sites at 5% FDR
table(meanA_CA > FDRs$shared_high_CA[FDRs$FDR_values==.05] & meanA_AR > FDRs$shared_high_AR[FDRs$FDR_values == .05])/length(meanA_CA)
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
  geom_density_2d(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
             aes(x = meanA_CA, y = meanA_AR), color = "blue") +
  geom_density_2d(data = data.frame(MVN_CA = meanA_MVNsim_CA_zero_bounded,
                                    MVN_AR = meanA_MVNsim_AR_zero_bounded), 
                  aes(x=MVN_CA, y=MVN_AR), color = "orange")
ggplot() +
  stat_density_2d(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
                  aes(x = meanA_CA, y = meanA_AR, fill = stat(level)), geom = "hex")

cor(meanA_CA, meanA_AR)
cor(meanA_MVNsim_AR_zero_bounded, meanA_MVNsim_CA_zero_bounded)
# I can use stat_contour if I have a z gridded for what I want to display


# make a new kind of density plot:
p_density <- ggplot(data = data.frame(CA = meanA_CA, AR = meanA_AR), #[c(T, rep(F, 10)), ],
       aes(x = CA, y = AR)) +
  #geom_point()
  #geom_bin2d(bins = 100) +
  #stat_density2d() +
  ggpointdensity::geom_pointdensity(adjust = 2, size = .5) + 
  # change method away from kde2d for original k neighbor plot, but may have to thin
  # takes way longer -- also need to rasterize
  scale_color_viridis(option = "magma") +
  scale_fill_viridis(option = "magma") +
  theme_classic() +
  xlab("N. America") +
  ylab("S. America") #+
  #labs(color = "No. Nearby Points")
plot(p_density +
       geom_density_2d(data = data.frame(MVN_CA = meanA_MVNsim_CA_zero_bounded,
                                         MVN_AR = meanA_MVNsim_AR_zero_bounded), 
                       aes(x=MVN_CA, y=MVN_AR), color = "orange"))

b <- as.mcmc(cbind(meanA_MVNsim_CA_zero_bounded, meanA_MVNsim_AR_zero_bounded))
a <- emdbook::HPDregionplot(b, vars = c("meanA_MVNsim_CA_zero_bounded", "meanA_MVNsim_AR_zero_bounded"), 
                            prob = c(0.99, .9, .75), 
                            n = 100, # number of grid points
                            #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
                            #col=c("lightgreen", "salmon", "lightblue", "yellow"), 
                            lwd = 3, 
                            h = .05,
                            add = TRUE)
b_data <- as.mcmc(cbind(meanA_CA, meanA_AR))
a_data <- emdbook::HPDregionplot(b_data, vars = c("meanA_CA", "meanA_AR"), 
                                 prob = c(0.99, .9, .75), 
                                 n = 100, # number of grid points
                                 #h = .1, let the function choose smoothing automatically (defaults to bandwidth.nrd())
                                 #col=c("lightgreen", "salmon", "lightblue", "yellow"), 
                                 lwd = 3, 
                                 h = .05,
                                 add = TRUE)
p_density2 <- ggplot() +
  geom_hex(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
           aes(x = meanA_CA, y = meanA_AR, color = ..count.., fill = ..count..), bins = 200) +
  scale_fill_viridis(option = "magma", end = 1, begin = 0.1) +
  scale_color_viridis(option = "magma", end = 1, begin = 0.1) +
  theme_classic() +
  xlab("N. America") +
  ylab("S. America") +
  geom_polygon(data = data.frame(a[[1]]), 
               aes(x = x, y = y),
               fill = NA, color = "orange", lwd = 1)
plot(p_density2)
ggsave("plots/A_CA_vs_AR_density.png",
       plot = p_density2,
       height = 5, width = 6, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/A_CA_vs_AR_density.pdf",
       plot = p_density2,
       height = 5, width = 6, units = "in", device = "pdf")
# thin points and plot nearest neighbors:
p_neighbors <- ggplot(data = data.frame(CA = meanA_CA, AR = meanA_AR)[c(T, rep(F, 100)), ],
                    aes(x = CA, y = AR)) +
  #geom_point()
  #geom_bin2d(bins = 100) +
  #stat_density2d() +
  ggpointdensity::geom_pointdensity(adjust = .1, size = .5) + 
  # change method away from kde2d for original k neighbor plot, but may have to thin
  # takes way longer -- also need to rasterize
  scale_color_viridis(option = "magma") +
  scale_fill_viridis(option = "magma") +
  theme_classic() +
  xlab("N. America") +
  ylab("S. America") #+
#labs(color = "No. Nearby Points")
plot(p_neighbors)
ggsave("plots/A_CA_vs_AR_neighbors.png",
       plot = p_density,
       height = 5, width = 6, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/A_CA_vs_AR_neighbors.png",
       height = 5, width = 6, units = "in", device = "png")



plot(meanA_CA, meanA_AR, col = alpha("grey", alpha = .3), pch = 20,
     main = "Ancestry frequencies in both hybrid zones across SNPs",
     xlab = "Mean African ancestry in California bees",
     ylab = "Mean African ancestry in Argentina bees")



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
sites_bed <- read.table("results/SNPs/combined_sept19/chr.var.sites.bed",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("scaffold", "start", "end", "snp_id"))
#table(paste(sites$scaffold, sites$pos, sep = "_") == sites_bed$snp) # match up, good
A_AR_CA <- sites %>%
  dplyr::select(chr, pos, cum_pos) %>%
  bind_cols(sites_bed, .) %>%
  mutate(CA = meanA_CA, AR = meanA_AR, combined = meanA) %>%
  # add FDRs
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
                FDR_shared_low, FDR_AR_low, FDR_CA_low, combined, chr, pos, cum_pos)
  

# write bed file with mean ancestry for Argentina and California included bees
# also include whether a site meets a FDR threshold for selection
A_AR_CA %>%
  write.table(., 
              "results/mean_ancestry_AR_CA.bed", 
              sep = "\t", quote = F, col.names = F, row.names = F)
save(A_AR_CA, file = "results/mean_ancestry_AR_CA.RData")

# I think I'm making some independence assumption here with the FDR rates
# It might be better to bin ancestry into broader windows, before calculating covariances.

# Also Harpur's study blasted the markers against the genome + 5kb to get coordinates for prior varroa-hygeine QTLs 

# plotting K covariance in ancestry matrix
melt(zAnc_bees$K) %>%
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
melt(cov2cor(zAnc_bees$K)) %>%
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

melt(zAnc_bees$K) %>% # filter out Avalon to get better color distinction between the others
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
melt(zAnc_bees$K) %>%
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
melt(cov2cor(zAnc_bees$K)) %>%
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
n_bees_per_pop <- sapply(names(zAnc_bees$alpha), function(p) meta.pop[meta.pop$population == p, "n_bees"])
sampling_var_diag <- ((2*n_bees_per_pop/(2*n_bees_per_pop - 1))*n_bees_per_pop*2*zAnc_bees$alpha*(1 - zAnc_bees$alpha))/(2*n_bees_per_pop)^2 
K_bees_minus_sampling <- zAnc_bees$K
for (i in 1:length(sampling_var_diag)){
  K_bees_minus_sampling[i,i] <- K_bees_minus_sampling[i,i] - sampling_var_diag[i]
}
K_bees_minus_sampling %>% # not quite
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


log(1/(1 - diag(zAnc_bees$K)/(zAnc_bees$alpha*(1 - zAnc_bees$alpha))))
log(diag(zAnc_bees$K)/(zAnc_bees$alpha*(1 - zAnc_bees$alpha)))

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
zAnc_bees$K
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
zAnc_bees$K
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
  label_short = c("High A: shared (5% FDR outlier)",
                  "High A: shared",
                  "Low A: shared",
                  "High A: S. America",
                  "Low A: S. America",
                  "High A: N. America"),
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
# plot all outliers in a grid
a %>%
  filter(snp_id %in% some_snp_outliers$snp_id) %>%
  left_join(., some_snp_outliers, by = "snp_id") %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
  geom_point() +
  geom_point(data = a_mean, shape = 2, color = "grey") +
  #ggtitle(paste0("cline across latitude for ", some_snp_outliers[i, "label"], ": ", some_snp_outliers[i, "snp_id"])) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  facet_grid(zone ~ label_short) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("plots/ancestry_clines_at_top_outlier_snps_all_pops_incl_mx.png",
       width = 8, height = 5, units = "in",
       device = "png")
# only plot the populations that are included (CA/AR 2014+2018),
# at the top outliers of each type:
a %>%
  filter(., population %in% c(AR_pops_included$population, CA_pops_included$population)) %>%
  filter(snp_id %in% some_snp_outliers$snp_id) %>%
  left_join(., some_snp_outliers, by = "snp_id") %>%
  filter(name != "a_shared_high_outlier_.05sig") %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = population)) +
  geom_point() +
  geom_point(data = filter(a_mean, population %in% c(AR_pops_included$population, CA_pops_included$population)),
             color = "grey", 
             shape = 2) +
  #ggtitle(paste0("cline across latitude for ", some_snp_outliers[i, "label"], ": ", some_snp_outliers[i, "snp_id"])) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  facet_grid(zone ~ label_short) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("plots/ancestry_clines_at_top_outlier_snps.png",
       width = 8, height = 5, units = "in",
       device = "png")
ggsave("../../bee_manuscript/figures/ancestry_clines_at_top_outlier_snps.png",
       width = 8, height = 5, units = "in",
       device = "png")

# new approach: get 1 SNP per outlier. So find top SNP in each outlier region:
# outlier_sets load outlier sets from rank_genes.R
# and outlier_setnames too
# here I have all the outliers of each type:
top_SNP_outliers_all
top_SNP_outliers
outlier_set_names
top_SNP_outlier_labels = c("High A: shared",
                "Low A: shared",
                "High A: S. America",
                "Low A: S. America",
                "High A: N. America")
for (i in 1:length(top_SNP_outlier_labels)){
a %>%
  filter(., population %in% c(AR_pops_included$population, CA_pops_included$population)) %>%
  filter(snp_id %in% top_SNP_outliers[[i]]$top_snp) %>%
  left_join(., top_SNP_outliers[[i]], by=c("snp_id"="top_snp", "scaffold")) %>%
    mutate(snp = paste(substr(scaffold, start = 6, stop = 100), start+1, sep = ":")) %>%
  ggplot(aes(x = abs(lat), y = A_ancestry, color = factor(FDR_region))) + #color = population, shape = factor(FDR_region))) +
  geom_point() +
  geom_point(data = filter(a_mean, population %in% c(AR_pops_included$population, CA_pops_included$population)),
             color = "grey", 
             shape = 2) +
  xlab("absolute latitude") +
  ylab("population mean African ancestry") +
  ggtitle(top_SNP_outlier_labels[i]) +
  facet_grid(zone ~ snp) +
  theme_bw() +
  #scale_shape_manual(values = c(20, 1))+
  scale_color_manual(values = c("red", "orange")) +
  theme(legend.position = "bottom") +
  #guides(color = F) +
  #labs(shape = "FDR") +
  labs(color = "FDR")
  ggsave(paste0("plots/ancestry_clines_at_all_outliers_top_snp_", outlier_set_names[i],".png"),
         width = 16, height = 5, units = "in",
         device = "png")
  ggsave(paste0("../../bee_manuscript/figures/ancestry_clines_at_all_outliers_top_snp_", outlier_set_names[i],".png"),
         width = 16, height = 5, units = "in",
         device = "png")
  }


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


