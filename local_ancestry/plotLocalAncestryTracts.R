# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(reshape2)

source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions


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


# Newest bee sequences
# load metadata
# load population ancestry frequencies:
pops <- read.table("../bee_samples_listed/byPop/pops_included.list", stringsAsFactors = F)$V1
popA <- lapply(pops, function(p) read.table(paste0("results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A <- do.call(cbind, popA)
colnames(A) <- pops
meanA <- apply(A, 1, mean)

sites0 <- read.table("results/SNPs/thin1kb_common3/included_scaffolds.pos", stringsAsFactors = F,
                    sep = "\t", header = F)
colnames(sites0) <- c("scaffold", "pos")
sites <- separate(sites0, scaffold, c("chr", "scaffold_n"), remove = F) %>%
  mutate(chr_n = substr(chr, 6, 100)) %>%
  mutate(snp_id = paste0("snp", chr_n, ".", scaffold, ".", pos))
#test_id <- read.table("results/SNPs/thin1kb_common3/included.snplist", stringsAsFactors = F,
#                      sep = "\t", header = F)$V1
#table(test_id == sites$snp_id) # good

meta <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  dplyr::select(c("population", "source", "year", "group")) %>%
  unique() %>%
  left_join(data.frame(population = pops, stringsAsFactors = F),
            ., by = "population")


# make plots - all loci, mean ancestry across all pops
png("plots/histogram_A_ancestry_all_loci.png", height = 6, width = 8, units = "in", res = 300)
hist(meanA, main = "all loci: mean ancestry across populations")
abline(v = max(meanA), col = "blue")
abline(v = min(meanA), col = "blue")
abline(v = quantile(meanA, .99), col = "orange")
abline(v = quantile(meanA, .01), col = "orange")
dev.off()


CA_A <- A[, meta$population[meta$group %in% c("N_CA", "S_CA", "CA_2018") & meta$year >= 2014]]
CA_A_earlier <- A[, meta$population[meta$group %in% c("N_CA", "S_CA") & meta$year < 2014]]
CA_A_2014 <- A[, meta$population[meta$group %in% c("N_CA", "S_CA") & meta$year == 2014]]
CA_A_2018 <- A[, meta$population[meta$group %in% c("CA_2018")]]
plot(apply(CA_A_2014, 1, mean), apply(CA_A_2018, 1, mean))
plot(apply(CA_A_earlier, 1, mean), apply(CA_A_2018, 1, mean))
AR_A <- A[, meta$population[meta$group == "AR_2018"]]
meanA_CA <- apply(CA_A, 1, mean)
meanA_AR <- apply(AR_A, 1, mean)

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

png("plots/k_matrix_all_pops.png", height = 6, width = 8, units = "in", res = 300)
melt(zAnc_bees$K) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees")
melt(cov2cor(zAnc_bees$K)) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees")
dev.off()

summary(meanA)
sd(meanA)

# simulate from MVN distribution:
AR_CA_K <- make_K_calcs(t(cbind(AR_A, CA_A)))
#sim_n <- 10^6 # slow
sim_n <- 10^5
MVNsim <- mvrnorm(n = sim_n, 
                  mu = AR_CA_K$alpha, 
                  Sigma = AR_CA_K$K, 
                  tol = 1e-6, 
                  empirical = FALSE, 
                  EISPACK = FALSE)
hist(apply(MVNsim, 1, mean))
summary(apply(MVNsim, 1, mean))
MVNsim_bounded <- MVNsim # sets bounds at 0 and 1
MVNsim_bounded[MVNsim < 0] <- 0
MVNsim_bounded[MVNsim > 1] <- 1
summary(apply(MVNsim_bounded, 1, mean)) # makes little difference
hist(apply(MVNsim_bounded, 1, mean))
MVNsim_AR <- MVNsim[ , colnames(AR_A)]
MVNsim_AR_bounded <- MVNsim_bounded[ , colnames(AR_A)]
MVNsim_CA <- MVNsim[ , colnames(CA_A)]
MVNsim_CA_bounded <- MVNsim_bounded[ , colnames(CA_A)]
plot(apply(MVNsim_AR_bounded, 1, mean), apply(MVNsim_CA_bounded, 1, mean))

png("plots/mean_A_ancestry_CA_vs_AR_grey.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = alpha("grey", .1),
     main = "comparison of A ancestry in California vs. Argentina zone")
dev.off()

png("plots/mean_A_ancestry_CA_vs_AR_MVNsim_outliers_blue.png", 
    height = 8, width = 12, res = 300, units = "in")
plot(meanA_CA, meanA_AR, pch = 20, col = ifelse((meanA_CA > quantile(apply(MVNsim_CA_bounded, 1, mean), .9999) | 
                                                   meanA_AR > quantile(apply(MVNsim_AR_bounded, 1, mean), .9999)) |
                                                  (meanA_CA < quantile(apply(MVNsim_CA_bounded, 1, mean), .0001) | 
                                                     meanA_AR < quantile(apply(MVNsim_AR_bounded, 1, mean), .0001)), 
                                                alpha("blue", .1), alpha("grey", .1)),
     main = "comparison of A ancestry in California vs. Argentina zone, blue = top 0.01% sim")
dev.off()
abline(v = quantile(apply(MVNsim_CA, 1, mean), .9999), col = "blue")
abline(h = quantile(apply(MVNsim_AR, 1, mean), .9999), col = "blue")
abline(v = quantile(apply(MVNsim_CA, 1, mean), .0001), col = "blue")
abline(h = quantile(apply(MVNsim_AR, 1, mean), .0001), col = "blue")
abline(v = quantile(apply(MVNsim_CA_bounded, 1, mean), .9999), col = "orange")
abline(h = quantile(apply(MVNsim_AR_bounded, 1, mean), .9999), col = "orange")
abline(v = quantile(apply(MVNsim_CA_bounded, 1, mean), .0001), col = "orange")
abline(h = quantile(apply(MVNsim_AR_bounded, 1, mean), .0001), col = "orange")

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
MVNsim_AR_bounded_quantile <- quantile(apply(MVNsim_AR_bounded, 1, mean), .99)
MVNsim_CA_bounded_quantile <- quantile(apply(MVNsim_CA_bounded, 1, mean), .99)
table(apply(MVNsim_AR_bounded, 1, mean) > MVNsim_AR_bounded_quantile & apply(MVNsim_CA_bounded, 1, mean) > MVNsim_CA_bounded_quantile)/sim_n
.01^2 # it's about six times as likely under the MVN to get a double outlier at top .01% than it would be under true independence between N and S America
# but still same order of magnitude. Also I maybe need to simulate a larger # of trials to get accuracy in the tails (only expect 1000 outliers each out of 100k, very small # overlap)
# when I simulated 1 million MVN's I got .000521 as the probability of outliers in both CA and AR at .01%
MVNsim_AR_bounded_quantile_low <- quantile(apply(MVNsim_AR_bounded, 1, mean), .01)
MVNsim_CA_bounded_quantile_low <- quantile(apply(MVNsim_CA_bounded, 1, mean), .01)
table(apply(MVNsim_AR_bounded, 1, mean) < MVNsim_AR_bounded_quantile_low & apply(MVNsim_CA_bounded, 1, mean) < MVNsim_CA_bounded_quantile_low)/sim_n
# on the low side it's similar, .00045, when I use 1 million simulations .. still only 450 observations


# these outliers are pretty surprising under a MVN framework
# how about under a poisson binomial framework? I should use means from the ancestry caller


melt(AR_CA_K$K) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  #scale_fill_gradient2(low = "white", mid = "grey", high = "darkblue") +
  ggtitle("K matrix - bees AR & CA")

r <- read.table("results/SNPs/thin1kb_common3/included_pos_on_Wallberg_Amel4.5_chr.r.bed",
                sep = "\t", stringsAsFactors = F) %>%
  mutate(snp_id = read.table("results/SNPs/thin1kb_common3/included_pos_on_Wallberg_Amel4.5_chr.map",
                             stringsAsFactors = F, header = F, sep = "\t")$V2) %>%
  rename(chr = V1) %>%
  rename(cM_Mb = V4) %>%
  mutate(pos_chr = V2 + 1) %>% # because .bed files are 0 indexed
  dplyr::select(c("snp_id", "chr", "pos_chr", "cM_Mb"))

d <- bind_cols(sites, A) %>%
  left_join(., r, by = c("chr", "snp_id"))
d$r_bin10 = with(d, cut(cM_Mb,  # note need to load map from scaffolds_to_chr.R
                             breaks = unique(quantile(c(0, map$r, 100), # I extend bounds so that every value gets in a bin
                                                      p = seq(0, 1, by = .1))),
                             right = T,
                             include.lowest = T))
table(d$r_bin10) # fewer SNPs represented in lowest recomb. bins
head(d[is.na(d$r_bin10),]) # good, should be none

d$r_bin5 = with(d, cut(cM_Mb,  # note need to load map from scaffolds_to_chr.R
                        breaks = unique(quantile(c(0, map$r, 100), # I extend bounds so that every value gets in a bin
                                                 p = seq(0, 1, by = .2))),
                        right = T,
                        include.lowest = T))
table(d$r_bin5) # fewer SNPs represented in lowest recomb. bins
head(d[is.na(d$r_bin5),]) # good, should be none


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

