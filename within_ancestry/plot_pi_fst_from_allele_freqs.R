library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
library(gridExtra)

# this script takes the population allele frequencies
# from combine_pop_allele_freqs.R
# and calculates pi and Fst between all populations

pops <- read.table("../bee_samples_listed/byPop/pops_included.list",
                   stringsAsFactors = F)$V1
ACM <- c("A", "C", "M")
ACM_pops <- c(ACM, pops)
ancestries <- c("AA", "CC", "MM")
ancestries_combined <- c(ancestries, "combined")

# get bee population meta data
bees <- do.call(rbind, 
                lapply(pops, function(p) data.frame(Bee_ID = read.table(paste0("../bee_samples_listed/byPop/", p, ".list"),
                                                                        stringsAsFactors = F)$V1, population = p, stringsAsFactors = F)))
meta.ind <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  left_join(bees, ., by = c("Bee_ID", "population")) 
meta.pop <- meta.ind %>%
  dplyr::select(c("population", "source", "year", "group", "lat", "long")) %>%
  dplyr::group_by(population, source, year, group) %>%
  dplyr::summarise(n_bees = n(),
                   lat = mean(lat),
                   long = mean(long)) %>%
  mutate(zone = ifelse(group == "AR_2018", "S. America", "N. America")) %>%
  left_join(data.frame(population = pops, stringsAsFactors = F),
            ., by = "population") %>%
  mutate(lat = ifelse(population == "Riverside_1999", .[.$population == "Riverside_2014", "lat"], lat))


# read in the SFS's that were successfully calculated (for comparison):
path_sfs_folded <- "results/folded_SFS_incorrect/non-outlier_regions/"
SFS <- lapply(pops, function(x) 
  read.table(paste0(path_sfs_folded, "/combined/", x, ".folded.sfs"),
             stringsAsFactors = F))
calc_pi_sfs_folded <- function(sfs){
  n <- 0:(length(sfs) - 1) # number of copies observed
  p <- n/((length(sfs) - 1)*2) # allele frequency
  f <- sfs/sum(sfs) # observation frequency of sfs bin
  het <- 2*p*(1-p)
  sum(het*f) # take weighted mean of pi across frequency bins
}

sfs1 <- unlist(read.table(paste0(path_sfs_folded, "/combined/AR01", ".folded.sfs"),
                   stringsAsFactors = F))
calc_pi_sfs_folded(sfs1)
lapply(SFS, calc_pi_sfs_folded)
thetas <- do.call(rbind,
                  lapply(ancestries_combined, function(a)
  data.frame(population = pops,
             ancestry = a,
             theta = unlist(lapply(pops, function(x)
    calc_pi_sfs_folded(read.table(paste0(path_sfs_folded, "/", a, "/", x, ".folded.sfs"),
                          stringsAsFactors = F)))))))
# something's not quite right. populations with less positions with data have higher diversity estimates:
plot(lapply(SFS, sum), thetas[thetas$ancestry=="combined", "theta"], main = "neg. relationship ANGSD estimated pi and number of sites in SFS",
     xlab = "total sites in SFS", ylab = "theta estimate (combined)")

thetas_ACM <- data.frame(population = ACM,
                               ancestry = "combined",
                               theta = unlist(lapply(ACM, function(x)
                                 calc_pi_sfs_folded(read.table(paste0(path_sfs_folded, "/", "combined", "/", x, ".folded.sfs"),
                                                               stringsAsFactors = F))))) %>%
  mutate(ancestry = paste0(population, population))

thetas_test <- data.frame(population = pops,
                      theta = unlist(lapply(SFS, calc_pi_sfs_folded)),
                     ancestry = "combined")

thetas %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = abs(lat), y = theta, color = ancestry, shape = factor(year))) +
  geom_point() +
  ggtitle("Allelic diversity at known SNPs within ancestries and combined across all ancestries") +
  xlab("Degrees latitude from the equator") +
  facet_grid(. ~ zone, scales = "free_x") +
  geom_abline(data = thetas_ACM, aes(intercept = theta, slope = 0, color = ancestry))
ggsave("plots/pi_by_latitude_from_folded_SFS.png", device = "png",
       width = 10, height = 5)
ggsave("../../bee_manuscript/figures/pi_by_latitude_from_folded_SFS.png", device = "png",
       width = 10, height = 5)



# read in 1/10th of the allele freq data:
n_snp = 10
freqs <- 1 - read.table(paste0("results/allele_freq_all/pops_included_plus_ACM_every", n_snp, "th_SNP.freqs.txt"),
                    stringsAsFactors = F, header = T)
# Q for later -- if these are all sites with SNPs (MAF > .05), why do about one third not have MAF > .05 in any of the included pops?
# these are the SNPs from Julie's set; I'll get a fresh set of SNP calls on HAv3.1 and redo this plot
table(apply(freqs, 1, max, na.rm = T) >= .05)
table(apply(freqs, 1, max, na.rm = T) >= .03)

het_mean <- apply(hets, 2, function(x) mean(x, na.rm = T))
table(complete.cases(hets)) # mostly no NAs
# how many alleles sampled at a site?
ns <- read.table(paste0("results/allele_freq_all/pops_included_plus_ACM_every", n_snp, "th_SNP.nInd"),
                 stringsAsFactors = F, header = T)*2 # x2 because diploid

# freqs and heterozygosity within AA CC MM ancestry
freqs_by_ancestry <- lapply(ancestries, function(a) 
  cbind(freqs[ , ACM], 
        1 - read.table(paste0("results/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/", 
                              a, "/allele_freq/pops_included_every", n_snp, "th_SNP.freqs.txt"),
             stringsAsFactors = F, header = T)))
hets_by_ancestry <- lapply(freqs_by_ancestry, function(f) 2*f*(1-f))
het_mean_by_ancestry <- lapply(hets_by_ancestry, function(h)  # filter for MAF >= .05 for at least one population
  apply(h, 2, function(x) mean(x, na.rm = T)))

# all AA heterozygosities are below ref bee pop A (makes sense):
hist(het_mean_by_ancestry[[1]][pops], xlim = c(0, 0.1),
     main = "histogram of pi within AA homozygous ancestry regions",
     xlab = "population pi for AA ancestry")
abline(v = het_mean["A"], col = "blue")


# how many alleles were sampled at each site?
ns_by_ancestry <- lapply(ancestries, function(a) 
  cbind(ns[ , ACM], read.table(paste0("results/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/", 
                                      a, "/allele_freq/pops_included_every", n_snp, "th_SNP.nInd"),
                                      stringsAsFactors = F, header = T)*2))


d_het <- data.frame(population = ACM_pops,
                    AA = het_mean_by_ancestry[[1]],
                    CC = het_mean_by_ancestry[[2]],
                    MM = het_mean_by_ancestry[[3]],
                    combined = het_mean,
                    ref_pop = c(rep(T, 3), rep(F, length(pops)))) %>%
  left_join(., meta.pop, by = "population")
ACM_het <- data.frame(ancestry = ancestries,
                      combined = het_mean[ACM],
                      ref_pop = T)

# plot heterozygosity overall and within ancestry as you move across the hybrid zone:
d_het %>%
  filter(!ref_pop) %>%
  tidyr::gather(., "ancestry", "pi", c(ancestries, "combined")) %>%
  ggplot(., aes(x = abs(lat), y = pi, color = ancestry, shape = factor(year))) +
  geom_point() +
  facet_wrap(~zone, scales = "free_x") +
  geom_abline(data = ACM_het, aes(intercept = combined, slope = 0, color = ancestry))
table(is.na(freqs_AA[, "AR01"]))

# heterozygosity corrected for small sample size:
het_small_sample_correction <- function(p, n, filter_under_2 = T){
  ifelse(filter_under_2 & n <= 2, NA, 2*(p - p^2*(n/(n-1)) + p*(1/(n-1))))}
hets_small_sample <- do.call(cbind, lapply(1:ncol(freqs), 
                                    function(i) het_small_sample_correction(p = freqs[, i],
                                                                                           n = ns[, i])))
colnames(hets_small_sample) <- ACM_pops

het_small_sample_mean <- apply(hets_small_sample, 2, function(x) mean(x, na.rm = T))

hets_small_sample_by_ancestry <- lapply(1:3, function(a) do.call(cbind, 
                                                                 lapply(1:ncol(freqs_by_ancestry[[a]]), 
                                                                        function(i) het_small_sample_correction(p = freqs_by_ancestry[[a]][ , i], 
                                                                                                                n = ns_by_ancestry[[a]][ , i])))) 
hets_small_sample_by_ancestry2 <- lapply(1:3, function(a) do.call(cbind, 
                                                                 lapply(1:ncol(freqs_by_ancestry[[a]]), 
                                                                        function(i) het_small_sample_correction(p = freqs_by_ancestry[[a]][ , i], 
                                                                                                                n = ns_by_ancestry[[a]][ , i],
                                                                                                                filter_under_2 = F)))) 

het_small_sample_mean_by_ancestry <- lapply(hets_small_sample_by_ancestry, function(h) 
  apply(h, 2, function(x) mean(x, na.rm = T)))
d_het_small_sample <- data.frame(population = ACM_pops,
                                 AA = het_small_sample_mean_by_ancestry[[1]],
                                 CC = het_small_sample_mean_by_ancestry[[2]],
                                 MM = het_small_sample_mean_by_ancestry[[3]],
                                 combined = het_small_sample_mean,
                                 ref_pop = c(rep(T, 3), rep(F, length(pops)))) %>%
  left_join(., meta.pop, by = "population")
ACM_het_small_sample <- data.frame(ancestry = ancestries,
                                   combined = het_small_sample_mean[ACM],
                                   ref_pop = T)

# what do we expect pi to be for these admixed populations?
# first I need admixture data for each population:
# get admixture data
prefix <- "CA_AR_MX_harpur_sheppard_kohn_wallberg"
# get ID's for PCA data (CAUTION - bam list order and admix results MUST MATCH!)
IDs <- read.table(paste0("../bee_samples_listed/", prefix, ".list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
K = 3 # 3 admixing populations
n = 250 # snps thinned to 1 every nth
prefix1 = paste0("ordered_scaffolds_", prefix, "_prunedBy", n)
name = paste0("K", K, "_", prefix1)
file = paste0("../global_ancestry/results/NGSAdmix/", name, ".qopt")
admix <- read.table(file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2", "anc3)
# get meta data for all individuals included in NGSadmix analysis (plus extras)
ids.pops <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                       header = T, sep = "\t") %>%
  dplyr::select(Bee_ID, population)
# label ancestries
anc_labels <- data.frame(ancestry = colnames(admix),
                         ancestry_label = sapply(colnames(admix), 
                                                 function(x) names(which.max(tapply(admix2[ , x], admix2$population, sum)))),
                         stringsAsFactors = F)
admix2 <- bind_cols(IDs, admix) %>%
  left_join(ids.pops, by = "Bee_ID") %>%
  data.table::setnames(c("Bee_ID", anc_labels$ancestry_label, "population"))
admix.pops <- group_by(admix2, population) %>%
  summarise(A = mean(A), M = mean(M), C = mean(C), n = n())
# what should the allele freqs be based on the reference pop allele freqs?
# I need to fix code below:
freqs_predicted <- lapply(pops, function(p) t(sapply(1:nrow(freqs[1:10,]), function(i)
  freqs[i, "A"]*admix.pops[admix.pops$population==p, "A"] +
    freqs[i, "M"]*admix.pops[admix.pops$population==p, "M"] +
    freqs[i, "C"]*admix.pops[admix.pops$population==p, "C"])))

# then I need an allele freq. estimate for each ancestry .. may not match reference panels freqs.



d_het_small_sample %>%
  filter(!ref_pop) %>%
  tidyr::gather(., "ancestry", "pi", c(ancestries, "combined")) %>%
  ggplot(., aes(x = abs(lat), y = pi, color = ancestry, shape = factor(year))) +
  geom_point() +
  xlab("Degrees latitude from the equator") +
  facet_grid(zone ~ ., scales = "free_x") +
  geom_abline(data = ACM_het_small_sample, aes(intercept = combined, slope = 0, color = ancestry))
ggsave("plots/pi_by_latitude_including_mexico.png", device = "png",
       width = 10, height = 5)
d_het_small_sample %>%
  filter(!ref_pop) %>%
  filter(population != "MX10") %>%
  tidyr::gather(., "ancestry", "pi", c(ancestries, "combined")) %>%
  ggplot(., aes(x = abs(lat), y = pi, color = ancestry, shape = factor(year))) +
  geom_point() +
  ggtitle("Allelic diversity at known SNPs within ancestries and combined across all ancestries") +
  xlab("Degrees latitude from the equator") +
  facet_grid(. ~ zone, scales = "free_x") +
  geom_abline(data = ACM_het_small_sample, aes(intercept = combined, slope = 0, color = ancestry))
ggsave("plots/pi_by_latitude.png", device = "png",
       width = 10, height = 5)
ggsave("../../bee_manuscript/figures/pi_by_latitude.png", device = "png",
       width = 10, height = 5)

# no real relationship with pop coverage = good sign
plot(lapply(SFS, sum), het_small_sample_mean[4:length(het_small_sample_mean)], main = "neg. relationship ANGSD estimated pi and number of sites in SFS",
     xlab = "total sites in SFS", ylab = "theta estimate (small sample corr.)")

  # plot Fst -- populations further apart should have higher Fst
calc_Fst <- function(v1, v2){
  exclude = is.na(v1) | is.na(v2)
  v_1 = v1[!exclude]
  v_2 = v2[!exclude]
  het_1 = 2*v_1*(1-v_1)
  het_2 = 2*v_2*(1-v_2)
  v_tot = (v_1 + v_2)/2
  het_tot = 2*v_tot*(1-v_tot)
  fst = 1 - (mean(het_1)+mean(het_2))/2/mean(het_tot)
  return(fst)
}
calc_dxy <- function(v1, v2){
  exclude = is.na(v1) | is.na(v2)
  v_1 = v1[!exclude]
  v_2 = v2[!exclude]
  dxy = mean(v_1*(1-v_2) + v_2*(1-v_1))
  return(dxy) 
}
# Hudson pairwise Fst estimator from Bhatia 2012 
# (I implement equation 10 from the supplement & take the ratio of the avg. numerator & denominator) 
calc_hudson_Fst <- function(v1, v2, n1, n2){ # allele frequencies and number of alleles sampled
  exclude = is.na(v1) | is.na(v2)
  p_1 = v1[!exclude]
  p_2 = v2[!exclude]
  n_1 = n1[!exclude]
  n_2 = n2[!exclude]

  N = (p_1 - p_2)^2 - p_1*(1 - p_1)/(n_1 - 1) - p_2*(1-p_2)/(n_2 - 1)
  D = p_1*(1 - p_2) + (1 - p_1)*p_2
  
  fst = mean(N)/mean(D)
  return(fst)
}

# sort populations by latitude
pops_order_lat <- meta.pop$population[order(meta.pop$lat)]
ACM_pops_order_lat <- c(ACM, pops_order_lat)

# get fst pairwise between all pops
fst_matrix <- matrix(0, length(ACM_pops), length(ACM_pops))
for (i in 1:length(ACM_pops)){
  for (j in 1:length(ACM_pops)){
    fst_matrix[i, j] <- calc_Fst(v1 = freqs[ , ACM_pops[i]], v2 = freqs[ , ACM_pops[j]])
  }
}
colnames(fst_matrix) <- ACM_pops
rownames(fst_matrix) <- ACM_pops


reshape2::melt(fst_matrix[ACM_pops_order_lat, ACM_pops_order_lat]) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 1, end = 0, direction = 1) +  
  #scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("Fst - bees")

# get fst by ancestry
# make all within ancestry plots:
fst_matrix_by_ancestry <- lapply(freqs_by_ancestry, function(a){

  fst_matrix_a <- matrix(0, length(ACM_pops), length(ACM_pops))
  for (i in 1:length(ACM_pops)){
    for (j in 1:length(ACM_pops)){
      fst_matrix_a[i, j] <- calc_Fst(v1 = a[ , ACM_pops[i]], v2 = a[ , ACM_pops[j]])
    }
  }
  colnames(fst_matrix_a) <- ACM_pops
  rownames(fst_matrix_a) <- ACM_pops
  return(fst_matrix_a)
  }
)
plots_fst_by_ancestry <- lapply(1:3, function(i) 
  reshape2::melt(fst_matrix_by_ancestry[[i]][ACM_pops_order_lat, ACM_pops_order_lat]) %>%
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
    scale_fill_viridis(begin = 1, end = 0, direction = 1) +  
    ggtitle(paste0("Fst ", ancestries[[i]], " ancestry - bees")))

plots_fst_by_ancestry


# Hudson Fst estimator:
# blind to ancestry:
# hudson estimator of Fst
fst_hudson_matrix <- matrix(0, length(ACM_pops), length(ACM_pops))
for (i in 1:length(ACM_pops)){
  for (j in 1:length(ACM_pops)){
    fst_hudson_matrix[i, j] <- calc_hudson_Fst(v1 = freqs[ , ACM_pops[i]], 
                                               v2 = freqs[ , ACM_pops[j]],
                                               n1 = ns[ , ACM_pops[i]],
                                               n2 = ns[ , ACM_pops[j]])
  }
}
colnames(fst_hudson_matrix) <- ACM_pops
rownames(fst_hudson_matrix) <- ACM_pops
# plot hudson Fst
plot_fst_hudson <- reshape2::melt(fst_hudson_matrix[ACM_pops_order_lat, ACM_pops_order_lat]) %>%
  filter(! Var1 == Var2) %>% # don't plot diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab("") +
  #scale_fill_viridis(begin = 1, end = 0, direction = 1) +  
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("Hudson Fst - bees")
plot_fst_hudson


# for each ancestry:
fst_hudson_matrix_by_ancestry <- lapply(1:3, function(a){
  fst_hudson_matrix_a <- matrix(0, length(ACM_pops), length(ACM_pops))
  for (i in 1:length(ACM_pops)){
    for (j in 1:length(ACM_pops)){
      fst_hudson_matrix_a[i, j] <- calc_hudson_Fst(v1 = freqs_by_ancestry[[a]][ , ACM_pops[i]], 
                                            v2 = freqs_by_ancestry[[a]][ , ACM_pops[j]],
                                            n1 = ns_by_ancestry[[a]][ , ACM_pops[i]],
                                            n2 = ns_by_ancestry[[a]][ , ACM_pops[j]]
                                            )
    }
  }
  colnames(fst_hudson_matrix_a) <- ACM_pops
  rownames(fst_hudson_matrix_a) <- ACM_pops
  return(fst_hudson_matrix_a)
}
)

plots_fst_hudson_by_ancestry <- lapply(1:3, function(i) 
  reshape2::melt(fst_hudson_matrix_by_ancestry[[i]][ACM_pops_order_lat, ACM_pops_order_lat]) %>%
    filter(! Var1 == Var2) %>% # don't plot diagonal
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab(NULL) +
    ylab(NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
    #scale_fill_viridis(begin = 1, end = 0, direction = 1) +  
    ggtitle(paste0("Hudson Fst ", ancestries[[i]], " ancestry - bees")))

plots_fst_hudson_by_ancestry

plots_fst_hudson_4 <- grid.arrange(plot_fst_hudson,
             plots_fst_hudson_by_ancestry[[1]],
             plots_fst_hudson_by_ancestry[[2]],
             plots_fst_hudson_by_ancestry[[3]],
             nrow = 2, 
             ncol = 2)
ggsave("plots/fst_overall_and_within_ancestry.png",
       plot = plots_fst_hudson_4,
       device = "png",
       height = 14,
       width = 14,
       units = "in")
ggsave("../../bee_manuscript/figures/fst_overall_and_within_ancestry.png",
       plot = plots_fst_hudson_4,
       device = "png",
       height = 14,
       width = 14,
       units = "in")

fst_matrix_by_ancestry[[1]]["AR01", "CA14"]
fst_matrix_by_ancestry[[1]]["AR01", "AR03"]
fst_matrix_by_ancestry[[1]]["CA12", "CA14"]
fst_matrix_by_ancestry[[1]]["CA12", "AR03"]
fst_matrix_by_ancestry[[1]]["CA03", "AR03"]

# plot dxy too
# combined across all ancestries
dxy_matrix <- matrix(0, length(ACM_pops), length(ACM_pops))
for (i in 1:length(ACM_pops)){
  for (j in 1:length(ACM_pops)){
    dxy_matrix[i, j] <- calc_dxy(v1 = freqs[ , ACM_pops[i]], v2 = freqs[ , ACM_pops[j]])
  }
}
colnames(dxy_matrix) <- ACM_pops
rownames(dxy_matrix) <- ACM_pops
plots_dxy <- reshape2::melt(dxy_matrix[ACM_pops_order_lat, ACM_pops_order_lat]) %>%
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
    ggtitle("dxy combined - bees")
plots_dxy

# within ancestry
dxy_matrix_by_ancestry <- lapply(freqs_by_ancestry, function(a){
  
  dxy_matrix_a <- matrix(0, length(ACM_pops), length(ACM_pops))
  for (i in 1:length(ACM_pops)){
    for (j in 1:length(ACM_pops)){
      dxy_matrix_a[i, j] <- calc_dxy(v1 = a[ , ACM_pops[i]], v2 = a[ , ACM_pops[j]])
    }
  }
  colnames(dxy_matrix_a) <- ACM_pops
  rownames(dxy_matrix_a) <- ACM_pops
  return(dxy_matrix_a)
}
)
plots_dxy_by_ancestry <- lapply(1:3, function(i) 
  reshape2::melt(dxy_matrix_by_ancestry[[i]][ACM_pops_order_lat, ACM_pops_order_lat]) %>%
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
    ggtitle(paste0("dxy ", ancestries[[i]], " ancestry - bees")))

plots_dxy_4 <- grid.arrange(plots_dxy,
  plots_dxy_by_ancestry[[1]],
             plots_dxy_by_ancestry[[2]],
             plots_dxy_by_ancestry[[3]],
             nrow = 2, 
             ncol = 2)
ggsave("plots/dxy_overall_and_within_ancestry.png",
       plot = plots_dxy_4,
       device = "png",
       height = 14,
       width = 14,
       units = "in")
ggsave("../../bee_manuscript/figures/dxy_overall_and_within_ancestry.png",
       plot = plots_dxy_4,
       device = "png",
       height = 14,
       width = 14,
       units = "in")

# can I plot ancestry heterozygosity too?

# TO DO: for outlier regions, get pi and Fst across the outliers for each ancestry,
# group samples: N. America (just include 2014), S. America, A, C, M. 
# Step 1: get list of outlier regions in a format ANGSD understands.
# I will also need some neutral outlier regions to get backgroun mean pi and Fst.
# Also plot PCA for within-ancestry diversity.

# look at pi for some outliers:
# this is just a CA selected high A outlier:
pi1_AA_NA <- read.table("results/outlier_regions/region_1/AA/NA_3.thetas.windows.pestPG",
                  header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("window", "chr", "pos", "theta_Wattersons", "theta_pairwise", "theta_FuLi", "theta_FayH", "theta_L", "tajima_D", "FuLi_F", "FuLi_D", "FayH", "ZengE", "n_effective_sites"))
pi1_NA <- read.table("results/outlier_regions/region_1/combined/NA_3.thetas.windows.pestPG",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("window", "chr", "pos", "theta_Wattersons", "theta_pairwise", "theta_FuLi", "theta_FayH", "theta_L", "tajima_D", "FuLi_F", "FuLi_D", "FayH", "ZengE", "n_effective_sites"))


pi1_AA_SA <- read.table("results/outlier_regions/region_1/AA/SA_3.thetas.windows.pestPG",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("window", "chr", "pos", "theta_Wattersons", "theta_pairwise", "theta_FuLi", "theta_FayH", "theta_L", "tajima_D", "FuLi_F", "FuLi_D", "FayH", "ZengE", "n_effective_sites"))
pi1_A <- read.table("results/outlier_regions/region_1/combined/A.thetas.windows.pestPG",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("window", "chr", "pos", "theta_Wattersons", "theta_pairwise", "theta_FuLi", "theta_FayH", "theta_L", "tajima_D", "FuLi_F", "FuLi_D", "FayH", "ZengE", "n_effective_sites"))


# undexpected that diversity would be higher in the CA samples
# except if they're not truly all A
with(pi1_A, plot(pos, theta_Wattersons/n_effective_sites, col = "blue",
                 ylim = c(0, 0.1)))
with(pi1_AA_NA, points(pos, theta_Wattersons/n_effective_sites))
with(pi1_NA, points(pos, theta_Wattersons/n_effective_sites, col = "orange"))


fstA1 <- read.table("results/outlier_regions/region_1/combined/A-NA_3.fst.windows",
                    stringsAsFactors = F, sep = "\t")
plot(fstA1$Nsites)

