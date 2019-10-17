library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
library(rethinking)
source("../colors.R") # for color palette
source("/media/erin/3TB/Documents/gitErin/covAncestry/forqs_sim/k_matrix.R") # import useful functions

# this script looks at the effect of recombination rate on the K matrix
# in an effort to disentangle demography and selection
# in creating ancestry covariances between populations, especially the negative covariances observed

# I need to load the local ancestry data from plotLocalAncestryTracts.R:
# objects: A, meta.ind, meta.pop, sites

# load recombination rates:
rmap <- read.table("../data/recomb_map/Wallberg_HAv3.1/map_rates_extended_10kb.bed",
                   header = T,
                   stringsAsFactors = F)

# get genomewide r quintiles:
rmap$r_bin5 <- cut(rmap$cM_Mb,
                   breaks = unique(quantile(c(0, rmap$cM_Mb, max(rmap$cM_Mb)),
                                            p = seq(0, 1, by = .2))),
                   right = T,
                   include.lowest = T)
table(rmap$r_bin5) # good, each bin includes 20% of windows
levels(rmap$r_bin5)

# map recombination rates (and bins/quintiles) onto sites
sites_r <- bedr(
  engine = "bedtools", 
  input = list(a = dplyr::mutate(sites, start = pos - 1, end = pos) %>%  # make 0 index for bedtools
                 dplyr::select(chr, start, end, pos, scaffold, chr_n, chr_start, cum_pos),
               b = rmap), 
  method = "map", 
  # cM_Mb is col 4 and r_bin5 is col 6
  params = "-g ../data/honeybee_genome/chr_names.lengths -c 4,6 -o collapse", 
  check.chr = F
) 
colnames(sites_r) <- c("chr", "start", "end", "pos", "scaffold", 
                       "chr_n", "chr_start", "cum_pos", "cM_Mb", "r_bin5")
sites_r$cM_Mb <- as.numeric(sites_r$cM_Mb)
sites_r$pos <- as.numeric(sites_r$pos)
sites_r$chr_n <- as.numeric(sites_r$chr_n)
sites_r$cum_pos <- as.numeric(sites_r$cum_pos)
sites_r$chr_start <- as.numeric(sites_r$chr_start)
sites_r$r_bin5_factor <- factor(sites_r$r_bin5, levels = levels(rmap$r_bin5),
                                ordered = T)
sites_r$r_bin5 <- as.numeric(sites_r$r_bin5_factor) # translate to 1-5 numbers
# negative correlation between being a tight snp and high recombination rate, even when we exclude Group11
cor(A_clines1$b_lat[sites_r$chr] <= quantile(A_clines1$b_lat, .01), sites_r$cM_Mb)
cor(A_clines1$b_lat[sites_r$chr != "Group11"] <= quantile(A_clines1$b_lat, .01), sites_r$cM_Mb[sites_r$chr != "Group11"])

str(sites_r)
table(sites_r$r_bin5_factor)
table(is.na(sites_r$r_bin5_factor))
levels(sites_r$r_bin5_factor)
save(sites_r, file = "results/sites_r.RData")

# look at mean A ancestry across recombination bins:
cbind(sites_r, meanA) %>%
  group_by(r_bin5_factor) %>%
  summarise(mean = mean(meanA))

# get mean A ancestry for each population for the 5 recombination bins
A_r <- cbind(sites_r, A) %>%
  tidyr::gather(., "population", "A", colnames(A)) %>%
  group_by(r_bin5, r_bin5_factor, population) %>%
  summarise(A = mean(A)) %>%
  data.frame() %>%
  left_join(., meta.pop, by = "population") %>%
  mutate(abs_lat_c = abs(lat) - mean(abs(lat)))
A_r$r_bin5_number = as.numeric(A_r$r_bin5_factor)
A_r$r_bin5_number_c = A_r$r_bin5_number - mean(A_r$r_bin5_number)
unique(A_r[ , c("r_bin5", "r_bin5_number")]) # in order low to high, good
str(A_r)
A_r %>%
  filter(r_bin5_factor %in% c("[0,2.92]", "(38.6,100]")) %>%
  ggplot(., aes(x = abs(lat), y = A, color = r_bin5_factor)) +
  geom_point() +
  facet_wrap(~zone)

# I test and find that mean cline isn't steeper for low recombination rate regions:
# nor does mean A change with recombination rate bin
m_lat_r_mean_cline <- map2stan( # quadratic approximation of the posterior MAP
  alist(
    A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat[r_bin5_number]*abs_lat_c + b_r*r_bin5_number_c,
    mu ~ dnorm(0, 5),
    theta ~ dcauchy(0, 2), # allows it to sample large #s as needed
    b_lat[r_bin5_number] ~ dnorm(0, 5),
    #b_lat[r_bin5_number] ~ dnorm(b_lat_mean, b_lat_theta), # hierarchical model for effect of latitude
    #b_lat_mean ~ dnorm(0, 5),
    #b_lat_theta ~ dcauchy(0, 2),
    b_r ~ dnorm(0,5)
  ),
  data = A_r,
  iter = 2000, warmup = 1000, chains = 2, cores = 2)
write.table(precis(m_lat_r_mean_cline, depth = 2)@output,
            "results/m_lat_r_mean_cline_model_summary.txt",
            row.names = T, col.names = T)
precis(m_lat_r_mean_cline, depth = 2)
pairs(m_lat_r_mean_cline)
plot(m_lat_r_mean_cline) # there are warnings from STAN about 1 divergent transition, but chains look ok



# do we see a significant effect of recombination rate on the shape of the cline?
# short answer: no
d_A_r <- cbind(sites_r, A) %>%
  tidyr::gather(., "population", "A", colnames(A)) %>%
  left_join(., meta.pop, by = "population") %>%
  mutate(abs_lat_c = abs(lat) - mean(abs(lat))) %>%
  mutate(cM_Mb_c = cM_Mb - mean(cM_Mb))
dim(d_A_r) # very large dataset, I'll do the fit with a smaller portion
d_A_r_test <- d_A_r[c(T, rep(F, 100)),]
dim(d_A_r_test)
table(d_A_r_test$population)

m_lat_r_bypop <- map( # quadratic approximation of the posterior MAP
  alist(
    A ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- mu + b_lat*abs_lat_c + b_r*cM_Mb_c + b_r_lat*cM_Mb_c*abs_lat_c,
    mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10),
    b_lat ~ dnorm(0, 5),
    b_r ~ dnorm(0,5),
    b_r_lat ~ dnorm(0, 5)
  ),
  data = d_A_r_test,
  start = list(mu = -1, theta = 2, b_lat = 0, b_r = 0, b_r_lat = 0))
precis(m_lat_r_bypop)
pairs(m_lat_r_bypop)


# make K matrix for each recombination rate bin
K_byr <- lapply(1:5, function(r)
  make_K_calcs(t(A[sites_r$r_bin5 == r, pops_by_lat])))

# plot K for each recombination rate bin
for (a in 1:5){
  melt(K_byr[[a]]$K) %>%
    filter(Var1 != Var2) %>% # don't plot diagonal
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1,
                       limits = c(-.005, .01)) +
    ggtitle(paste0("K matrix covariances rbin ", a, ": ", levels(rmap$r_bin5)[a]))
  ggsave(paste0("plots/k_matrix_covariance_r", a, ".png"), 
         height = 6, width = 8, 
         units = "in", device = "png")
  melt(cov2cor(K_byr[[a]]$K)) %>%
    filter(Var1 != Var2) %>% # don't plot diagonal
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1,
                       limits = c(-.35, .55)) +
    ggtitle(paste0("K matrix correlations rbin ", a, ": ", levels(rmap$r_bin5)[a]))
  ggsave(paste0("plots/k_matrix_correlation_r", a, ".png"), 
         height = 6, width = 8, 
         units = "in", device = "png")
}
# try plotting these again but using the genomewide alpha rather than alpha
# for that recombination rate bin
pop_alpha_genomewide <- apply(A[ , pops_by_lat], 2, mean)
K_byr2 <- lapply(1:5, function(r)
  calcK(ancFreqMatrix = t(A[sites_r$r_bin5 == r, pops_by_lat]),
        alpha = pop_alpha_genomewide))
# make new plots:
for (a in 1:5){
  melt(K_byr2[[a]]) %>%
    filter(Var1 != Var2) %>% # don't plot diagonal
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1,
                       limits = c(-.005, .01)) +
    ggtitle(paste0("K matrix covariances rbin ", a, ": ", levels(rmap$r_bin5)[a], " genomewide alpha"))
  ggsave(paste0("plots/k_matrix_covariance_r", a, "_genomewide_alpha.png"), 
         height = 6, width = 8, 
         units = "in", device = "png")
  melt(cov2cor(K_byr2[[a]])) %>%
    filter(Var1 != Var2) %>% # don't plot diagonal
    ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_viridis(begin = 0, end = 1, direction = 1,
                       limits = c(-.35, .55)) +
    ggtitle(paste0("K matrix correlations rbin ", a, ": ", levels(rmap$r_bin5)[a], " genomewide alpha"))
  ggsave(paste0("plots/k_matrix_correlation_r", a, "_genomewide_alpha.png"), 
         height = 6, width = 8, 
         units = "in", device = "png")
}
# hmm.. make a side-by-side plot for the genomewide ones. Ask Graham how to compare. 
# Take a stab at interpretation and add it to your manuscript first. Also add whatever you have
# on how the slopes don't appear to get steeper from low to high recombination.
# also try excluding + outliers again. And if you have it, your steepest sloped loci, to see if that's driving
# the effect -- mean slope doesn't get steeper, but top 1% are in low r regions

# calculate mean correlations for different recombination rate bins:
# function from plotLocalAncestryTracts.R:
pair_types <- data.frame(label = c("Warm SA", "Cold SA vs. Warm SA", "Cold SA", "Cold NA vs. Warm SA", "Cold NA vs. Cold SA", "Cold NA"),
                         type = unique(mean_corrs$type)[order(unique(mean_corrs$type))],
                         stringsAsFactors = F) %>%
  mutate(label = factor(label, levels = c("Warm SA", "Cold SA", "Cold NA", "Cold NA vs. Cold SA", "Cold SA vs. Warm SA", "Cold NA vs. Warm SA"),
                        ordered = T))
mean_corrs <- do.call(rbind,
                      lapply(1:5, function(k) get_mean_corr_from_K(K_byr[[k]]$K) %>%
                       mutate(r_bin5 = k))) %>%
  left_join(., distinct(dplyr::select(sites_r, c("r_bin5", "r_bin5_factor"))), by = "r_bin5") %>%
  left_join(., pair_types, by = "type")
mean_corrs %>%
  write.table(., "results/mean_anc_corr_grouped_by_r.txt",
              col.names = T, row.names = F, quote = F, sep = "\t")
mean_corrs %>%
  ggplot(., aes(x = label, y = mean_anc_corr, color = r_bin5_factor)) +
  geom_point(alpha = .75) +
  geom_point(data = left_join(mean_corr_k, pair_types, by = "type"), 
             aes(x = label, y = mean_anc_corr), color = "black", pch = 4) +
  ylab("Mean Ancestry Correlation") +
  theme_classic() +
  xlab("") +
  scale_color_viridis_d(name = "Recombination bin (cM/Mb)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("plots/mean_k_corr_by_groups_and_r.png"), 
       height = 4, width = 6, 
       units = "in", device = "png")
ggsave(paste0("../../bee_manuscript/figures/mean_k_corr_by_groups_and_r.pdf"), 
       height = 4, width = 6, 
       units = "in", device = "pdf")

