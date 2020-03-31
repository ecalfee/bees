library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
source("../colors.R") # get color palette
library(devtools)
#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("SNPRelate")
#BiocManager::install("SeqArray")
#devtools::install_github("kegrinde/STEAM")
library(STEAM)
PREFIX = "combined_sept19"
# bees
pops <- read.table("../bee_samples_listed/byPop/combined_sept19_pops.list",
                   stringsAsFactors = F)$V1
# ancestries
ACM <- c("A", "C", "M")
bees <- do.call(rbind, 
                lapply(pops, function(p) data.frame(Bee_ID = read.table(paste0("../bee_samples_listed/byPop/", p, ".list"),
                                                                        stringsAsFactors = F)$V1, population = p, stringsAsFactors = F)))

# metadata
meta.ind <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  left_join(bees, ., by = c("Bee_ID", "population")) 
# get admixture data
load("../global_ancestry/results/NGSAdmix/ACM_K3_combined_sept19_chr_prunedBy250.rData") #d_admix_ACM_labelled

# load measurements and metadata for reference bee populations:
acm.meta <- read.table("results/reference_samples_scut_car_lig_iber_mell.csv",
                       header = T, sep = ",", stringsAsFactors = F)
acm.measurements <- read.table("results/results_museum_bees_1_per_pop_wing_lengths_1.29.20.csv",
                               header = T, sep = ",", stringsAsFactors = F) %>%
  rename(length_pixels = Length) %>%
  mutate(id = stringr::str_replace(Label, ".bmp", "")) %>%
  mutate(type = ifelse(Label == "Masstab_1107-Sud-Afrika-21.bmp", "ruler", "forewing"))
acm.ruler <- filter(acm.measurements, type == "ruler") %>%
  mutate(pixels_per_cm = length_pixels)
subspecies_codes <- data.frame(subspecies_code = c(1, 2, 6, 9, 21), 
                               subspecies = c("mellifera", "iberiensis", "carnica", "ligustica", "scutellata"),
                               ACM_group = c("M", "M", "C", "C", "A"),
                               stringsAsFactors = F)
acm.wings <- filter(acm.measurements, type == "forewing") %>%
  separate(., col = "id", into = c("pop_id", "bee_id", "subspecies_code")) %>%
  mutate(pixels_per_cm = acm.ruler$pixels_per_cm,
         wing_cm = length_pixels/pixels_per_cm,
         pop_id = as.integer(pop_id),
         bee_id = as.integer(bee_id),
         subspecies_code = as.integer(subspecies_code)) %>%
  left_join(., subspecies_codes, by = "subspecies_code") %>%
  left_join(., acm.meta, by = c("pop_id"="CASEID")) %>%
  filter(!is.na(ACM))
table(acm.wings$ACM == acm.wings$ACM_group)

# groups separate by length_cm
acm.summary <- 
  acm.wings %>%
  group_by(ACM) %>%
  summarise(mean_cm = mean(wing_cm),
            n = n())
# plot reference bees by latitude
acm.wings %>%
  ggplot(aes(x = abs(BREITDEZ), y = wing_cm, color = ACM)) +
  geom_point()



# load unscaled wing measurements
measurements <- read.table("results/results_wings_rulers_1.28.20.csv", 
                    sep = ",", header = T) %>%
  mutate(id = stringr::str_replace(Label, "[.+]JPG", "")) %>%
         #id = stringr::str_replace(id)) %>%
  tidyr::separate(., id, c("wing", "letter", "Bee_ID", "month", "day", "year"), 
                  convert = F, sep = "-") %>%
  tidyr::unite(., "date", c("month", "day", "year"), sep = ".") %>%
  mutate(Bee_ID = stringr::str_replace(Bee_ID, "[_]", ""),
         type = ifelse(X <= 293, "forewing", "ruler")) %>%
  rename(length_pixels = Length) %>%
  filter(., Label != "UNK_RW-G-AR1103-05-30-19.JPG") # label error; it's CA1103 but I already have a good image of this bee wing
# extract ruler/scale measurements
rulers <- measurements %>%
  filter(type == "ruler") %>%
  mutate(length_cm = ifelse(length_pixels < 900, 0.5, 1), # all images have 0.5cm or 1cm length measurement
         pixels_per_cm = length_pixels/length_cm) %>%
  arrange(desc(X)) %>%
  filter(!duplicated(Label)) # take last measurement only -- some rulers were remeasured
rulers %>%
  ggplot(., aes(x = X, y = pixels_per_cm)) +
  geom_point()

wings <- measurements %>%
  filter(type == "forewing") %>%
  filter(!duplicated(Label)) %>%
  left_join(., rulers[ , c("Label", "pixels_per_cm")], by = "Label") %>%
  mutate(wing_cm = length_pixels/pixels_per_cm) # scale wing measurements
wings %>% # one bee was imaged twice
  filter(Bee_ID == "AR0605") # same bee diff images, .003cm diff. in measurement error
wings <- filter(wings, Label != "RW-D-AR0605-01-25-20.JPG") # just use the original higher res. image
rulers <- filter(rulers, Label != "RW-D-AR0605-01-25-20.JPG") # duplicated pictures, just use other higher res one
wings %>%
  ggplot(., aes(x = Bee_ID, y = wing_cm)) +
  geom_point()


# combine data
wings.meta <- wings %>%
  left_join(., meta.ind, by = "Bee_ID") %>%
  left_join(., d_admix_ACM_labelled[ , c("Bee_ID", ACM)], by = "Bee_ID")

# write data
write.table(wings[ , c("Bee_ID", "Label", "type", "wing_cm", "length_pixels", "pixels_per_cm", "wing", "date")],
            file = "results/bee_wing_lengths_cm_1.28.20.txt",
            col.names = T, row.names = F, sep = "\t")
rulers %>%
  dplyr::select(Bee_ID, Label, pixels_per_cm) %>%
  write.table(., file = "results/bee_scale_1cm_rulers_1.28.20.txt", col.names = T, row.names = F, sep = "\t")
# note: there are more rulers than wing lengths because some wings were broken for full length but are otherwise ok for vein measurements etc.
  


# wing length cline across latitude
wings.meta %>%
  ggplot(., aes(x = abs(lat), y = wing_cm, color = geographic_location)) + 
  geom_point() +
  theme_classic()+
  geom_smooth() +
  geom_hline(data = acm.summary, aes(yintercept = mean_cm))

cor(wings.meta$wing_cm, wings.meta$A)

m_wing <- with(wings.meta, lm(wing_cm ~ A))
summary(m_wing)
glance(m_wing)
tidy(m_wing)

m_wing_SA <- with(wings.meta %>%
                    mutate(SA = as.numeric(geographic_location == "Argentina")), 
                  lm(wing_cm ~ A*SA))
summary(m_wing_SA)
glance(m_wing_SA)
tidy(m_wing_SA)

m_wing2 <- with(d_A, lm(wing_cm ~ alpha)) # same as m_wing
m_wing_temp <- with(d_A, lm(wing_cm ~ alpha + AnnualMeanTemp))
m_wing_lat <- with(d_A, lm(wing_cm ~ alpha + lat))
m_wing_min_cold <- with(d_A, lm(wing_cm ~ alpha + MinTempColdestMonth))
m_wing_mean_cold <- with(d_A, lm(wing_cm ~ alpha + MeanTempColdestQuarter))
m_wing_precip <- with(d_A, lm(wing_cm ~ alpha + AnnualPrecip))
summary(m_wing_lat)
summary(m_wing_temp)
summary(m_wing_min_cold)
summary(m_wing_mean_cold)
summary(m_wing_precip)
glance(m_wing2)
glance(m_wing_temp)
# plot
# wing length predicted by mean genomewide ancestry
wings.meta %>%
  mutate(zone = ifelse(geographic_location == "Argentina", "S. America", "N. America")) %>%
  ggplot(., aes(x = A, y = wing_cm*10)) + 
  # background
  annotate("rect", xmin = 0, xmax = 1, # linear fit on range 0,1
                ymin = 7.5, ymax = 9.5, fill = "grey", alpha = .2) +
  
  #geom_smooth(data = acm.wings %>% # reference bees
  #             left_join(., data.frame(ACM = c("A", "C", "M"), 
  #                                     A = c(1, 0, 0), # African ancestry ref pops.
  #                                     by = "ACM", stringsAsFactors = F),
  #                       by = "ACM"), method = "lm", type = 2, color = "black") +
  geom_point(size = 1, alpha = .75, aes(color = zone)) + # hybrid zone bees
  theme_classic() +
  #geom_abline(aes(intercept = 10*coefficients(m_wing)["(Intercept)"],
  #                slope = 10*coefficients(m_wing)["wings.meta$A"])) +
  #geom_segment(aes(x = 0, xend = 1, # linear fit on range 0,1
  #                 y = 10*coefficients(m_wing)["(Intercept)"],
  #                 yend = 10*(coefficients(m_wing)["(Intercept)"] + 
  #                              coefficients(m_wing)["wings.meta$A"])),
  #             lwd = .1) +
  geom_line(data = data.frame(A = c(0, 1), # linear fit on range 0,1
                   wing_cm = c(coefficients(m_wing)["(Intercept)"],
                   coefficients(m_wing)["(Intercept)"] + 
                                coefficients(m_wing)["A"])),
               color = "black") +
  #ylim(c(0,10)) +
  #xlim(c(0,1)) +
  xlab("African ancestry proportion") +
  ylab("Wing length (mm)") +
  geom_jitter(data = acm.wings %>% # reference bees
               left_join(., data.frame(ACM = c("A", "C", "M"), 
                                       A = c(1.03, -0.03, -.03), # African ancestry ref pops.
                                       by = "ACM", stringsAsFactors = F),
                         by = "ACM"), width = 0.015, size = 1, alpha = 0.75, aes(color = ACM)) +
  
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "")
  
ggsave("plots/wing_length_by_A_ancestry.png", height = 4, width = 5.2)
ggsave("../../bee_manuscript/figures/wing_length_by_A_ancestry.png", height = 4, width = 5.2)
ggsave("../../bee_manuscript/figures_supp/wing_length_by_A_ancestry.tiff", height = 4, width = 5.2)


# add model predictions and residuals to wings dataframe
wings.meta$model_predicted_cm <- stats::predict(m_wing)
wings.meta$model_residual_cm <- stats::residuals(m_wing)

d.wings <- wings.meta %>%
  arrange(abs(lat)) %>%
  dplyr::select(Bee_ID, wing_cm, model_predicted_cm, model_residual_cm)
save(list = c("d.wings", "m_wing", "m_wing_SA"), file = "results/wing_fits.RData")


# get individual ancestry data for bees with wing data:
A_wings_wide <- do.call(cbind, 
                        lapply(d.wings$Bee_ID, function(id) 
  read.table(paste0("../local_ancestry/results/ancestry_hmm/combined_sept19/posterior/anc/", 
                    id, ".A.anc"), stringsAsFactors = F)))
colnames(A_wings_wide) <- d.wings$Bee_ID
A_wings_wide <- as.matrix(A_wings_wide)
A_wings_mean <- apply(A_wings_wide, 2, mean)
#dplyr::select(meta.ind$Bee_ID)
save(A_wings_wide, file = "results/A_wings_wide.RData")
# counts, rather than marginalizing over the posterior:

A_wings_wide_counts <- do.call(cbind, 
                        lapply(d.wings$Bee_ID, function(id) 
                          read.table(paste0("../local_ancestry/results/ancestry_hmm/combined_sept19/posterior/anc_count/", 
                                            id, ".A.count"), stringsAsFactors = F)))
colnames(A_wings_wide_counts) <- d.wings$Bee_ID
A_wings_wide_counts <- as.matrix(A_wings_wide_counts)
save(A_wings_wide_counts, file = "results/A_wings_wide_counts.RData")
load("results/A_wings_wide_counts.RData")
# alt. initial regression model uses mean of local ancestry genomewide rather than global ancestry estimate;
# choice of global estimate or mean of local ancestry estimate makes no difference, which is expected:

A_wings_counts_mean <- apply(A_wings_wide_counts, 2, mean)
cor(d.wings$wing_cm, A_wings_mean)
m_wing2 <- lm(d.wings$wing_cm ~ A_wings_mean)
summary(m_wing2)
summary(m_wing)
plot(d.wings$model_predicted_cm, predict(m_wing2))
abline(0,1,col="blue")
plot(d.wings$model_residual_cm, resid(m_wing2))
abline(0,1,col="blue") 


# admixture mapping of residuals after fitting prediction of wing length from global ancestry
fits <- t(sapply(1:nrow(A_wings_wide), function(i){
                  m <- lm(d.wings$model_residual_cm ~ A_wings_wide[i, ])
                  cf <- summary(m)$coefficients
                  results <- c(cf[2,], cf[1,])
                  names(results) <- c("b", "se", "t.value", "p.value",
                                      "intercept", "intercept_se", 
                                      "intercept_t.value", "intercept_p.value")
                  return(results)
                }))

save(fits, file = "results/ancestry_mapping_wing_length.RData")
summary(fits[ , "p.value"]) # lowest p-value is .001
hist(fits[ , "p.value"])
load("results/ancestry_mapping_wing_length.RData")

# a better way to do this is to fit to counts, because using the posterior MAP
# doesn't have undesirable shrinkage like using the posterior mean
# should be otherwise very similar
fits_counts <- t(sapply(1:nrow(A_wings_wide_counts), function(i){
  m <- lm(d.wings$model_residual_cm ~ A_wings_wide_counts[i, ])
  cf <- summary(m)$coefficients
  results <- c(cf[2,], cf[1,])
  names(results) <- c("b", "se", "t.value", "p.value",
                      "intercept", "intercept_se", 
                      "intercept_t.value", "intercept_p.value")
  return(results)
}))

save(fits_counts, file = "results/ancestry_mapping_wing_length_counts.RData")
summary(fits_counts[ , "p.value"])
hist(fits_counts[ , "p.value"])
load("results/ancestry_mapping_wing_length_counts.RData")

# plot results across genome
# get SNP sites where ancestry was called
sites0 <- read.table("../local_ancestry/results/SNPs/combined_sept19/chr.var.sites", stringsAsFactors = F,
                     sep = "\t", header = F)[ , 1:2]
colnames(sites0) <- c("scaffold", "pos")
chr_lengths <- cbind(read.table("../data/honeybee_genome/chr.names", stringsAsFactors = F),
                     read.table("../data/honeybee_genome/chr.lengths", stringsAsFactors = F)) %>%
  data.table::setnames(c("chr", "scaffold", "chr_length")) %>%
  mutate(chr_n = 1:16) %>%
  mutate(chr_end = cumsum(chr_length)) %>%
  mutate(chr_start = chr_end - chr_length) %>%
  mutate(chr_mid = (chr_start + chr_end)/2)

# recombination bin
# sites_rpos.RData created in K_by_recom_rate.R
load("../local_ancestry/results/sites_rpos.RData") # sites_rpos

sites <- left_join(sites0, chr_lengths[ , c("chr", "scaffold", "chr_n", "chr_start")], by = "scaffold") %>%
  mutate(cum_pos = pos + chr_start) %>%
  left_join(., sites_rpos[ , c("scaffold", "pos", "pos_cM", "cM_Mb", "r_bin5", "r_bin5_factor")],
            by = c("scaffold", "pos"))
sites_map <- sites[ , c("pos_cM", "chr_n")] %>%
  rename(cM = pos_cM, chr = chr_n)


# what is the significance threshold after
# correcting for multiple-testing?
length_genome_cm <- sites %>%
  group_by(chr_n) %>%
  summarise(length_cM = diff(range(pos_cM)))
L = sum(length_genome_cm$length_cM) # total length cM
delta = L/nrow(sites) # marker density (mean cM spacing between markers)
C = 16 # number of chromosomes
fits_counts <- data.frame(fits_counts, stringsAsFactors = F) %>%
  mutate(z = b/se)
hist(fits_counts$z)
summary(fits_counts$z)

# quick adjustment for multiple testing -- number of blocks:
# i.e. back of the envelope bonferroni correction:
n_blocks = 237*22/100*30 # 237Mb*22cM/Mb*1M/100cM*30gen
n_blocks = L/100*30 # Morgans * min number gen
.05/n_blocks # p-value needed to get alpha = 0.05 w/ bonferroni
p.adjust(p = 0.05, method = "bonferroni", n = 16)
p.adjust(p = .05/n_blocks, method = "bonferroni", n = n_blocks)


# use R package from Grinde 2019 to get threshold sig. p-value
pval_thresh <- STEAM::get_thresh_analytic(g = 30,
                           type = "pval",
                           map = sites_map,
                           alpha = 0.05)

# calculate on my own from eq. 2 Grinde 2019:
# checks out. 
mult_test_p <- function(C, L, delta, g, z, alpha){
  beta = 0.01*g
  v <- function(y){
    (2/y)*(pnorm(y/2) - 0.5)/
      ((y/2)*pnorm(y/2) + dnorm(y/2))
  }
  1 - exp(-2*C*(1 - pnorm(z)) - 
           2*beta*L*z*dnorm(z)*v(z*sqrt(2*beta*delta)))
  }
mult_test_p_diff <- function(C, L, delta, g, z, alpha){
  mult_test_p(C, L, delta, g, z) - alpha
}


z_thresh <- uniroot(mult_test_p_diff,
        interval = c(1.96,10), 
        g = 30, 
        delta = L/nrow(sites), 
        L = L, 
        C = 16, 
        alpha = 0.05)$root
2*pnorm(z_thresh, lower.tail = F)


# range of thresholds for a range of admixture generations g
gs <- c(22, 47.65, 60.4, 150) # min, median, mean, max
pval_threshs <- sapply(gs, function(g) STEAM::get_thresh_analytic(g = g,
                                          type = "pval",
                                          map = sites_map,
                                          alpha = 0.05))
pval_threshs
-log10(pval_threshs)
# show admixture mapping results with
# threshold based on 47.6 gen. admixture (median across pops)
cbind(sites, data.frame(fits_counts, stringsAsFactors = F)) %>%
  mutate(even_chr = (chr_n %% 2 == 1)) %>%
  #mutate(row_n = 1:nrow(.)) %>%  filter(-log10(p.value) > 1 | row_n %% 5 == 1) %>% # only plot every 10th point under 2 
  ggplot(., aes(x = cum_pos, y = -log10(p.value), color = even_chr)) +
  scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  geom_point(cex = .2) +
  theme_classic() +
  xlab("Chromosome") +
  ylab(expression('-log'[10]*'('*italic('P')*')')) +
  geom_hline(yintercept = -log10(pval_threshs[[2]]), 
             color = "red", lty = "dashed") +
  scale_color_manual(values = brewer.pal(4,"Paired")) +
  theme(legend.position = "none")
ggsave("plots/admixture_mapping_wing_length.png", 
       height = 3, width = 5.2)
ggsave("../../bee_manuscript/figures/admixture_mapping_wing_length.png", 
       height = 3, width = 5.2, dpi = 600)
ggsave("../../bee_manuscript/figures_supp/admixture_mapping_wing_length.tiff", 
       height = 3, width = 5.2, dpi = 600)

cbind(sites, data.frame(fits, stringsAsFactors = F)) %>%
  mutate(even_chr = (chr_n %% 2 == 1)) %>%
  #mutate(row_n = 1:nrow(.)) %>%  filter(-log10(p.value) > 2 | row_n %% 10 == 1) %>% # only plot every 10th point under 2 
  ggplot(., aes(x = cum_pos, y = -log10(p.value), color = even_chr), 
         size = .2) +
  scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  geom_point() +
  theme_classic() +
  xlab("Chromosome") +
  ylab(expression('-log'[10]*'('*italic('P')*')')) +
  scale_color_manual(values = brewer.pal(4,"Paired")) +
  theme(legend.position = "none")

# any overlap between admixture mapping outliers and ancestry outliers?
# no
load("../local_ancestry/results/mean_ancestry_AR_CA.RData")
table(A_AR_CA_fits$scaffold[fits_counts$p.value <= .01])
table(A_AR_CA$FDR_shared_high[fits_counts$p.value <= .01])
table(A_AR_CA$FDR_CA_high[fits_counts$p.value <= .01])
table(A_AR_CA$FDR_AR_high[fits_counts$p.value <= .01])
table(A_AR_CA$FDR_AR_low[fits_counts$p.value <= .01])
min(fits_counts$p.value[!is.na(A_AR_CA$FDR_shared_high)])
min(fits_counts$p.value[!is.na(A_AR_CA$FDR_AR_high)])
min(fits_counts$p.value[!is.na(A_AR_CA$FDR_CA_high)])
min(fits_counts$p.value[!is.na(A_AR_CA$FDR_AR_low)])
#plot(log10(fits_counts$p.value), clines$params[clines$params$term == "b", "estimate"])
#cor(log10(fits_counts$p.value), log10(clines$params[clines$params$term == "b", "p.value"]))     




# make output file of all samples for supplementary information:
german_names <- data.frame(LANDNAM = c("Italien", "Tanzania", "England", "Norwegen", "Rhodesien", 
                                       "Frankreich", "Kenya", "Osterreich", "JugoslawienSlow", 
                                       "JugoslawienSerb", "Rumanien", "Ungarn", "Burundi", 
                                       "Rusland", "JugoslawienKroa", "Spanien", "Sud-Afrika"), 
                           location = c("Italy", "Tanzania", "England", "Norway", "Rhodesia (Zimbabwe)", 
                                        "France", "Kenya", "Austria", "Yugoslavia (Slovenia)", 
                                        "Yugoslavia (Serbia)", "Romania", "Hungary", "Burundi", 
                                        "Russia", "Yugoslavia (Croatia)", "Spain", "South Africa"))
sources <- data.frame(source = c("Ramirez", "Sheppard", "Harpur", "Calfee"),
                      publication = c("Published in Cridland et al. 2018, NCBI PRJNA385500",
                                 "Published in Cridland et al. 2017, collected by Sheppard et al., NCBI PRJNA294105", 
                                 "Published in Harpur et al. 2014, NCBI PRJNA216922",
                                 "This study"))



enjambre <- data.frame(Bee_ID = 
                         c("AR1002", 
                           "AR1403", 
                           "AR1511", 
                           "AR2002", 
                           "AR2201",
                           "AR2512",
                           "AR2603",
                           "AR2913"),
                       collection_note =
                         c("at entrance of a feral colony in an old gas tank",
                           "foraging on water within 5m of a feral colony",
                           "on ground within 5m of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a utility pole",
                           "at entrance of a feral colony in a tree cavity"),
                       stringsAsFactors = F)

bees_all <- d_admix_ACM_labelled %>%
  dplyr::select(Bee_ID, geographic_location_short, population, Date, Time, enjambre, lat, long, est_coverage, A, C, M, source) %>%
  left_join(., full_join(wings, 
                         filter(measurements, type == "ruler") %>%
                           filter(!duplicated(Label)) %>%
                           filter(., Label != "RW-D-AR0605-01-25-20.JPG") %>%
                           dplyr::select("Bee_ID", "Label"), 
                         by = "Bee_ID") %>%
              mutate(Label = ifelse(is.na(Label.x), as.character(Label.y), as.character(Label.x))) %>%
              dplyr::select(., Bee_ID, Label, wing_cm), 
            by = "Bee_ID") %>%
  mutate(location = ifelse(source == "Ramirez", "California", 
                           ifelse(geographic_location_short == "Pag Island, Croatia", "Croatia", geographic_location_short))) %>%
  left_join(., sources, by = "source") %>%
  bind_rows(., acm.wings %>%
  mutate(Bee_ID = paste0("Pop", pop_id, "_Ind", bee_id)) %>% 
  rename(population = ACM) %>%
  left_join(., german_names, by = "LANDNAM") %>%
  dplyr::select(Bee_ID, Label, wing_cm, location, population) %>%
  mutate(publication = "Published in Ruttner 1988, Oberursel Collection, Germany")) %>%
  rename(wing_image_file = Label) %>%
  rename(wing_length_cm = wing_cm) %>%
  rename(collection_date = Date,
         collection_time = Time,
         feral_nest = enjambre) %>%
  mutate(collection_date = sapply(1:nrow(.), function(i)
                                  ifelse(location[i] == "Argentina", paste(strsplit(collection_date[i], 
                                                                          split = "[.]")[[1]][c(2,1,3)], 
                                                                          collapse = "/"),
                                  collection_date[i]))) %>% # fix dates reversed for Argentina samples
  left_join(., enjambre, by = "Bee_ID") %>%
  mutate(collection_note = ifelse(is.na(collection_note) & publication == "This study", 
                                   "foraging on vegetation", 
                                   collection_note)) %>%
  rename(., latitude = lat, longitude = long) %>%
  dplyr::select(Bee_ID, location, population, latitude, longitude, collection_date, collection_time, collection_note,
                feral_nest, publication, wing_image_file, wing_length_cm, est_coverage, A, C, M) %>%
  arrange(latitude, population, location) %>%
  arrange(publication == "Published in Ruttner 1988, Oberursel Collection, Germany",
          publication != "This study")
#View(bees_all)
write.table(bees_all, "../../bee_manuscript/files_supp/S1_Table_Sample_information.txt",
            col.names = T, row.names = F, sep = "\t", quote = F) # write out supplementary table with sample information
