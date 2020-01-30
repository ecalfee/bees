library(dplyr)
library(tidyr)
library(ggplot2)
source("../colors.R") # get color palette

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
load("../global_ancestry/results/NGSAdmix/ACM_K3_combined_sept19_chr_prunedBy250.rData")

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
  left_join(., acm.meta, by = c("pop_id"="CASEID"))
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
  filter(., Label != "RW-D-AR0605-01-25-20.JPG") %>% # duplicated pictures, just use other higher res one
  dplyr::select(Bee_ID, Label, pixels_per_cm) %>%
  write.table(., file = "results/bee_scale_1cm_rulers_1.28.20.txt", col.names = T, row.names = F, sep = "\t")
# note: there are more rulers than wing lengths because some wings were broken for full length but are otherwise ok for vein measurements etc.
  


# wing length cline across latitude
wings.meta %>%
  ggplot(., aes(x = abs(lat), y = wing_cm, color = geographic_location)) + 
  geom_point() +
  theme_classic()+
  geom_smooth() +
  geom_hline(data = acm.wing.lengths, aes(yintercept = mean_cm))

cor(wings.meta$wing_cm, wings.meta$A)
m_wing <- lm(wings.meta$wing_cm ~ wings.meta$A)
summary(m_wing)
# plot
# wing length predicted by mean genomewide ancestry
wings.meta %>%
  mutate(zone = ifelse(geographic_location == "Argentina", "S. America", "N. America")) %>%
  ggplot(., aes(x = A, y = wing_cm*10)) + 
  geom_point(data = acm.wings %>% # reference bees
               left_join(., data.frame(ACM = c("A", "C", "M"), 
                                       A = c(1, 0, 0), # African ancestry ref pops.
                                       by = "ACM", stringsAsFactors = F),
                         by = "ACM"), size = 1, alpha = 0.75, aes(color = ACM)) +
  #geom_smooth(data = acm.wings %>% # reference bees
  #             left_join(., data.frame(ACM = c("A", "C", "M"), 
  #                                     A = c(1, 0, 0), # African ancestry ref pops.
  #                                     by = "ACM", stringsAsFactors = F),
  #                       by = "ACM"), method = "lm", type = 2, color = "black") +
  geom_point(size = 1, alpha = .75, aes(color = zone)) + # hybrid zone bees
  theme_classic()+
  geom_abline(aes(intercept = 10*coefficients(m_wing)["(Intercept)"],
                  slope = 10*coefficients(m_wing)["wings.meta$A"])) +
  #ylim(c(0,10)) +
  #xlim(c(0,1)) +
  xlab("African ancestry proportion") +
  ylab("Wing length (mm)") +
  scale_color_manual(values = c(col_NA_SA_both, col_ACM), name = "")
  
ggsave("plots/wing_length_by_A_ancestry.png", height = 3, width = 6)
ggsave("../../bee_manuscript/figures/wing_length_by_A_ancestry.tiff", height = 3, width = 6)

# add model predictions and residuals to wings dataframe
wings.meta$model_predicted_cm <- stats::predict(m_wing)
wings.meta$model_residual_cm <- stats::residuals(m_wing)

d.wings <- wings.meta %>%
  arrange(abs(lat)) %>%
  dplyr::select(Bee_ID, wing_cm, model_predicted_cm, model_residual_cm)
save(list = "d.wings", file = "results/wing_fits.RData")


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
save(A_wings_wide_counts, file = "results/A_wings_wide_counts.RData")

# alt. initial regression model uses mean of local ancestry genomewide rather than global ancestry estimate;
# choice of global estimate or mean of local ancestry estimate makes no difference, which is expected:
A_wings_wide_counts <- as.matrix(A_wings_wide_counts)
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
                  summary(m)$coefficients[2, c(1,4)]
                }))

save(fits, file = "results/ancestry_mapping_wing_length.RData")
summary(fits[ , 2]) # lowest p-value is .001
hist(fits[ , 2])

# fit to counts
fits_counts <- t(sapply(1:nrow(A_wings_wide_counts), function(i){
  m <- lm(d.wings$model_residual_cm ~ A_wings_wide_counts[i, ])
  summary(m)$coefficients[2, c(1,4)]
}))

save(fits_counts, file = "results/ancestry_mapping_wing_length_counts.RData")
summary(fits_counts[ , 2])
hist(fits_counts[ , 2])

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

sites <- left_join(sites0, chr_lengths[ , c("chr", "scaffold", "chr_n", "chr_start")], by = "scaffold") %>%
  mutate(cum_pos = pos + chr_start)

cbind(sites, data.frame(fits_counts, stringsAsFactors = F) %>% 
        data.table::setnames(c("b", "p.value"))) %>%
  mutate(even_chr = (chr_n %% 2 == 1)) %>%
  #mutate(row_n = 1:nrow(.)) %>%  filter(-log10(p.value) > 1 | row_n %% 5 == 1) %>% # only plot every 10th point under 2 
  ggplot(., aes(x = cum_pos, y = -log10(p.value), color = even_chr)) +
  scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  geom_point(cex = .2) +
  theme_classic() +
  xlab("Chromosome") +
  ylab(expression('-log'[10]*'('*italic('P')*')')) +
  scale_color_manual(values = brewer.pal(4,"Paired")) +
  theme(legend.position = "none")
ggsave("plots/admixture_mapping_wing_length.png", height = 2, width = 6)
ggsave("../../bee_manuscript/figures/admixture_mapping_wing_length.tiff", height = 2, width = 6)

cbind(sites, data.frame(fits, stringsAsFactors = F) %>% 
        data.table::setnames(c("b", "p.value"))) %>%
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

