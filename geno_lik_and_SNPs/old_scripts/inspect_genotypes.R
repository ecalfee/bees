# inspect genotype calls from all SNPs for A/C/M bees:
library(viridis)
library(viridisLite)
library(dplyr)
library(tidyr)
library(ggplot2)

# start by looking at chr1:
ACM <- c("A", "C", "M")
sites <- read.table("results/combined_sept19/variant_sites/Group1.rpos",
                    header = F, stringsAsFactors = F, sep = "\t") %>%
  data.table::setnames(c("scaffold", "pos", "major", "minor", "rpos"))
ACM_GL <- lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/genotypes/", a, "/Group1.mafs.gz"),
                header = T, stringsAsFactors = F, sep = "\t") %>%
    dplyr::select(chromo, position, knownEM, nInd) %>%
    data.table::setnames(c("scaffold", "pos", paste0("freq_", a), paste0("nInd_", a))) %>%
    left_join(sites, ., by = c("scaffold", "pos")))
GL <- cbind(ACM_GL[[1]], ACM_GL[[2]][ , 6:7], ACM_GL[[3]][ , 6:7])

ACM_counts <- lapply(ACM, function(a) 
  read.table(paste0("../within_ancestry/results/combined_sept19/allele_freq/combined/Group1/",
                    a, ".mafs.gz"),
                      header = T, stringsAsFactors = F, sep = "\t") %>%
    dplyr::select(chromo, position, phat, nInd) %>%
    data.table::setnames(c("scaffold", "pos", paste0("freq_", a), paste0("nInd_", a))) %>%
    left_join(sites, ., by = c("scaffold", "pos")))
counts <- cbind(ACM_counts[[1]], ACM_counts[[2]][ , 6:7], ACM_counts[[3]][ , 6:7])

ACM_genos0 <- lapply(ACM, function(a) 
  read.table(paste0("results/combined_sept19/genotypes/", a, "/Group1.geno.gz"),
             header = F, stringsAsFactors = F, sep = "\t", na.strings = "-1") %>%
    rename(scaffold = V1, pos = V2) %>%
    left_join(sites[ , c("scaffold", "pos")], ., by = c("scaffold", "pos"))) # so doesn't skip SNPs with no genotypes called
ACM_genos1 <- lapply(ACM_genos0, function(acm_genos){ # count allele from genotype calls
  acm_genos[ , c("scaffold", "pos")] %>%
    mutate(n_tot = apply(acm_genos[ , -c(1,2)], 
                       1, 
                       function(row) 2*(sum(!is.na(row))))) %>%
    mutate(n_minor = apply(acm_genos[ , -c(1,2)], 1, function(row) 
      sum(row, na.rm = T))) %>%
    mutate(n_major = n_tot - n_minor)
})
CMA_genos <- data.frame(sites, 
                        C_major = ACM_genos1[[2]][ , "n_major"],
                        C_minor = ACM_genos1[[2]][ , "n_minor"],
                        M_major = ACM_genos1[[3]][ , "n_major"],
                        M_minor = ACM_genos1[[3]][ , "n_minor"],
                        A_major = ACM_genos1[[1]][ , "n_major"],
                        A_minor = ACM_genos1[[1]][ , "n_minor"])


ggplot(data = data.frame(GL = GL$freq_A, counts = counts$freq_A, nInd = GL$nInd_A),
       aes(x = counts, y = GL, color = nInd)) +
  geom_point() +
  scale_color_viridis(option = "viridis")
cor(GL$freq_A, counts$freq_A, use = "complete.obs")
cor(GL$freq_C, counts$freq_C, use = "complete.obs")
cor(GL$freq_M, counts$freq_M, use = "complete.obs")
cor(GL$freq_M, GL$freq_A, use = "complete.obs")
cor(GL$freq_C, GL$freq_A, use = "complete.obs")
cor(GL$freq_M, GL$freq_C, use = "complete.obs")
ggplot(data = data.frame(A = GL$freq_A, 
                         C = GL$freq_C, 
                         minInd = sapply(1:nrow(GL), function(i) min(GL$nInd_A[i], GL$nInd_C[i]))),
       aes(x = A, y = C, color = minInd)) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  ggtitle("A vs. C")
ggplot(data = data.frame(A = GL$freq_A, 
                         M = GL$freq_M, 
                         minInd = sapply(1:nrow(GL), function(i) min(GL$nInd_A[i], GL$nInd_M[i]))),
       aes(x = A, y = M, color = minInd)) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  ggtitle("A vs. M")
mean(GL$freq_A, na.rm = T)
mean(GL$freq_C, na.rm = T)
mean(GL$freq_M, na.rm = T)

# what is the approximate discriminatory power if I take a basic approach and 
# start by filtering SNPs for # individuals with data and 1 pop with > 0.3 MAF:
GL1 <- filter(GL, complete.cases(GL)) %>%
  mutate(freq_tot = (freq_A + freq_C + freq_M)/3) %>%
  mutate(het_tot = 2*freq_tot*(1-freq_tot)) %>%
  filter(het_tot != 0) %>%
  mutate(fst_A = 1 - (2 * freq_A * (1 - freq_A))/het_tot,
         fst_C = 1 - (2 * freq_C * (1 - freq_C))/het_tot,
         fst_M = 1 - (2 * freq_M * (1 - freq_M))/het_tot)
GL2 <- GL1 %>%
  filter(., nInd_A >= 6, nInd_C >= 6, nInd_M >= 6)
GL3 <- GL2 %>%
  filter(., (freq_A >= .3 | freq_C >= .3 | freq_M >= .3) &
           (freq_A <= .7 | freq_C <= .7 | freq_M <= .7))
GL4 <- GL2 %>%
  filter(., fst_A > 0)
GL5 <- GL2 %>%
  filter(., abs(freq_A - freq_C) >= 0.3 | abs(freq_A - freq_M) >= 0.3)
GL6 <- GL2 %>%
  filter(., abs(freq_A - (freq_C + freq_M)/2) >= 0.3)
GL7 <- GL2 %>%
  filter(., (freq_A >= .4 | freq_C >= .4 | freq_M >= .4) &
           (freq_A <= .6 | freq_C <= .6 | freq_M <= .6))
GL8 <- GL2 %>%
  filter(., fst_A > -.1)
GLs <- list(GL1, GL2, GL3, GL4, GL5, GL6, GL7, GL8)
lapply(c(GL, GLs), dim)
lapply(GLs, function(x) mean(x$fst_A))
lapply(GLs, function(x) mean(x$fst_C))
lapply(GLs, function(x) mean(x$fst_M))
lapply(GLs, function(x) mean(x$nInd_A))
lapply(GLs, function(x) mean(x$nInd_C))
lapply(GLs, function(x) mean(x$nInd_M))

table(GL$nInd_C >= 8 & GL$nInd_A >= 8 & GL$nInd_M >= 8)/nrow(GL)
lapply(1:16, function(i) cor(GL$freq_A[GL$nInd_A == i], 
                             counts$freq_A[GL$nInd_A == i], 
                             use = "complete.obs"))

# If I use GL7, at least one pop with freq >= 0.4, what is the spacing like 
# between snps?
GL7_spacing <- diff(GL7$rpos)
summary(GL7_spacing)
summary(diff(GL2$rpos))
summary(diff(GL3$rpos))
min(diff(GL3$rpos))
# all GL3 sites:
GL3_CMA_genos <- left_join(GL3[ , c("scaffold", "pos")],
                           CMA_genos, by = c("scaffold", "pos")) %>%
  mutate(rdiff = c(1, format(round(diff(rpos)/100, 12), scientific = F))) # 12 decimals, units Morgans, no scientific number
# rdiff starts at 1 for chr position 1 (there is not difference in Morgans to calc)

# get counts from external files:
ids <- read.table("../local_ancestry/results/ancestry_hmm/AR11.ploidy", header = F,
                  sep = "\t", stringsAsFactors = F)$V1
# each id has a major and minor allele counts:
ids_major_minor <- unlist(lapply(ids, function(i) paste(i, c("major", "minor"), sep = "_")))
ids_counts <- cbind(GL3[ , c("scaffold", "pos")],
                    data.frame(do.call(cbind, 
                                       lapply(ids, function(i) 
                                         read.table(paste0("../local_ancestry/results/SNPs/TEST/countsMajMin/", 
                                                           i, ".counts.txt"), 
                                                    header = F, sep = "\t", stringsAsFactors = F)))) %>%
                      data.table::setnames(ids_major_minor))

# allele counts for ref A/C/M pops for ancestry hmm
write.table(GL3_CMA_genos[ , c("scaffold", "pos", 
                               "C_major", "C_minor", 
                               "M_major", "M_minor", 
                               "A_major", "A_minor", "rdiff")], 
            "../local_ancestry/results/SNPs/TEST/Group1_highLD_CMA.genos",
            col.names = F, row.names = F, quote = F, sep = "\t")
# variant sites file
write.table(GL3_CMA_genos[ , 1:4], "../local_ancestry/results/SNPs/TEST/Group1_highLD_CMA.var.sites",
            col.names = F, row.names = F, quote = F, sep = "\t")
# ancestry hmm input file for pop AR11
write.table(left_join(GL3_CMA_genos[ , c("scaffold", "pos", 
                               "C_major", "C_minor", 
                               "M_major", "M_minor", 
                               "A_major", "A_minor", "rdiff")],
                      ids_counts, by = c("scaffold", "pos")),
            "../local_ancestry/results/ancestry_hmm/TEST/high_LD/Group1_highLD_CMA.counts",
            col.names = F, row.names = F, quote = F, sep = "\t")




# thin SNPs:
spacing = .005 # min cM spacing
keep_snp <- rep(F, nrow(GL3))
last_rpos <- -100
for (i in 1:nrow(GL3)) {
  this_rpos <- GL3[i, "rpos"]
  if (this_rpos - last_rpos >= spacing){
    keep_snp[i] <- T
    last_rpos <- this_rpos 
  } # else do nothing  
}
GL3_CMA_genos_low <- left_join(GL3[keep_snp, c("scaffold", "pos")],
                           CMA_genos, by = c("scaffold", "pos")) %>%
  mutate(rdiff = c(1, format(round(diff(rpos)/100, 12), scientific = F))) # 12 decimals, units Morgans, no scientific number
write.table(GL3_CMA_genos_low[ , c("scaffold", "pos", 
                               "C_major", "C_minor", 
                               "M_major", "M_minor", 
                               "A_major", "A_minor", "rdiff")], 
            "../local_ancestry/results/SNPs/TEST/Group1_lowLD_CMA.genos",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(GL3_CMA_genos_low[ , 1:4], "../local_ancestry/results/SNPs/TEST/Group1_lowLD_CMA.var.sites",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(left_join(GL3_CMA_genos_low[ , c("scaffold", "pos", 
                                         "C_major", "C_minor", 
                                         "M_major", "M_minor", 
                                         "A_major", "A_minor", "rdiff")],
                      ids_counts, by = c("scaffold", "pos")),
            "../local_ancestry/results/ancestry_hmm/TEST/low_LD/Group1_lowLD_CMA.counts",
            col.names = F, row.names = F, quote = F, sep = "\t")



# more moderate LD filtering to .001cM
# thin SNPs:
spacing2 = .002 # min cM spacing
keep_snp2 <- rep(F, nrow(GL3))
last_rpos2 <- -100
for (i in 1:nrow(GL3)) {
  this_rpos2 <- GL3[i, "rpos"]
  if (this_rpos2 - last_rpos2 >= spacing2){
    keep_snp2[i] <- T
    last_rpos2 <- this_rpos2 
  } # else do nothing  
}
GL3_CMA_genos_med <- left_join(GL3[keep_snp2, c("scaffold", "pos")],
                               CMA_genos, by = c("scaffold", "pos")) %>%
  mutate(rdiff = c(1, format(round(diff(rpos)/100, 12), scientific = F))) # 12 decimals, units Morgans, no scientific number

write.table(GL3_CMA_genos_med[ , c("scaffold", "pos", 
                                   "C_major", "C_minor", 
                                   "M_major", "M_minor", 
                                   "A_major", "A_minor", "rdiff")], 
            "../local_ancestry/results/SNPs/TEST/Group1_medLD_CMA.genos",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(GL3_CMA_genos_med[ , 1:4], "../local_ancestry/results/SNPs/TEST/Group1_medLD_CMA.var.sites",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(left_join(GL3_CMA_genos_med[ , c("scaffold", "pos", 
                                         "C_major", "C_minor", 
                                         "M_major", "M_minor", 
                                         "A_major", "A_minor", "rdiff")],
                      ids_counts, by = c("scaffold", "pos")),
            "../local_ancestry/results/ancestry_hmm/TEST/med_LD/Group1_medLD_CMA.counts",
            col.names = F, row.names = F, quote = F, sep = "\t")

# very low LD:
# thin SNPs:
spacing3 = .01 # min cM spacing
keep_snp3 <- rep(F, nrow(GL3))
last_rpos3 <- -100
for (i in 1:nrow(GL3)) {
  this_rpos3 <- GL3[i, "rpos"]
  if (this_rpos3 - last_rpos3 >= spacing3){
    keep_snp3[i] <- T
    last_rpos3 <- this_rpos3 
  } # else do nothing  
}
GL3_CMA_genos_verylow <- left_join(GL3[keep_snp3, c("scaffold", "pos")],
                               CMA_genos, by = c("scaffold", "pos")) %>%
  mutate(rdiff = c(1, format(round(diff(rpos)/100, 12), scientific = F))) # 12 decimals, units Morgans, no scientific number
write.table(GL3_CMA_genos_verylow[ , c("scaffold", "pos", 
                                   "C_major", "C_minor", 
                                   "M_major", "M_minor", 
                                   "A_major", "A_minor", "rdiff")], 
            "../local_ancestry/results/SNPs/TEST/Group1_verylowLD_CMA.genos",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(GL3_CMA_genos_verylow[ , 1:4], "../local_ancestry/results/SNPs/TEST/Group1_verylowLD_CMA.var.sites",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(left_join(GL3_CMA_genos_verylow[ , c("scaffold", "pos", 
                                             "C_major", "C_minor", 
                                             "M_major", "M_minor", 
                                             "A_major", "A_minor", "rdiff")],
                      ids_counts, by = c("scaffold", "pos")),
            "../local_ancestry/results/ancestry_hmm/TEST/verylow_LD/Group1_verylowLD_CMA.counts",
            col.names = F, row.names = F, quote = F, sep = "\t")



# compare results across low med and high LD runs of ancestry_hmm
anc <- lapply(pop_ids, function(id)
  calc_anc_from_post(post = read_post(id = id, dir = dir_input)))
verylow2 <- read_post("AR1103", dir = "../local_ancestry/results/ancestry_hmm/TEST/verylow_LD/")
anc_verylow2 <- cbind(verylow2, calc_anc_from_post(post = verylow2))
low2 <- read_post("AR1103", dir = "../local_ancestry/results/ancestry_hmm/TEST/low_LD/")
anc_low2 <- cbind(low2, calc_anc_from_post(post = low2))
high2 <- read_post("AR1103", dir = "../local_ancestry/results/ancestry_hmm/TEST/high_LD/")
anc_high2 <- cbind(high2, calc_anc_from_post(post = high2))
med2 <- read_post("AR1103", dir = "../local_ancestry/results/ancestry_hmm/TEST/med_LD/")
anc_med2 <- cbind(med2, calc_anc_from_post(post = med2))

anc_low2_match <- inner_join(verylow2[ , c("chrom", "position")],
                                      anc_low2,
                                      by = c("chrom", "position"))
anc_verylow2_match <- left_join(anc_low2_match[ , c("chrom", "position")],
                             anc_verylow2,
                             by = c("chrom", "position"))
anc_high2_match <- left_join(verylow2[ , c("chrom", "position")],
                            anc_high2,
                            by = c("chrom", "position"))
anc_high2_match_low <- left_join(low2[ , c("chrom", "position")],
                             anc_high2,
                             by = c("chrom", "position"))
anc_high2_match_med <- left_join(med2[ , c("chrom", "position")],
                                 anc_high2,
                                 by = c("chrom", "position"))
anc_low2_match_med <- inner_join(med2[ , c("chrom", "position")],
                             anc_low2,
                             by = c("chrom", "position"))
anc_med2_match_low <- inner_join(low2[ , c("chrom", "position")],
                                 anc_med2,
                                 by = c("chrom", "position"))
anc_verylow2_match_med <- inner_join(med2[ , c("chrom", "position")],
                                 anc_verylow2,
                                 by = c("chrom", "position"))
anc_med2_match_verylow <- inner_join(verylow2[ , c("chrom", "position")],
                                 anc_med2,
                                 by = c("chrom", "position"))
cor(anc_verylow2$A, anc_high2_match$A)
plot(anc_verylow2$A, anc_high2_match$A)
cor(anc_verylow2_match$A, anc_low2_match$A)
plot(anc_verylow2_match$A, anc_low2_match$A)
cor(anc_verylow2_match_med$A, anc_med2_match_verylow$A)
plot(anc_verylow2_match_med$A, anc_med2_match_verylow$A)
cor(anc_low2_match_med$A, anc_med2_match_low$A)
plot(anc_low2_match_med$A, anc_med2_match_low$A)
cor(anc_high2_match_med$A, anc_med2$A)
plot(anc_high2_match_med$A, anc_med2$A)

cor(anc_low2$A, anc_high2_match_low$A)
plot(anc_low2$A, anc_high2_match_low$A)
cor(anc_low2$C, anc_high2_match_low$C)
plot(anc_low2$C, anc_high2_match_low$C)
cor(anc_low2$M, anc_high2_match_low$M)
plot(anc_low2$M, anc_high2_match_low$M)
anc_low2[which(anc_low2$A < .2 & anc_high2_match_low$A > .6),]
anc_low2[which(anc_low2$C < .2 & anc_high2_match_low$C > .6),]
expLowLD <- rexp(n = 10000, rate = 60)
