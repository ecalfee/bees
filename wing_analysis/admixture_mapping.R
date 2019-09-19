# admixture mapping for wing traits, starting with length:
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmtest)

# get wing lengths for all bees
wings <- read.table("ALL_Wings_CA_AR.csv",
                      sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(Bee_ID = sapply(Label, function(x) strsplit(x, split = "-")[[1]][3]))

# note: there are some duplicates. one must be mislabelled. may remove the full pair for now.
View(filter(wings, Bee_ID %in% wings$Bee_ID[duplicated(wings$Bee_ID)]))

# load bee data
# Bee IDs
IDs <- read.table(paste0("../bee_samples_listed/byPop/pops_included.IDs"), stringsAsFactors = F,
                  header = F) %>%
  data.table::setnames("Bee_ID")


# get metadata
meta.ind <- read.table("../bee_samples_listed/all.meta", header = T, stringsAsFactors = F, sep = "\t") %>%
  left_join(IDs, ., by = "Bee_ID") %>%
  filter(Bee_ID %in% wings$Bee_ID) %>%
  left_join(., wings, by = "Bee_ID")

# local ancestry calls - A ancestry
# get ancestry frequencies for each population across the genome
dir_results <- "../local_ancestry/Amel4.5_results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot"
#a <- read.table(paste0(dir_results, "/anc/CA1406.A.anc"), stringsAsFactors = F)
ind_A <- lapply(IDs$Bee_ID, function(p) read.table(paste0(dir_results, "/anc/", p, ".A.anc"),
                                            stringsAsFactors = F))
A_wide <- do.call(cbind, ind_A) %>% # combine and sort to match meta.ind
  data.table::setnames(IDs$Bee_ID) #%>%
  #dplyr::select(meta.ind$Bee_ID)
save(A_wide, file = "A_wide.RData")

# genomewide mean ancestry for each individual
meta.ind$alpha <- apply(A_wide, 2, mean)
alpha <- data.frame(Bee_ID = IDs$Bee_ID, 
                    alpha = apply(A_wide, 2, mean),
                    stringsAsFactors = F)

# get sites information:
sites0 <- read.table("../local_ancestry/Amel4.5_results/SNPs/thin1kb_common3/included_scaffolds.pos", stringsAsFactors = F,
                     sep = "\t", header = F)
colnames(sites0) <- c("scaffold", "pos")
sites <- tidyr::separate(sites0, scaffold, c("chr", "scaffold_n"), remove = F) %>%
  mutate(scaffold_n = as.numeric(scaffold_n)) %>%
  mutate(chr_n = as.numeric(substr(chr, 6, 100))) %>%
  mutate(snp_id = paste0("snp", chr_n, ".", scaffold, ".", pos))

A <- cbind(sites, A_wide) %>%
      tidyr::gather(., "Bee_ID", "A", IDs$Bee_ID) %>%
  left_join(., alpha, by = "Bee_ID") %>%
  left_join(., wings, by = "Bee_ID") %>%
  rename(wing_length = Length) %>%
  dplyr::select(., c("chr", "pos", "snp_id", "Bee_ID", "A", "alpha", "wing_length")) %>%
  filter(!is.na(wing_length))

# save as R data object so I can retrieve it more easily in the future
save(A, file = "A.RData")

a_snp <- filter(A, snp_id == "snp1.Group1.1.2714")
b_snp <- filter(A, snp_id == "snp16.Group16.8.1214672")
lma1 <- lm(wing_length ~ alpha + A, data = a_snp)
lma0 <- lm(wing_length ~ alpha, data = a_snp)
lrta.byhand <- -2*(as.numeric(logLik(lma0)) - as.numeric(logLik(lma1)))
lrta <- lrtest(lma0, lma1)
some_snps <- head(unique(A$snp_id))
A_small <- filter(A, snp_id %in% some_snps)
calc_sig <- function(Y, A_snp, A_genome){
  lm1 <- lm(Y ~ A_snp + A_genome)
  lm0 <- lm(Y ~ A_genome)
  lik_ratio_test <- lrtest(lm0, lm1)
  return(c(coef(lm1)[2:3], chisq = lik_ratio_test[[4]][2], pval = lik_ratio_test[[5]][2]))
}
calc_sig(Y = A_small$wing_length[A_small$snp_id == some_snps[1]],
         A_snp = A_small$A[A_small$snp_id == some_snps[1]],
         A_genome = A_small$alpha[A_small$snp_id == some_snps[1]])
A_small

