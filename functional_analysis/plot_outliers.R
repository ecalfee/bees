library(dplyr)
library(ggplot2)
library(bedr)
source("../colors.R") # for color palette
load("../local_ancestry/results/A.RData")
load("../local_ancestry/results/C.RData")
load("../local_ancestry/results/M.RData")

# population means:
A_mean = apply(A, 2, mean) %>%
  data.frame(A = ., population = names(.), stringsAsFactors = F)

# this script ranks genes & plots them by evidence for selection
# significance text
sig_text = data.frame(FDR = c(0.01, 0.05, 0.1, NA),
                      stars = c("***", "**", "*", "n.s."))

# load qtls
qtl = read.table("results/Approximate_QTL_positions_HAv3.1.txt", header = T, sep = "\t") %>%
  filter(!is.na(scaffold)) # position found

# load pop/ind metadata for bees
load("../local_ancestry/results/meta.RData")


# load genome metadata
# Compute chromosome sizes
chr_lengths <- read.table("../data/honeybee_genome/chr.list", sep = "\t", stringsAsFactors = F)
colnames(chr_lengths) <- c("scaffold", "length", "chr_group", "chr_lg")
chr_lengths <- chr_lengths %>%
  mutate(chr = as.numeric(substr(chr_lg, 3, 100))) %>%
  mutate(chr_n = chr) %>%
  arrange(chr) %>%
  mutate(chr_start = cumsum(length) - length,
         chr_end = cumsum(length)) %>%
  mutate(chr_midpoint = (chr_start + chr_end)/2)

# #-------------------------------------HAv3.1 Analysis------------------------------------------------------------

# load mean ancestry all snps with ancestry calls genomewide
load("../local_ancestry/results/mean_ancestry_AR_CA.RData") # loads A_AR_CA

# get mean ancestry for all genes
gene_file <- "results/HAv3.1_genes_only.chr.gff3.sorted.bed"
gene_file_columns <- c("scaffold", "start", "end", "source", "gene", "gene_info")


# group outliers into contiguous regions
threshold_extend_region <- 0.1

# start with high A shared outliers
# high shared outliers defined by BOTH being subthreshold 10% FDR in CA and AR
threshold_keep_region = 0.1
high.shared2 <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_high", "FDR_CA_high", "FDR_shared_high")) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high),
         FDR_CA_high = as.numeric(FDR_CA_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_high) & !is.na(FDR_CA_high)) %>%
  filter(FDR_AR_high <= threshold_extend_region & FDR_CA_high <= threshold_extend_region)
#table(is.valid.region(high.shared2, check.chr = F))
#length(bedr.merge.region(high.shared2, check.chr = F))
high.shared.outliers2 <- bedr(
  engine = "bedtools",
  input = list(i = high.shared2),
  method = "merge",
  params = "-d 10000 -c 7,8 -o min", # merge if within 10kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR_AR", "min_FDR_CA", "region")) %>%
  filter(min_FDR_AR <= threshold_keep_region & min_FDR_CA <= threshold_keep_region)
#high.shared.outliers2

# low shared outliers (NONE):
low.shared <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_shared_low")) %>%
  mutate(FDR_shared_low = as.numeric(FDR_shared_low)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_shared_low)) %>%
  filter(FDR_shared_low <= threshold_extend_region)
#dim(low.shared)

# high CA outliers:
high.CA <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_CA_high")) %>%
  mutate(FDR_CA_high = as.numeric(FDR_CA_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_CA_high)) %>%
  filter(FDR_CA_high <= threshold_extend_region)

# high AR outliers:
high.AR <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_high")) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_high)) %>%
  filter(FDR_AR_high <= threshold_extend_region)

# low AR outliers (there are no low CA outliers -- underpowered. Also I found no low shared outliers.)
low.AR <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_low")) %>%
  mutate(FDR_AR_low = as.numeric(FDR_AR_low)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_low)) %>%
  filter(FDR_AR_low <= threshold_extend_region)

low.AR.only <- low.AR # no low CA outliers

# merge into contiguous regions
low.AR.only.outliers <- bedr(
  engine = "bedtools",
  input = list(i = low.AR.only),
  method = "merge",
  params = "-d 10000 -c 7 -o min", # merge if within 10kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "region")) %>%
  filter(min_FDR <= threshold_keep_region)

# don't filter out shared sites, just ID them:
high.CA.intersect <- bedr(
  engine = "bedtools",
  input = list(a = high.CA[ , c("chr", "start", "end")],
               b = high.shared.outliers2[ , c("chr", "start", "end")]),
  method = "intersect",
  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths", # exclude if it's in a 'shared outlier' region
  check.chr = F
) %>%
  left_join(high.CA, ., by = c("chr", "start", "end")) %>%
  mutate(bp_shared_outliers = V7) %>%
  dplyr::select(c(colnames(high.CA), bp_shared_outliers))

# merge into contiguous regions
high.CA.outliers <- bedr(
  engine = "bedtools",
  input = list(i = high.CA.intersect),
  method = "merge",
  params = "-d 10000 -c 7,8 -o min,sum", # merge if within 10kb
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "bp_shared_outliers")) %>%
  mutate(region = rownames(.)) %>%
  mutate(start = as.integer(start),
         end = as.integer(end),
         min_FDR = as.numeric(min_FDR),
         bp_shared_outliers = as.integer(bp_shared_outliers),
         percent_bp_shared = bp_shared_outliers/(end - start)) %>%
  dplyr::select(c("chr", "start", "end", "min_FDR", "region", "bp_shared_outliers", "percent_bp_shared")) %>%
  filter(min_FDR <= threshold_keep_region) # keep only outliers meeting threshold

high.AR.intersect <- bedr(
  engine = "bedtools",
  input = list(a = high.AR[ , c("chr", "start", "end")],
               b = high.shared.outliers2[ , c("chr", "start", "end")]),
  method = "intersect",
  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths", # exclude if it's in a 'shared outlier' region
  check.chr = F
) %>%
  left_join(high.AR, ., by = c("chr", "start", "end")) %>%
  mutate(bp_shared_outliers = V7) %>%
  dplyr::select(c(colnames(high.AR), bp_shared_outliers))

# merge into contiguous regions
high.AR.outliers <- bedr(
  engine = "bedtools",
  input = list(i = high.AR.intersect),
  method = "merge",
  params = "-d 10000 -c 7,8 -o min,sum", # merge if within 10kb
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "bp_shared_outliers")) %>%
  mutate(region = rownames(.)) %>%
  mutate(start = as.integer(start),
         end = as.integer(end),
         min_FDR = as.numeric(min_FDR),
         bp_shared_outliers = as.integer(bp_shared_outliers),
         percent_bp_shared = bp_shared_outliers/(end - start)) %>%
  dplyr::select(c("chr", "start", "end", "min_FDR", "region", "bp_shared_outliers", "percent_bp_shared")) %>%
  filter(min_FDR <= threshold_keep_region) # keep only outliers meeting threshold


# merge into contiguous regions
low.AR.outliers <- low.AR.only.outliers # b/c no shared low outliers

# write files with outliers regions:
outlier_columns <- c("chr", "start", "end", "region", "min_FDR", "bp_shared_outliers", "percent_bp_shared", "min_FDR_AR", "min_FDR_CA")
high.shared.outliers3 <- high.shared.outliers2 %>%
  mutate(min_FDR_AR = as.numeric(min_FDR_AR),
         min_FDR_CA = as.numeric(min_FDR_CA),
         min_FDR = ifelse(min_FDR_AR < min_FDR_CA, min_FDR_CA, min_FDR_AR),
         bp_shared_outliers = end - start,
         percent_bp_shared = 1) %>%
  dplyr::select(outlier_columns)
low.AR.outliers3 <- low.AR.outliers %>%
  mutate(min_FDR = as.numeric(min_FDR), min_FDR_AR = min_FDR, min_FDR_CA = NA,
         bp_shared_outliers = 0, percent_bp_shared = 0) %>%
  dplyr::select(outlier_columns)
high.AR.outliers3 <- high.AR.outliers %>%
  mutate(min_FDR_AR = min_FDR, min_FDR_CA = NA) %>%
  dplyr::select(outlier_columns)
high.CA.outliers3 <- high.CA.outliers %>%
  mutate(min_FDR_CA = min_FDR, min_FDR_AR = NA) %>%
  dplyr::select(outlier_columns)

outlier_sets <- list(high.shared.outliers3, low.AR.outliers3, high.AR.outliers3, high.CA.outliers3) # only 3 categories of outliers at 5% FDR.
outlier_set_names <- c("high_shared2", "low_AR", "high_AR", "high_CA")

for (i in 1:length(outlier_sets)){
  write.table(outlier_sets[[i]],
              paste0("results/outlier_regions/", outlier_set_names[i], ".bed"),
              quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(outlier_sets[[i]],
              paste0("results/outlier_regions/", outlier_set_names[i], ".noHeader.bed"),
              quote = F, col.names = F, row.names = F, sep = "\t")
}

# write one file with all types of outliers
# NOTE: These merge outlier regions if within 10kb
outliers_all <- do.call(bind_rows,
                        lapply(1:length(outlier_sets), function(i)
   return(mutate(outlier_sets[[i]], outlier_type = outlier_set_names[i])))) %>%
  rename(bp_overlap_shared_outlier = bp_shared_outliers) %>%
  mutate(., length = end - start)

outliers_all_genome_sort <- bedr(
  engine = "bedtools",
  input = list(i = outliers_all),
  method = "sort",
  params = "-faidx ../data/honeybee_genome/Amel_HAv3.1.fasta.fai",
  check.chr = F
) %>%
  data.table::setnames(colnames(outliers_all)) %>%
  mutate(region_n = 1:nrow(.))

# NOTE: These merge outlier regions if within 10kb
write.table(outliers_all_genome_sort,
            paste0("results/outlier_regions/all.bed"),
            quote = F, col.names = T, row.names = F, sep = "\t")
write.table(outliers_all_genome_sort,
            paste0("results/outlier_regions/all.noHeader.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")

# add 20kb buffer on either side of the outlier region & put in genome order
outliers_all_buffer <- bedr(
  engine = "bedtools",
  input = list(i = outliers_all_genome_sort),
  method = "slop",
  params = "-b 20000 -g ../data/honeybee_genome/chr.lengths",
  check.chr = F
) %>%
  data.table::setnames(colnames(outliers_all_genome_sort)) %>%
  mutate(region_w_buffer = paste0(chr, ":", start, "-", end))
# write output files:
write.table(outliers_all_buffer,
            paste0("results/outlier_regions/all.plus20kb.bed"),
            quote = F, col.names = T, row.names = F, sep = "\t")
write.table(outliers_all_buffer,
            paste0("results/outlier_regions/all.plus20kb.noHeader.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")


#******************************************************** WHICH GENES FALL IN OUTLIER REGIONS?
# simple case: which genes
genes_high_AR <- bedr(
    engine = "bedtools",
    input = list(a = gene_file,
                 b = A_AR_CA %>%
                   mutate(FDR_AR_high = as.numeric(FDR_AR_high)) %>%
                   filter(!is.na(FDR_AR_high)) %>%
                   dplyr::select(c("scaffold", "start", "end", "FDR_AR_high")) %>%
                   rename(chr = scaffold)),
    method = "map",
    params = "-g ../data/honeybee_genome/chr.lengths -c 4, -o min",
    check.chr = F
  ) %>%
  data.table::setnames(c(gene_file_columns, "FDR_AR_high")) %>%
  filter(FDR_AR_high != ".")
genes_low_AR <- bedr(
  engine = "bedtools",
  input = list(a = gene_file,
               b = A_AR_CA %>%
                 mutate(FDR_AR_low = as.numeric(FDR_AR_low)) %>%
                 filter(!is.na(FDR_AR_low)) %>%
                 dplyr::select(c("scaffold", "start", "end", "FDR_AR_low")) %>%
                 rename(chr = scaffold)),
  method = "map",
  params = "-g ../data/honeybee_genome/chr.lengths -c 4, -o min",
  check.chr = F
) %>%
  data.table::setnames(c(gene_file_columns, "FDR_AR_low")) %>%
  filter(FDR_AR_low != ".")
genes_high_CA <- bedr(
  engine = "bedtools",
  input = list(a = gene_file,
               b = A_AR_CA %>%
                 mutate(FDR_CA_high = as.numeric(FDR_CA_high)) %>%
                 filter(!is.na(FDR_CA_high)) %>%
                 dplyr::select(c("scaffold", "start", "end", "FDR_CA_high")) %>%
                 rename(chr = scaffold)),
  method = "map",
  params = "-g ../data/honeybee_genome/chr.lengths -c 4, -o min",
  check.chr = F
) %>%
  data.table::setnames(c(gene_file_columns, "FDR_CA_high")) %>%
  filter(FDR_CA_high != ".")



genes_high_shared <- inner_join(genes_high_CA, genes_high_AR, by = gene_file_columns)
genes_high_CA_only <- left_join(genes_high_CA, genes_high_AR, by = gene_file_columns) %>%
  filter(is.na(FDR_AR_high))
genes_high_AR_only <- left_join(genes_high_AR, genes_high_CA, by = gene_file_columns) %>%
  filter(is.na(FDR_CA_high))
genes_combined <- bind_rows(mutate(genes_high_shared, outlier_type = "high_shared"),
               mutate(genes_high_CA_only, outlier_type = "high_CA_only"),
               mutate(genes_high_AR_only, outlier_type = "high_AR_only"),
               mutate(genes_low_AR, outlier_type = "low_AR_only")) %>%
  dplyr::arrange(scaffold, start)
write.table(genes_combined, "results/genes_0.1FDR_combined.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
# # note: all subsequent analyses use 'genes_combined" from txt file above


outlier_types <- c("shared_high", "AR_high", "CA_high", # no CA_low outliers
                   "shared_low", "AR_low")

# get approximate cumulative positions for all points with ancestry calls
pretty_label_zone = data.frame(zone = c("CA", "AR"),
                               zone_pretty = c("N. America", "S. America"),
                               stringsAsFactors = F)
A_AR_CA_cumulative <- A_AR_CA %>%
  mutate(chr_n = as.numeric(substr(chr, 6, 100))) %>% # turn Group11 into 11
  arrange(chr_n) %>% # sort by chromosome order
  tidyr::gather(., "zone", "A_ancestry", c("CA", "AR")) %>%
  mutate(FDR = sapply(1:nrow(.), function(i) ifelse(.$zone[i] == "CA",
                                                    min(.$FDR_CA_high[i], .$FDR_CA_low[i], na.rm = T),
                                                    min(.$FDR_AR_high[i], .$FDR_AR_low[i], na.rm = T))),
         FDR = ifelse(FDR == Inf, NA, FDR)) %>%
  mutate(color_by = ifelse(is.na(FDR), ifelse((chr_n %% 2 == 0), # even chromosomes different color
                                              "n.s. - even chr", "n.s. - odd chr"), FDR)) %>%
  left_join(., pretty_label_zone, by = "zone")



#********************************************************************************************************************
# all 3 ancestries together
ACM_AR_CA <- A_AR_CA[ , c("scaffold", "start", "end", "chr", "pos", "cum_pos", "snp_id",
                          "FDR_shared_high", "FDR_AR_high", "FDR_CA_high", "FDR_AR_low")] %>%
  mutate(A_AR = meanA_AR, A_CA = meanA_CA,
         M_AR = meanM_AR, M_CA = meanM_CA,
         C_AR = meanC_AR, C_CA = meanC_CA) %>%
  tidyr::pivot_longer(data = ., cols = c("A_AR", "A_CA", "C_AR", "C_CA", "M_AR", "M_CA"),
               values_to = "ancestry_freq", names_to = "ancestry") %>%
  tidyr::separate(., "ancestry", c("ancestry", "zone")) %>%
  left_join(., pretty_label_zone, by = "zone")

ACM_means <- ACM_AR_CA %>%
  group_by(ancestry, zone, zone_pretty) %>%
  summarise(ancestry_sd = sd(ancestry_freq),
            ancestry_freq = mean(ancestry_freq))

# Ancestry informative markers (AIMs) for all 3 ancestries
get_aim_freq <- function(zone, ancestry, chr_n){
  # sites
  s <-  read.table(paste0("../clines/results/AIMs/", ancestry, "/Group", chr_n, ".ACM.freqs"),
                   header = T)
  # get maf all pops
  mafs <- do.call(cbind,
                  lapply(meta.pop$population[meta.pop$zone == zone],
                         function(pop) # for all sites, get pop freq
                           left_join(s[ , c("scaffold", "pos", "major", "minor")],
                                     read.table(paste0("../clines/results/AIMs/", ancestry, "/Group", chr_n,
                                                       "/", pop, ".mafs.gz"), header = T),
                                     by = c("scaffold"="chromo", "pos"="position", "major", "minor")) %>%
                           dplyr::select(phat)))
  # summarise
  s %>%
    mutate(freq_with_some_NA = apply(mafs, 1, function(x)
      sum(x * meta.pop$n_bees[meta.pop$zone == zone], na.rm = T)/
        sum(ifelse(is.na(x), 0 , 1) * meta.pop$n_bees[meta.pop$zone == zone], na.rm = T)),
      n_NA = apply(mafs, 1, function(x) sum(is.na(x))),
      freq = apply(mafs, 1, function(x)
        sum(x * meta.pop$n_bees[meta.pop$zone == zone])/
          sum(meta.pop$n_bees[meta.pop$zone == zone])),
      zone = zone,
      AIM_ancestry = ancestry,
      chr_n = chr_n)
}
# get N and S American allele freqs at all AIMs
aims0 <- do.call(rbind,
                   lapply(c("N. America", "S. America"),
                          function(z) do.call(rbind,
                                              lapply(ACM, function(a)
                                                do.call(rbind,
                                                        lapply(1:16, function(i)
                                                          get_aim_freq(zone = z,
                                                                       ancestry = a,
                                                                       chr_n = i)))))))
aims0$flip = sapply(1:nrow(aims0), function(i) aims0[i, paste0("freq_", aims0$AIM_ancestry[i])] < 0.5) # flip if freq M for an M aim is low (not high)
aims0$freq_polarized <- ifelse(aims0$flip, 1 - aims0$freq, aims0$freq)
aims0$chr = paste0("Group", aims0$chr_n)
aims <- aims0 %>% # add filter for whether all pops have some coverage/data
  pivot_wider(data = ., id_cols = c("scaffold", "pos", "AIM_ancestry"),
              names_from = "zone", values_from = "freq_polarized") %>%
  mutate(all_pops_have_data = !is.na(`N. America` + `S. America`)) %>%
  dplyr::select(scaffold, pos, AIM_ancestry, all_pops_have_data) %>%
  left_join(aims0, ., by = c("scaffold", "pos", "AIM_ancestry")) %>%
  left_join(., A_AR_CA[ , c("scaffold", "pos", "snp_id")], # identify snps also used for ancestry_hmm
            by = c("scaffold", "pos")) %>%
  mutate(hmm_marker = !is.na(snp_id))
# nrow(aims)/2 # ~38k

AIMS_ACM_AR_CA <- aims %>%
  rename(., ancestry = AIM_ancestry, zone_pretty = zone) %>%
  filter(all_pops_have_data) %>% # filter out any aims without data for all pops in N and S America
  rename(frequency = freq_polarized) %>%
  dplyr::select(-freq) %>%
  mutate(type = "freq_polarized") %>%
  bind_rows(.,
            rename(ACM_AR_CA, frequency = ancestry_freq) %>%
              mutate(type = "ancestry_freq")) %>%
  left_join(., data.frame(type = c("ancestry_freq", "freq_polarized"),
                          type_pretty = factor(c("Local ancestry", "AIMs"),
                                               levels = c("Local ancestry", "AIMs"),
                                               order = T),
                          stringsAsFactors = F), by = "type") %>%
  mutate(pos_Mb = pos/10^6)
# AIMS_ACM_AR_CA %>%
#   group_by(ancestry, type, hmm_marker) %>%
#   summarise(n = n())
# aims %>%
#   filter(all_pops_have_data) %>%
#   group_by(AIM_ancestry, hmm_marker) %>%
#   summarise(n = n())

# zoom in on chr11:
p11_outliers <- AIMS_ACM_AR_CA %>%
  filter(chr == "Group11") %>%
  filter(., pos > 1.3*10^7 & pos < 1.6*10^7) %>%
  filter(., type == "ancestry_freq") %>%
  ggplot(.) +
  geom_hline(data = ACM_means, aes(yintercept = ancestry_freq, color = ancestry),
             linetype = "solid") + # dashed
  geom_point(
    aes(x = pos_Mb, y = frequency,
        color = ancestry, size = type)) +
  xlab("Chromosome 11 (Mb)") +
  ylab("Ancestry frequency") +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "low_AR") %>%
              filter(chr == 11) %>%
              mutate(zone = "SA", zone_pretty = "S. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  scale_color_manual(values = col_ACM, name = "Ancestry") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  scale_size_manual(values = c(ancestry_freq = 0.05, freq_polarized = 0.05), guide = F) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linetype = "blank"))) +
  theme_classic() +
  facet_grid(zone_pretty ~ .)
#p11_outliers
ggsave("plots/ACM_frequency_plot_AR_CA_FDR_chr11_outlier_ACM.png",
       p11_outliers,
       height = 5, width = 7.5, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/ACM_frequency_plot_AR_CA_FDR_chr11_outlier_ACM.png",
       p11_outliers,
       height = 5, width = 7.5, units = "in", dpi = 600, device = "png")


# now with AIMs too:
p11_aims <- AIMS_ACM_AR_CA %>%
  filter(chr == "Group11") %>%
  filter(., pos > 1.3*10^7 & pos < 1.6*10^7) %>%
  ggplot(.) +
  geom_hline(data = ACM_means, aes(yintercept = ancestry_freq, color = ancestry),
             linetype = "solid") + # dashed
  geom_point(
             aes(x = pos_Mb, y = frequency,
                 color = ancestry, size = type,
                 shape = ifelse((hmm_marker | type == "ancestry_freq"),
                                                               "overlap", "indep"))) +
  xlab("Chromosome 11 (Mb)") +
  ylab("Frequency") +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                                chr_lengths, by = "scaffold") %>%
                 filter(outlier_type == "low_AR") %>%
                 filter(chr == 11) %>%
              mutate(zone = "SA", zone_pretty = "S. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +

  scale_color_manual(values = col_ACM, name = "Ancestry") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  scale_shape_manual(values = c(indep = 4, overlap = 1),
                     labels = c("AIM only", "HMM marker"),
                     name = "Marker overlap") +
  scale_size_manual(values = c(ancestry_freq = 0.05, freq_polarized = 0.2), guide = F) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linetype = "blank"))) +
  theme_classic() +
  facet_grid(type_pretty ~ zone_pretty)

ggsave("plots/ACM_frequency_plot_AR_CA_FDR_chr11_outlier_ACM2.png",
       p11_aims,
       height = 5, width = 7.5, units = "in", device = "png")
ggsave("../../bee_manuscript/figures_supp/ACM_frequency_plot_AR_CA_FDR_chr11_outlier2.tiff",
       p11_aims,
       height = 5, width = 7.5, units = "in", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")
ggsave("../../bee_manuscript/figures/ACM_frequency_plot_AR_CA_FDR_chr11_outlier2.png",
       p11_aims,
       height = 5, width = 7.5, units = "in", dpi = 600, device = "png")

# what percentile is the peak for M ancestry in N. America?
ACM_AR_CA %>%
  filter(chr == "Group11") %>%
  filter(pos > 1.3*10^7 & pos < 1.6*10^7) %>%
  filter(ancestry == "M") %>%
  group_by(zone) %>%
  summarise(mean = mean(ancestry_freq, na.rm = T),
            max = max(ancestry_freq, na.rm = T),
            min = min(ancestry_freq, na.rm = T))

# largest peak of M in N. America within selected region in S. America
ACM_AR_CA %>%
  filter(zone == "AR" & ancestry == "M") %>%
  summarise(percent_over_0.402 = sum(ancestry_freq >= 0.402, na.rm = T)/sum(!is.na(ancestry_freq)))
# 2nd largest peak of M in N. America in the selected region in S. America
ACM_AR_CA %>%
  filter(zone == "CA" & ancestry == "M", chr == "Group11", pos == 14366137) %>%
  dplyr::select(scaffold, chr, pos, ancestry, zone, ancestry_freq, zone_pretty)
ACM_AR_CA %>%
  filter(zone == "AR" & ancestry == "M") %>%
  summarise(percent_over_0.402 = sum(ancestry_freq >= 0.381, na.rm = T)/sum(!is.na(ancestry_freq)))

# shared high outliers. where are they?
filter(outliers_all, outlier_type == "high_shared2") %>%
  group_by(chr, outlier_type) %>%
  summarise(min = min(start),
            max = max(end),
            diff = max - min,
            n = n())
filter(outliers_all, outlier_type == c("high_shared2", "high_CA", "high_AR")) %>%
  group_by(chr, outlier_type) %>%
  summarise(min = min(start),
            max = max(end),
            diff = max - min,
            n = n())
filter(outliers_all) %>%
  group_by(outlier_type) %>%
  summarise(n = n())
filter(outliers_all) %>%
  filter(outlier_type != "low_AR") %>%
  group_by(chr) %>%
  summarise(n = n())
filter(outliers_all) %>%
  filter(outlier_type == "low_AR") %>%
  group_by(chr) %>%
  summarise(n = n())


# chr1
p1_outliers <- AIMS_ACM_AR_CA %>%
  filter(chr == "Group1") %>%
  filter(., pos > 1.025*10^7 & pos < 1.225*10^7) %>%
  filter(., type == "ancestry_freq") %>%
  ggplot(.) +
  geom_hline(data = ACM_means, aes(yintercept = ancestry_freq, color = ancestry),
             linetype = "solid") + # dashed
  geom_point(
    aes(x = pos_Mb, y = frequency,
        color = ancestry), size = 0.1) +
  xlab("Chromosome 1 (Mb)") +
  ylab("Ancestry frequency") +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_shared2") %>%
              filter(chr == 1) %>%
              filter(., end > 1.025*10^7 & start < 1.225*10^7) %>%
              mutate(zone = "AR", zone_pretty = "S. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = 0) +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_shared2") %>%
              filter(chr == 1) %>%
              filter(., end > 1.025*10^7 & start < 1.225*10^7) %>%
              mutate(zone = "CA", zone_pretty = "N. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = 0) +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_AR") %>%
              filter(chr == 1) %>%
              filter(., end > 1.025*10^7 & start < 1.225*10^7) %>%
              mutate(zone = "AR", zone_pretty = "S. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = 0.2) +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_CA") %>%
              filter(chr == 1) %>%
              filter(., end > 1.025*10^7 & start < 1.225*10^7) %>%
              mutate(zone = "CA", zone_pretty = "N. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = 0.2) +
  scale_color_manual(values = col_ACM, name = "Ancestry") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  scale_size_manual(values = c(ancestry_freq = 0.05, freq_polarized = 0.05), guide = F) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linetype = "blank"))) +
  theme_classic() +
  facet_grid(zone_pretty ~ .)
#p1_outliers
ggsave("plots/ACM_frequency_plot_AR_CA_FDR_chr1_outlier_ACM.png",
       height = 5, width = 7.5, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/ACM_frequency_plot_AR_CA_FDR_chr1_outlier_ACM.png",
       height = 5, width = 7.5, units = "in", dpi = 600, device = "png")

# now plot with AIMs
p1_aims <- AIMS_ACM_AR_CA %>%
  filter(chr == "Group1") %>%
  filter(., pos > 1.025*10^7 & pos < 1.225*10^7) %>%
  ggplot(.) +
  geom_hline(data = ACM_means, aes(yintercept = ancestry_freq, color = ancestry),
             linetype = "solid") + # dashed
  geom_point(
    aes(x = pos_Mb, y = frequency,
        color = ancestry, size = type,
        shape = ifelse((hmm_marker | type == "ancestry_freq"),
                       "overlap", "indep"))) +
  xlab("Chromosome 1 (Mb)") +
  ylab("Frequency") +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_shared2") %>%
              filter(chr == 1) %>%
              filter(., end > 1.025*10^7 & start < 1.225*10^7) %>%
              mutate(zone = "AR", zone_pretty = "S. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr),
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_shared2") %>%
              filter(chr == 1) %>%
              filter(., end > 1.025*10^7 & start < 1.225*10^7) %>%
              mutate(zone = "CA", zone_pretty = "N. America"),
            aes(xmin = start/10^6, xmax = end/10^6,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  scale_color_manual(values = col_ACM, name = "Ancestry") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  scale_shape_manual(values = c(indep = 4, overlap = 1),
                     labels = c("AIM only", "HMM marker"),
                     name = "Marker overlap") +
  scale_size_manual(values = c(ancestry_freq = 0.05, freq_polarized = 0.2), guide = F) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linetype = "blank"))) +
  theme_classic() +
  facet_grid(type_pretty ~ zone_pretty)

ggsave("plots/ACM_frequency_plot_AR_CA_FDR_chr1_outlier_ACM2.png",
       p1_aims,
       height = 5, width = 7.5, units = "in", device = "png")
ggsave("../../bee_manuscript/figures_supp/ACM_frequency_plot_AR_CA_FDR_chr1_outlier2.tiff",
       p1_aims,
       height = 5, width = 7.5, units = "in", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")
ggsave("../../bee_manuscript/figures/ACM_frequency_plot_AR_CA_FDR_chr1_outlier2.png",
       p1_aims,
       height = 5, width = 7.5, units = "in", dpi = 600, device = "png")

mean_genomewide <- data.frame(A_ancestry = c(mean(meanA_CA), mean(meanA_AR)), zone_pretty = c("N. America", "S. America"))


# simple plot of outliers, whole genome
p_outliers_genomewide <- ggplot() +
  geom_point(data = (A_AR_CA_cumulative %>% # plot grey points first
                       filter(is.na(FDR))),
             aes(x = cum_pos, y = A_ancestry,
                 color = color_by), size = .01) +
  geom_point(data = A_AR_CA_cumulative %>%
               filter(!is.na(FDR)), # then plot sig points on top
             aes(x = cum_pos, y = A_ancestry,
                 color = color_by), size = .01) +
  xlab("Chromosome") +
  ylab("A ancestry") +
  scale_color_manual(name = NULL,
                      values = col_FDR,
                      limits = c("0.01", "0.05", "0.1"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.10 FDR")
  ) +
  # add in means for reference
  geom_hline(data = mean_genomewide, aes(yintercept = A_ancestry), col = "black", linetype = "dashed") +
  scale_x_continuous(label = chr_lengths$chr, breaks = chr_lengths$chr_mid) +
  theme_classic() +
  facet_grid(zone_pretty ~ .) +
  theme(legend.position = "top", legend.margin = margin(t = 0, unit='cm')) +
  guides(colour = guide_legend(override.aes = list(size = 2, shape = 15)))

#p_outliers_genomewide
ggsave("plots/A_frequency_plot_AR_CA_FDR_whole_genome_wide.png",
       p_outliers_genomewide,
       height = 3, width = 7.5, units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures/A_frequency_plot_AR_CA_FDR_whole_genome_wide.png",
      p_outliers_genomewide,
      height = 3, width = 7.5, units = "in", dpi = 600, device = "png")

# put genomewide plot together with 2 outlier regions: p_outliers_genomewide
p_outliers_combined <- arrangeGrob(p_outliers_genomewide + ggtitle("A"),
                                     p1_outliers + ggtitle("B") +
                                     guides(color = guide_legend(override.aes = list(size=2, shape = 15, linetype = 0))),
                                     p11_outliers + ggtitle("C") + theme(legend.position = "none"),
                                   layout_matrix = rbind(c(1,1),
                                                         c(2,3)),
                                   widths = c(5,3))

plot(p_outliers_combined)

ggsave("plots/A_outliers_grob.png",
       plot = p_outliers_combined,
       height = 6, width = 7.5, units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures/A_outliers_grob.png",
       plot = p_outliers_combined,
       height = 6, width = 7.5, units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_main/A_outliers_grob.tiff",
       plot = p_outliers_combined,
       height = 6, width = 7.5, units = "in", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")

######################## QTL ########################
# buffer for visibility only:
buffer_4_visibility = 25000 # add buffer for visibility only

p_small <- ggplot() + # need to fix
   geom_rect(data = qtl %>%
                  left_join(., chr_lengths, by = "scaffold"),
                aes(xmin = start + chr_start - buffer_4_visibility, # add 50kb for visibility only
                    xmax = end + chr_start + buffer_4_visibility,
                    ymin = 0.2, ymax = 0.7, fill = phenotype),
                alpha = 0.75) +
   geom_point(data = A_AR_CA_cumulative %>%
                filter(!is.na(FDR)), # just plot sig points
              aes(x = cum_pos, y = A_ancestry,
                  color = color_by), size = .01) +

   xlab("Chromosome") +
   ylab("Mean African ancestry") +
   scale_color_manual(name = NULL,
                      values = col_FDR,
                      limits = c("0.01", "0.05", "0.1"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.10 FDR")) +
   # add in means for reference
  geom_hline(data = mean_genomewide, aes(yintercept = A_ancestry), col = "black", linetype = "dashed") +
   scale_x_continuous(label = chr_lengths$chr, breaks = chr_lengths$chr_mid) +
   theme_classic() +
   facet_grid(zone_pretty ~ .) +
   theme(legend.position = "bottom", #legend.margin = margin(t = 0, unit='cm'),
         legend.box = "vertical") +
   guides(colour = guide_legend(override.aes = list(size = 2, shape = 15))) +
   labs(fill = element_blank())

p_small + ggtitle("Ancestry outliers and QTLs (+/- 25kb)")
ggsave("plots/A_outliers_plus_QTLs_whole_genome.png",
        plot = p_small + ggtitle("Ancestry outliers and QTLs (+/- 25kb)"),
        height = 4.5, width = 7.5, units = "in", dpi = 600, device = "png")

View(qtl)
qtl_sorted <- qtl %>%
  mutate(start0 = start,
         end0 = end, #orient scaffolds
         start = as.numeric(ifelse(start < end, start0, end0)),
         end = as.numeric(ifelse(start < end, end0, start0)),
         chr = as.character(scaffold)) %>%
  arrange(scaffold, start)
qtl_overlap <- do.call(rbind,
                       lapply(unique(outliers_all$outlier_type), function(t)
  bedr(
  engine = "bedtools",
  input = list(a = filter(outliers_all, outlier_type == t)[ , c("chr", "start", "end", "outlier_type")],
               b = qtl_sorted[ , c("chr", "start", "end", "QTL_name")]),
  method = "map",
  params = "-sorted -o distinct -c 4 -g ../data/honeybee_genome/chr.lengths",
  check.chr = F
) %>% data.table::setnames(., c("chr", "start", "end", "outlier_type", "QTLs")) %>%
  filter(QTLs != ".")))
View(qtl_overlap)


###-------------------------------------------------------###
# read in genes list with beebase IDs:
# first all outlier genes (file made above)
#genes_combined <- read.table("results/genes_0.1FDR_combined.txt", header = T, stringsAsFactors = F)

# then add DAVID gene name results
DAVID_results <- read.table("results/DAVID_results_gene_names.txt",
                            sep = "\t", header = T, stringsAsFactors = F)

remove_id <- function(x, link) stringr::str_replace(x, paste0("(.+)[", link, "]"), "") # turns ID=3253 into just 3253
get_id <- function(x, link) stringr::str_replace(stringr::str_extract(x, paste0("(.+)[", link, "]")), paste0("[", link, "]"), "") # returns "ID"
get_gene_id_info <- function(x, split = ";", link = "="){ # takes in dictionary and returns named vector
  v = strsplit(x, split = split)[[1]] # turns a dictionary string e.g. "ID=89;Gene=20" into a vector, split by ;
  y = remove_id(v, link = link) # returns just the values: 89, 20
  names(y) = get_id(v, link = link) # returns just the IDs: ID, Gene
  return(y)
}

# #gene_id_cols <- c("ID", "Dbxref", "Name", "gbkey", "gene", "gene_biotype")
# get_gene_id_info(x = "BEEBASE:2,GeneID:LOC213", split = ",", link = ":")
genes_combined2 = do.call(bind_rows,
                          lapply(genes_combined$gene_info,
                                 function(x) get_gene_id_info(x, split = ";", link = "="))) %>%
  cbind(dplyr::select(genes_combined, -c(gene_info, gene)), .) %>%
  cbind(.,
        do.call(bind_rows,
                lapply(.$Dbxref, function(x)
                  get_gene_id_info(x, split = ",", link = ":")))) %>%
  left_join(., DAVID_results[ , c("BEEBASE_ID", "Name")] %>%
              rename(DAVID_gene_name = Name),
            by = c("BEEBASE"="BEEBASE_ID")) %>%
  dplyr::arrange(desc(outlier_type), scaffold, start) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high),
         FDR_CA_high = as.numeric(FDR_CA_high),
         FDR_AR_low = as.numeric(FDR_AR_low))

# table(is.na(genes_combined2$BEEBASE))
# DAVID_results$BEEBASE_ID[!(DAVID_results$BEEBASE_ID %in% genes_combined2$BEEBASE)]
# genes_combined2$BEEBASE[!(genes_combined2$BEEBASE %in% DAVID_results$BEEBASE_ID) & !is.na(genes_combined2$BEEBASE)]
 
 
# write genes and DAVID functional information to file
write.table(dplyr::select(genes_combined2, -c(Dbxref, gene)),
            "results/genes_0.1FDR_combined_with_DAVID_gene_names_3.7.20.txt",
            sep = "\t", quote = F,
            col.names = T, row.names = F)

# write out file for supplement
genes_combined2 %>%
  rename(High_A_outlier_North_America_FDR = FDR_CA_high,
         High_A_outlier_South_America_FDR = FDR_AR_high,
         Low_A_outlier_South_America_FDR = FDR_AR_low) %>%
  mutate(Low_A_outlier_North_America_FDR = NA) %>%
  dplyr::select(scaffold, start, end, source, Name, BEEBASE, gene_biotype, DAVID_gene_name, ends_with("FDR")) %>%
  #View()
  write.table(.,
            "../../bee_manuscript/files_supp/Ancestry_outlier_genes.txt",
            sep = "\t", quote = F,
            col.names = T, row.names = F)
 
A_AR_CA_bed <- A_AR_CA %>%
  dplyr::select(-chr) %>%
  rename(chr = scaffold) %>%
  filter(!(is.na(FDR_AR_low) & is.na(FDR_CA_low) & is.na(FDR_AR_high) & is.na(FDR_CA_high))) %>%
  arrange(chr, pos)
mapped_A_AR_CA_high_shared <- bedr(
  engine = "bedtools",
  input = list(a = A_AR_CA_bed,
               b = high.shared.outliers3),
  method = "map",
  params = "-g ../data/honeybee_genome/chr.lengths -c 4,5 -o collapse,min",
  check.chr = F,
  check.sort = F
) %>%
  data.table::setnames(c(colnames(A_AR_CA_bed), "outlier_region_high_shared", "outlier_region_min_FDR_high_shared"))
mapped_A_AR_CA_high_CA <- bedr(
  engine = "bedtools",
  input = list(a = A_AR_CA_bed,
               b = high.CA.outliers3),
  method = "map",
  params = "-g ../data/honeybee_genome/chr.lengths -c 4,5 -o collapse,min",
  check.chr = F,
  check.sort = F
) %>%
  data.table::setnames(c(colnames(A_AR_CA_bed), "outlier_region_high_CA", "outlier_region_min_FDR_high_CA"))
mapped_A_AR_CA_high_AR <- bedr(
  engine = "bedtools",
  input = list(a = A_AR_CA_bed,
               b = high.AR.outliers3),
  method = "map",
  params = "-g ../data/honeybee_genome/chr.lengths -c 4,5 -o collapse,min",
  check.chr = F,
  check.sort = F
) %>%
  data.table::setnames(c(colnames(A_AR_CA_bed), "outlier_region_high_AR", "outlier_region_min_FDR_high_AR"))
mapped_A_AR_CA_low_AR <- bedr(
  engine = "bedtools",
  input = list(a = A_AR_CA_bed,
               b = low.AR.outliers3),
  method = "map",
  params = "-g ../data/honeybee_genome/chr.lengths -c 4,5 -o collapse,min",
  check.chr = F,
  check.sort = F
) %>%
  data.table::setnames(c(colnames(A_AR_CA_bed), "outlier_region_low_AR", "outlier_region_min_FDR_low_AR"))

mapped_A_AR_CA_regions <- cbind(A_AR_CA_bed,
                                mapped_A_AR_CA_high_shared[ , c("outlier_region_high_shared", "outlier_region_min_FDR_high_shared")],
                                mapped_A_AR_CA_high_CA[ , c("outlier_region_high_CA", "outlier_region_min_FDR_high_CA")],
                                mapped_A_AR_CA_high_AR[ , c("outlier_region_high_AR", "outlier_region_min_FDR_high_AR")],
                                mapped_A_AR_CA_low_AR[ , c("outlier_region_low_AR", "outlier_region_min_FDR_low_AR")])
mapped_A_AR_CA_regions[mapped_A_AR_CA_regions == "."] <- NA

# plot single SNP per shared high region with highest combined frequency
p_clines_high_A_shared <- cbind(A, dplyr::select(A_AR_CA, c("scaffold", "pos"))) %>%
  rename(chr = scaffold) %>%
  left_join(mapped_A_AR_CA_regions %>%
              filter(!is.na(outlier_region_high_shared)) %>%
              group_by(outlier_region_high_shared) %>%
              summarise(snp_id = snp_id[which.max(combined)],
                        pos = pos[which.max(combined)],
                        chr = chr[which.max(combined)],
                        FDR_AR_high = FDR_AR_high[which.max(combined)],
                        FDR_CA_high = FDR_CA_high[which.max(combined)]),
            ., by = c("chr", "pos")) %>%
  pivot_longer(cols = colnames(A), names_to = "population", values_to = "A") %>%
  left_join(., meta.pop, by = "population") %>%
  mutate(FDR = ifelse(zone == "S. America", FDR_AR_high, FDR_CA_high)) %>%
  ggplot(data = ., aes(x = abs(lat), y = A)) +
  geom_abline(intercept = 0.5, slope = 0, color = "grey", linetype = "dashed") +
  geom_abline(intercept = c(0, 1), slope = 0, color = "grey") +
  geom_point(aes(col = factor(outlier_region_high_shared), shape = factor(FDR))) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  theme_classic() +
  theme(strip.text.y = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  facet_grid(outlier_region_high_shared~zone, scales = "free_x") + # make shape FDR
  geom_point(data = left_join(meta.pop, A_mean, by = "population") %>%
               mutate(type = "Genomewide mean"), 
             aes(x = abs(lat), y = A), color = "black", pch = 1) +
  labs(shape = "FDR", color = "Shared high A outlier region") +
  ylab("A ancestry frequency") +
  xlab("Degrees latitude from the equator") +
  scale_shape_manual(values = shapes_sig)
#p_clines_high_A_shared
ggsave("plots/outliers_high_shared_top_snp_clines.png", device = "png",
       plot = p_clines_high_A_shared,
       height = 7.5, width = 5.4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/outliers_high_shared_top_snp_clines.png", device = "png",
       plot = p_clines_high_A_shared,
       height = 7.5, width = 5.4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/outliers_high_shared_top_snp_clines.tif", device = "tiff",
       plot = p_clines_high_A_shared,
       height = 7.5, width = 5.4, units = "in", dpi = 600,
       compression = "lzw", type = "cairo")

# low AR only:
p_clines_low_A_AR <- cbind(A, dplyr::select(A_AR_CA, c("scaffold", "pos"))) %>%
  rename(chr = scaffold) %>%
  left_join(mapped_A_AR_CA_regions %>%
              filter(!is.na(outlier_region_low_AR)) %>%
              group_by(outlier_region_low_AR) %>%
              summarise(snp_id = snp_id[which.min(AR)],
                        pos = pos[which.min(AR)],
                        chr = chr[which.min(AR)],
                        FDR_AR_low = FDR_AR_low[which.min(AR)],
                        FDR_CA_low = FDR_CA_low[which.min(AR)]),
            ., by = c("chr", "pos")) %>%
  pivot_longer(cols = colnames(A), names_to = "population", values_to = "A") %>%
  left_join(., meta.pop, by = "population") %>%
  mutate(FDR = ifelse(zone == "S. America", FDR_AR_low, FDR_CA_low)) %>%
  mutate(FDR = ifelse(is.na(FDR), "n.s.", FDR)) %>% # not significant category
  ggplot(data = ., aes(x = abs(lat), y = A)) +
  geom_abline(intercept = 0.5, slope = 0, color = "grey", linetype = "dashed") +
  geom_abline(intercept = c(0, 1), slope = 0, color = "grey") +
  geom_point(aes(col = factor(outlier_region_low_AR), shape = factor(FDR))) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  theme_classic() +
  theme(strip.text.y = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  facet_grid(outlier_region_low_AR~zone, scales = "free_x") + # make shape FDR
  geom_point(data = left_join(meta.pop, A_mean, by = "population") %>%
               mutate(type = "Genomewide mean"), 
             aes(x = abs(lat), y = A), color = "black", pch = 1) +
  labs(shape = "FDR", color = "Low A in S. America outlier region") +
  ylab("A ancestry frequency") +
  xlab("Degrees latitude from the equator") +
  scale_shape_manual(values = shapes_sig)
p_clines_low_A_AR
ggsave("plots/outliers_low_AR_top_snp_clines.png", device = "png",
       plot = p_clines_low_A_AR,
       height = 4, width = 5.4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/outliers_low_AR_top_snp_clines.png", device = "png",
       plot = p_clines_low_A_AR,
       height = 4, width = 5.4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/outliers_low_AR_top_snp_clines.tiff", device = "tiff",
       plot = p_clines_low_A_AR,
       height = 4, width = 5.4, units = "in", dpi = 600)


# both clines together:
top2_snp_clines <- cbind(A, dplyr::select(A_AR_CA, c("scaffold", "pos"))) %>%
  rename(chr = scaffold) %>%
  left_join(mapped_A_AR_CA_regions, ., by = c("chr", "pos")) %>%
  filter(., combined == max(combined) | AR == min(AR)) %>% # top outliers
  pivot_longer(cols = colnames(A), names_to = "population", values_to = "A") %>%
  left_join(., meta.pop, by = "population") %>%
  rename(scaffold = chr) %>%
  left_join(., chr_lengths[ , c("scaffold", "chr", "chr_n")], by = "scaffold") %>%
  mutate(type = paste0("Chr", chr_n, ":", pos)) %>%
  ggplot(data = ., aes(x = abs(lat), y = A, color = type, shape = type)) +
  geom_abline(intercept = 0.5, slope = 0, color = "grey", linetype = "dashed") +
  geom_abline(intercept = c(0, 1), slope = 0, color = "grey") +
  geom_point() +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.position = "top") +
  facet_grid(zone~snp_id, scales = "free_x") + # make shape FDR
  geom_point(data = left_join(meta.pop, A_mean, by = "population") %>%
               mutate(type = "Genomewide mean")) +
  labs(shape = "", fill = "", color = "") +
  ylab("A ancestry frequency") +
  xlab("Degrees latitude from the equator") +
  scale_shape_manual(values = shapes_top2_outliers) +
  scale_color_manual(values = col_top2_outliers)
#top2_snp_clines
ggsave("plots/top_snp_clines.png", device = "png",
       plot = top2_snp_clines,
       height = 5, width = 5.4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/top_snp_clines.png", device = "png",
       plot = top2_snp_clines,
       height = 5, width = 5.4, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/Fig7.tif", device = "tiff",
       plot = top2_snp_clines,
       height = 5, width = 5.4, units = "in", dpi = 600,
       compression = "lzw", type = "cairo")


# print out all outlier ancestry windows
outliers_all2 <- read.table("results/outlier_regions/all.bed", header = T, stringsAsFactors = F) %>%
  left_join(., data.frame(outlier_type = c("low_AR", "high_AR", "high_CA", "high_shared2"),
                          outlier_group = c("Low_A_South_America", "High_A_South_America", "High_A_North_America", "High_A_Shared")),
            by = "outlier_type") %>%
  rename(FDR_South_America = min_FDR_AR,
         FDR_North_America = min_FDR_CA) %>%
  mutate(outlier_type = as.factor(outlier_group, ordered = T,
                                  levels = c("High_A_Shared", "High_A_South_America", "High_A_North_America", "Low_A_South_America"))) %>%
  dplyr::select(chr, start, end, region, outlier_type, FDR_South_America, FDR_North_America, bp_overlap_shared_outlier, percent_bp_shared) %>%
  dplyr::arrange(chr, start, outlier_type)
#View(outliers_all2)
# note: the version uploaded as S3 is the same except some arbitrary ordering between listing a shared outlier or ind outlier first when they start at the same position
write.table(outliers_all2, "../../bee_manuscript/files_supp/Ancestry_outlier_regions.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
