library(dplyr)
library(ggplot2)
library(bedr)
library(ggrastr)# to rasterize points but not legends in ggplot graphs
source("../colors.R") # for color palette
load("../local_ancestry/results/A.RData")
load("../local_ancestry/results/C.RData")
load("../local_ancestry/results/M.RData")
# this script ranks genes & plots them by evidence for selection
# significance text
sig_text = data.frame(FDR = c(0.01, 0.05, 0.1, NA),
                      stars = c("***", "**", "*", "n.s."))

# load qtls
qtl = read.csv("../data/QTL_markers_Hunt_lab/AMEL_QTLS.csv", sep = ",", stringsAsFactors = F)

# load genome metadata
# Compute chromosome sizes
chr_lengths <- read.table("../data/honeybee_genome/chr.list", sep = "\t", stringsAsFactors = F)
head(chr_lengths)
colnames(chr_lengths) <- c("scaffold", "length", "chr_group", "chr_lg")
chr_lengths <- chr_lengths %>%
  mutate(chr = as.numeric(substr(chr_lg, 3, 100))) %>%
  arrange(chr) %>%
  mutate(chr_start = cumsum(length) - length,
         chr_end = cumsum(length)) %>%
  mutate(chr_midpoint = (chr_start + chr_end)/2)

#-------------------------------------HAv3.1 Analysis------------------------------------------------------------

# load mean ancestry all snps with ancestry calls genomewide
load("../local_ancestry/results/mean_ancestry_AR_CA.RData") # loads A_AR_CA

# load false-discovery-rate cutoffs based on MVN simulations
FDRs <- read.table("../local_ancestry/results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", 
                   stringsAsFactors = F, header = T, sep = "\t")

# get mean ancestry for all genes
gene_file <- "results/HAv3.1_genes_only.chr.gff3.sorted.bed"
gene_file_columns <- c("scaffold", "start", "end", "source", "gene", "gene_info")
local_ancestry_file <- "../local_ancestry/results/mean_ancestry_AR_CA.bed"
# local_ancestry_file columns colnames(A_AR_CA)
genes <- bedr(
  engine = "bedtools", 
  input = list(a = gene_file,
               b = local_ancestry_file), 
  method = "map", 
  params = "-g ../data/honeybee_genome/chr.lengths -c 5,6,13 -o mean",
  check.chr = F
) %>%
  data.table::setnames(c(gene_file_columns, colnames(A_AR_CA[c(5,6,13)])))
nrow(genes)


# group outliers into contiguous regions
threshold_extend_region <- 0.1
threshold_keep_region <- 0.1 # only keep contiguous regions that meet FDR 0.05 cutoff, but extend them out to the 0.1 cutoff

# start with high A shared outliers
high.shared <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_shared_high")) %>%
  mutate(FDR_shared_high = as.numeric(FDR_shared_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_shared_high)) %>%
  filter(FDR_shared_high <= threshold_extend_region)
high.shared2 <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_high", "FDR_CA_high", "FDR_shared_high")) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high),
         FDR_CA_high = as.numeric(FDR_CA_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_high) & !is.na(FDR_CA_high)) %>%
  filter(FDR_AR_high <= threshold_extend_region & FDR_CA_high <= threshold_extend_region)
table(is.valid.region(high.shared, check.chr = F))
length(bedr.merge.region(high.shared, check.chr = F))
high.shared.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = high.shared), 
  method = "merge", 
  params = "-d 10000 -c 7 -o min", # merge if within 10kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "region")) %>%
  filter(min_FDR <= threshold_keep_region)
high.shared.outliers

# high shared outliers defined by BOTH being subthreshold 10% FDR in CA and AR
threshold_keep_region2 = 0.1
high.shared2 <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_high", "FDR_CA_high", "FDR_shared_high")) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high),
         FDR_CA_high = as.numeric(FDR_CA_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_high) & !is.na(FDR_CA_high)) %>%
  filter(FDR_AR_high <= threshold_extend_region & FDR_CA_high <= threshold_extend_region)
table(is.valid.region(high.shared2, check.chr = F))
length(bedr.merge.region(high.shared2, check.chr = F))
high.shared.outliers2 <- bedr(
  engine = "bedtools", 
  input = list(i = high.shared2), 
  method = "merge", 
  params = "-d 10000 -c 7,8 -o min", # merge if within 10kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR_AR", "min_FDR_CA", "region")) %>%
  filter(min_FDR_AR <= threshold_keep_region2 & min_FDR_CA <= threshold_keep_region2)
high.shared.outliers2
table(high.shared.outliers$min_FDR)
table(high.shared.outliers2[, c("min_FDR_CA", "min_FDR_AR")])
table(high.shared.outliers$chr)
table(high.shared.outliers2$chr)

# low shared outliers -- NONE
low.shared <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_shared_low")) %>%
  mutate(FDR_shared_low = as.numeric(FDR_shared_low)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_shared_low)) %>%
  filter(FDR_shared_low <= threshold_extend_region)
table(is.valid.region(low.shared, check.chr = F))
length(bedr.merge.region(low.shared, check.chr = F))
low.shared.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = low.shared), 
  method = "merge", 
  params = "-d 10000 -c 7 -o min", # merge if within 10kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "region"))

# high CA outliers:
high.CA <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_CA_high")) %>%
  mutate(FDR_CA_high = as.numeric(FDR_CA_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_CA_high)) %>%
  filter(FDR_CA_high <= threshold_extend_region)

#high.CA.only <- high.CA # use for now

# high in CA only (not a shared outlier)
# exclude if it's in a 'shared outlier' region
high.CA.only <- bedr(
  engine = "bedtools", 
  input = list(a = high.CA[ , c("chr", "start", "end")],
               b = high.shared.outliers2[ , c("chr", "start", "end")]), 
  method = "intersect", 
  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths", # exclude if it's in a 'shared outlier' region
  check.chr = F
) %>%
  left_join(high.CA, ., by = c("chr", "start", "end")) %>%
  filter(V7 == 0) %>%
  dplyr::select(colnames(high.CA))


# merge into contiguous regions
high.CA.only.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = high.CA.only), 
  method = "merge", 
  params = "-d 10000 -c 7 -o min", # merge if within 1kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "region")) %>%
  filter(min_FDR <= threshold_keep_region)

# high AR outliers:
high.AR <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_high")) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_high)) %>%
  filter(FDR_AR_high <= threshold_extend_region)

#high outliers AR only  -- NONE
high.AR.only <- bedr( # exclude regions that overlap with high shared outliers
  engine = "bedtools", 
  input = list(a = high.AR[ , c("chr", "start", "end")],
               b = high.shared.outliers2[ , c("chr", "start", "end")]), 
  method = "intersect", 
  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths", # exclude if it's in a 'shared outlier' region
  check.chr = F
) %>%
  #View() # I did a visual check.
  left_join(high.AR, ., by = c("chr", "start", "end")) %>%
  filter(V7 == 0) %>%
  dplyr::select(colnames(high.AR))

# merge into contiguous regions
high.AR.only.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = high.AR.only), 
  method = "merge", 
  params = "-d 10000 -c 7 -o min", # merge if within 1kb
  check.chr = F
) %>%
  mutate(region = rownames(.)) %>%
  data.table::setnames(c("chr", "start", "end", "min_FDR", "region")) %>%
  filter(min_FDR <= threshold_keep_region)

# low AR outliers (there are no low CA outliers -- underpowered. Also I found no low shared outliers.)
low.AR <- A_AR_CA %>%
  dplyr::select(c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_AR_low")) %>%
  mutate(FDR_AR_low = as.numeric(FDR_AR_low)) %>%
  rename(chr = scaffold) %>%
  filter(!is.na(FDR_AR_low)) %>%
  filter(FDR_AR_low <= threshold_extend_region)

#low.AR.only <- low.AR # use for now
# NOTE: Group 1.23 has a low outlier that is only partially overlaps with the shared.low but
# but probably shouldn't really be considerd an 'AR' low only
low.AR.only <- low.AR
  #bedr( # exclude regions that overlap with high shared outliers
  #engine = "bedtools", 
  #input = list(a = low.AR[ , c("chr", "start", "end")],
  #             b = low.shared.outliers[ , c("chr", "start", "end")]), 
  #method = "intersect", 
  #params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths", # exclude if it's in a 'shared outlier' region
  #check.chr = F
#) %>%
  #View() # I did a visual check.
  #left_join(., low.AR, by = c("chr", "start", "end")) %>%
  #filter(V7 == 0) %>%
  #dplyr::select(colnames(low.AR))

# merge into contiguous regions
low.AR.only.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = low.AR.only), 
  method = "merge", 
  params = "-d 10000 -c 7 -o min", # merge if within 1kb
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
  #View() # I did a visual check.
  left_join(high.CA, ., by = c("chr", "start", "end")) %>%
  mutate(bp_shared_outliers = V7) %>%
  dplyr::select(c(colnames(high.CA), bp_shared_outliers))

# merge into contiguous regions
high.CA.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = high.CA.intersect), 
  method = "merge", 
  params = "-d 10000 -c 7,8 -o min,sum", # merge if within 1kb
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
               b = high.shared.outliers[ , c("chr", "start", "end")]), 
  method = "intersect", 
  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths", # exclude if it's in a 'shared outlier' region
  check.chr = F
) %>%
  #View() # I did a visual check.
  left_join(high.AR, ., by = c("chr", "start", "end")) %>%
  mutate(bp_shared_outliers = V7) %>%
  dplyr::select(c(colnames(high.AR), bp_shared_outliers))

# merge into contiguous regions
high.AR.outliers <- bedr(
  engine = "bedtools", 
  input = list(i = high.AR.intersect), 
  method = "merge", 
  params = "-d 10000 -c 7,8 -o min,sum", # merge if within 1kb
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

# NA -- no low shared outliers
#low.AR.intersect <- 
#  bedr( 
#  engine = "bedtools", 
#  input = list(a = low.AR[ , c("chr", "start", "end")],
#               b = low.shared.outliers[ , c("chr", "start", "end")]), 
#  method = "intersect", 
#  params = "-sorted -wao -g ../data/honeybee_genome/chr.lengths",
#  check.chr = F
#) %>%
  #View() # I did a visual check.
#  left_join(., low.AR, by = c("chr", "start", "end")) %>%
#  mutate(bp_shared_outliers = V7) %>%
#  dplyr::select(c(colnames(low.AR), bp_shared_outliers))

# merge into contiguous regions
low.AR.outliers <- low.AR.only.outliers # b/c no shared low outliers 
  #bedr(
  #engine = "bedtools", 
  #input = list(i = low.AR.intersect), 
  #method = "merge", 
  #params = "-d 10000 -c 7,8 -o min,sum", # merge if within 1kb
  #check.chr = F
#) %>%
  #data.table::setnames(c("chr", "start", "end", "min_FDR", "bp_shared_outliers")) %>%
  #mutate(region = rownames(.)) %>%
  #dplyr::select(c("chr", "start", "end", "min_FDR", "region", "bp_shared_outliers")) %>%
  #filter(min_FDR <= threshold_keep_region)

# write files with outliers regions:
#outlier_sets <- list(high.shared.outliers, low.shared.outliers, high.AR.outliers, low.AR.outliers, high.CA.outliers)
#outlier_set_names <- c("high_shared", "low_shared", "high_AR", "low_AR", "high_CA")
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
table(genes_combined$outlier_type) # mostly in the low AR category (very large region)
table(genes_combined[ , c("scaffold", "outlier_type")])
write.table(genes_combined, "results/genes_0.1FDR_combined.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
nrow(genes)


gene_outliers <- bedr(
  engine = "bedtools", 
  input = list(a = gene_file,
               b = "results/outlier_regions/all.noHeader.bed"), 
  method = "map", 
  params = "-g ../data/honeybee_genome/chr.lengths -c 5,6,5 -o distinct,distinct,count_distinct",
  check.chr = F
) %>%
  data.table::setnames(c(gene_file_columns, colnames(outliers_all_genome_sort[c(5,6)]), "overlaps_n_regions")) %>%
  filter(overlaps_n_regions > 0) %>%
  left_join(., genes, by = colnames(genes)[1:6])
dim(gene_outliers)
View(gene_outliers)
# Note: several genes are within multiple outliers regions, e.g. high_shared and high_CA or nearby high_shared and high_shared regions
table(gene_outliers$scaffold) # 288 genes, mostly in two large outlier regions
head(gene_outliers)
gene_outliers2 <- gene_outliers %>%
  tidyr::separate(., region, paste0("region", 1:max(.$overlaps_n_regions)), sep = ",") %>%
  tidyr::gather(., "which_region", "region", paste0("region", 1:max(.$overlaps_n_regions))) %>%
  filter(!is.na(region)) %>% # gets rid of duplicates
  dplyr::select(-which_region, -outlier_type) %>%
  left_join(., outliers_all_genome_sort[ , c("region", "min_FDR", "outlier_type", "bp_overlap_shared_outlier", "length")],
            by = "region") %>%
  arrange(scaffold, start)
# write list of genes 2 file:  
write.table(gene_outliers2, "results/genes_within_outlier_regions.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")

# TO DO: make plots of local ancestry showing outliers and genes within outlier regions. And add QTLs.
# TO DO: note any interesting candidate genes/pathways








# plots
# plot just scaffolds with high A shared outliers (later add genes):
for (x in unique(high.shared.outliers$chr)){
  custom_colors <- c("red", "orange", "skyblue", 
                     grey.colors(nrow(filter(high.shared.outliers, 
                                             chr == x)), end = .7),
                     "darkgrey")
  names(custom_colors) <- c("0.01", "0.05", "0.1", 
                            filter(high.shared.outliers, chr == x)$region, 
                            "NA")
  A_AR_CA %>%
    tidyr::gather(., "location", "A_ancestry", c("AR", "CA")) %>%
    filter(scaffold == x) %>%
    ggplot(.) +
    geom_point(aes(x = pos, y = A_ancestry, 
                   color = factor(FDR_shared_high)), size = .3) +
    xlab("bp position on scaffold") +
    ylab("mean African ancestry") +
    scale_colour_manual(values = custom_colors, drop = F) + 
    facet_wrap(~location, nrow = 2, ncol = 1) +
    ggtitle(paste0("A ancestry in both hybrid zones, ", x)) +
    geom_segment(data = filter(high.shared.outliers, chr == x), 
                 aes(x=start, xend=end, y=0.68, yend=0.68, color = region),
                 lwd = 4) +
    ggsave(paste0("plots/shared_high/A_frequency_plot_shared_high_outliers_AR_CA_FDR_", x, ".png"),
           height = 5, width = 10, units = "in", device = "png")
  
}



# high A in CA but not AR:
for (x in unique(high.CA.outliers$chr)){
  custom_colors <- c("red", "orange", "skyblue", 
                     grey.colors(nrow(filter(high.CA.outliers, 
                                             chr == x)), end = .7),
                     "darkgrey")
  names(custom_colors) <- c("0.01", "0.05", "0.1", 
                            filter(high.CA.outliers, chr == x)$region, 
                            "NA")
  ggplot() +
    geom_point(data = filter(A_AR_CA_cumulative, scaffold == x), #%>%
               #mutate(FDR_CA_high = ifelse(zone == "AR", "NA", FDR_CA_high)), 
               aes(x = (start+end)/2, y = A_ancestry, 
                   color = FDR_CA_high), size = .3) +
    
    xlab("bp position on scaffold") +
    ylab("mean African ancestry") +
    scale_colour_manual(values = custom_colors, drop = F) + 
    facet_wrap(~zone, nrow = 2, ncol = 1) +
    ggtitle(paste0("A ancestry in both hybrid zones, ", x)) +
    geom_segment(data = filter(high.CA.outliers, chr == x), 
                 aes(x=start, xend=end, y=0.68, yend=0.68, color = region),
                 lwd = 4) +
    ggsave(paste0("plots/CA_only/A_frequency_plot_CA_high_outliers_AR_CA_FDR_", x, ".png"),
           height = 5, width = 10, units = "in", device = "png")
  
}
# high in AR but not CA:
for (x in unique(high.AR.outliers$chr)){
  custom_colors <- c("red", "orange", "skyblue", 
                     grey.colors(nrow(filter(high.AR.outliers, 
                                             chr == x)), end = .7),
                     "darkgrey")
  names(custom_colors) <- c("0.01", "0.05", "0.1", 
                            filter(high.AR.outliers, chr == x)$region, 
                            "NA")
  ggplot() +
    geom_point(data = filter(A_AR_CA_cumulative, scaffold == x), #%>%
               #mutate(FDR_CA_high = ifelse(zone == "CA", "NA", FDR_AR_high)), 
               aes(x = (start+end)/2, y = A_ancestry, 
                   color = FDR_AR_high), size = .3) +
    
    xlab("bp position on scaffold") +
    ylab("mean African ancestry") +
    scale_colour_manual(values = custom_colors, drop = F) + 
    facet_wrap(~zone, nrow = 2, ncol = 1) +
    ggtitle(paste0("A ancestry in both hybrid zones, ", x)) +
    geom_segment(data = filter(high.AR.outliers, chr == x), 
                 aes(x=start, xend=end, y=0.68, yend=0.68, color = region),
                 lwd = 4) +
    ggsave(paste0("plots/AR_only/A_frequency_plot_AR_high_outliers_AR_CA_FDR_", x, ".png"),
           height = 5, width = 10, units = "in", device = "png")
  
}
# low in AR only
for (x in unique(low.AR.outliers$chr)){
  custom_colors <- c("red", "orange", "skyblue", 
                     grey.colors(nrow(filter(low.AR.outliers, 
                                             chr == x)), end = .7),
                     "darkgrey")
  names(custom_colors) <- c("0.01", "0.05", "0.1", 
                            filter(low.AR.outliers, chr == x)$region, 
                            "NA")
  ggplot() +
    geom_point(data = filter(A_AR_CA_cumulative, scaffold == x), #%>%
               #mutate(FDR_CA_low = ifelse(zone == "CA", "NA", FDR_AR_low)), 
               aes(x = (start+end)/2, y = A_ancestry, 
                   color = FDR_AR_low), size = .3) +
    
    xlab("bp position on scaffold") +
    ylab("mean African ancestry") +
    scale_colour_manual(values = custom_colors, drop = F) + 
    facet_wrap(~zone, nrow = 2, ncol = 1) +
    ggtitle(paste0("A ancestry in both hybrid zones, ", x)) +
    geom_segment(data = filter(low.AR.outliers, chr == x), 
                 aes(x=start, xend=end, y=0.68, yend=0.68, color = region),
                 lwd = 4) +
    ggsave(paste0("plots/AR_only/A_frequency_plot_AR_low_outliers_AR_CA_FDR_", x, ".png"),
           height = 5, width = 10, units = "in", device = "png")
  
}





# get overlap with genes and outliers
genes1 <- read.table("results/amel_OGSv3.2_genes_only.noGroupUn.sorted.bed",
                     stringsAsFactors = F, header = F, sep = "\t") %>%
  data.table::setnames(c("chr", "start", "end", "gene_list", "dff3_type", "gene_ID")) %>%
  mutate(ID = substr(gene_ID, 4, 100)) %>%# take off the ID= in the gene ID
  dplyr::select(chr, start, end, gene_list, ID)

# get mean ancestry for all genes
genes2 <- bedr(
  engine = "bedtools", 
  input = list(a = genes1,
               b = dplyr::select(A_AR_CA, c("scaffold", "start", "end", "snp_id", "AR", "CA")) %>%
                 rename(., chr = scaffold)), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 5,6 -o mean", # merge if within 1kb
  check.chr = F
) %>% 
  data.table::setnames(c("chr", "start", "end", "gene_list", "ID", "AR", "CA")) %>%
  mutate(AR = as.numeric(AR),
         CA = as.numeric(CA)) %>%
  mutate(combined = (AR*21 + CA*17)/(17+21)) # combined average across all pops included
table(is.na(genes2$combined)) # not all genes have ancestry calls

table(is.valid.region(genes2, check.chr = F))

# ID genes that overlap with high shared A ancestry
high.shared.genes0 <- bedr(
  engine = "bedtools", 
  input = list(a = genes2,
               b = high.shared.outliers), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 5,4 -o collapse",
  check.chr = F
) %>%
  data.table::setnames(c(colnames(genes2), "region", "FDR_region")) %>%
  filter(region != ".") %>%
  mutate(FDR_region = as.numeric(FDR_region))

# assign FDR for each gene
high.shared.genes <- bedr(
  engine = "bedtools", 
  input = list(a = high.shared.genes0,
               b = dplyr::select(A_AR_CA, c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_shared_high")) %>%
                 rename(., chr = scaffold) %>%
                 filter(., FDR_shared_high != "NA")), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 7 -o min", # merge if within 1kb
  check.chr = F
) %>% 
  data.table::setnames(c(colnames(genes2), "region", "FDR_region", "FDR_gene")) %>%
  mutate(FDR_gene = as.numeric(FDR_gene)) %>% # some genes individually won't be significant
  left_join(., sig_text, by = c("FDR_gene"="FDR"))


# shared low genes (1 region)
# ID genes that overlap with high shared A ancestry
low.shared.genes0 <- bedr(
  engine = "bedtools", 
  input = list(a = genes2,
               b = low.shared.outliers), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 5,4 -o collapse",
  check.chr = F
) %>%
  data.table::setnames(c(colnames(genes2), "region", "FDR_region")) %>%
  filter(region != ".") %>%
  mutate(FDR_region = as.numeric(FDR_region))

# assign FDR for each gene
low.shared.genes <- bedr(
  engine = "bedtools", 
  input = list(a = low.shared.genes0,
               b = dplyr::select(A_AR_CA, c("scaffold", "start", "end", "snp_id", "AR", "CA", "FDR_shared_low")) %>%
                 rename(., chr = scaffold) %>%
                 filter(., FDR_shared_low != "NA")), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 7 -o min", # merge if within 1kb
  check.chr = F
) %>% 
  data.table::setnames(c(colnames(genes2), "region", "FDR_region", "FDR_gene")) %>%
  mutate(FDR_gene = as.numeric(FDR_gene)) %>% # some genes individually won't be significant
  left_join(., sig_text, by = c("FDR_gene"="FDR"))

# just CA or just AR:
ind.pop.names <- c("FDR_CA_high", "FDR_AR_high", "FDR_AR_low")
ind.pop.outliers <- list(high.CA.outliers, high.AR.outliers, low.AR.outliers)
outliers.ind.genes0 <- lapply(ind.pop.outliers, 
                         function(x) bedr(
  engine = "bedtools", 
  input = list(a = genes2,
               b = x), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 5,4 -o collapse",
  check.chr = F
) %>%
  data.table::setnames(c(colnames(genes2), "region", "FDR_region")) %>%
  filter(region != ".") %>%
  mutate(FDR_region = as.numeric(FDR_region)))

# assign FDR for each gene
outliers.ind.genes <- lapply(1:3, 
                             function(i)
                               bedr(
  engine = "bedtools", 
  input = list(a = outliers.ind.genes0[[i]],
               b = dplyr::select(A_AR_CA, c("scaffold", "start", "end", "snp_id", "AR", "CA", ind.pop.names[i])) %>%
                 rename(., chr = scaffold) %>%
                 filter(., .[ , 7] != "NA")), 
  method = "map", 
  params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 7 -o min", # merge if within 1kb
  check.chr = F
) %>% 
  data.table::setnames(c(colnames(genes2), "region", "FDR_region", "FDR_gene")) %>%
  mutate(FDR_gene = as.numeric(FDR_gene)) %>% # some genes individually won't be significant
  left_join(., sig_text, by = c("FDR_gene"="FDR")))





# print lists of candidate genes
# keep clustered by region but sort priority within scaffold based on combined ancestry
# high
high.shared.genes %>%
  mutate(FDR_gene_stars = paste(FDR_gene, stars)) %>%
  dplyr::select(region, FDR_region, ID, chr, start, end, gene_list, AR, CA, combined, FDR_gene_stars) %>%
  data.table::setnames(c("Outlier Region", "FDR region", "Gene ID", "Scaffold Amel4.5", "Gene Start (bp)", "Gene End (bp)",
                         "Gene List", "Mean A Ancestry Argentina", "Mean A Ancestry California", "Mean A Ancestry Combined Sample", 
                         "FDR gene")) %>%
  write.table(., "results/outlier_genes_list_high_shared_A_ancestry.txt", na = "n.s.",
              sep = "\t", quote = F, col.names = T, row.names = F)

low.shared.genes %>%
  mutate(FDR_gene_stars = paste(FDR_gene, stars)) %>%
  dplyr::select(region, FDR_region, ID, chr, start, end, gene_list, AR, CA, combined, FDR_gene_stars) %>%
  data.table::setnames(c("Outlier Region", "FDR region", "Gene ID", "Scaffold Amel4.5", "Gene Start (bp)", "Gene End (bp)",
                         "Gene List", "Mean A Ancestry Argentina", "Mean A Ancestry California", "Mean A Ancestry Combined Sample", 
                         "FDR gene")) %>%
  write.table(., "results/outlier_genes_list_low_shared_A_ancestry.txt", na = "n.s.",
              sep = "\t", quote = F, col.names = T, row.names = F)

# single population outliers:
high.CA.genes <- outliers.ind.genes[[1]]
high.CA.genes$also_shared <- sapply(high.CA.genes$ID, function(x) x %in% high.shared.genes$ID)
high.CA.genes$outlier_type <- "high CA only"
table(high.CA.genes$also_shared)
high.AR.genes <- outliers.ind.genes[[2]]
high.AR.genes$also_shared <- sapply(high.AR.genes$ID, function(x) x %in% high.shared.genes$ID)
high.AR.genes$outlier_type <- "high AR only"
table(high.AR.genes$also_shared)
low.AR.genes <- outliers.ind.genes[[3]]
low.AR.genes$also_shared <- sapply(low.AR.genes$ID, function(x) x %in% low.shared.genes$ID)
low.AR.genes$outlier_type <- "low AR only"
table(low.AR.genes$also_shared)
# write file with all individual zone outlier genes:
bind_rows(high.CA.genes, high.AR.genes, low.AR.genes) %>%
  mutate(FDR_gene_stars = paste(FDR_gene, stars)) %>%
  dplyr::select(outlier_type, region, FDR_region, ID, chr, start, end, gene_list, AR, CA, combined, FDR_gene_stars) %>%
  data.table::setnames(c("Outlier Type", "Outlier Region", "FDR region", "Gene ID", "Scaffold Amel4.5", "Gene Start (bp)", "Gene End (bp)",
                         "Gene List", "Mean A Ancestry Argentina", "Mean A Ancestry California", "Mean A Ancestry Combined Sample", 
                         "FDR gene")) %>%
  write.table(., "results/outlier_genes_list_one_zone_only_A_ancestry.txt", na = "n.s.",
              sep = "\t", quote = F, col.names = T, row.names = F)

high.shared.genes$outlier_type <- "high shared"
low.shared.genes$outlier_type <- "low shared"
bind_rows(high.shared.genes, low.shared.genes, high.CA.genes, high.AR.genes, low.AR.genes) %>%
  group_by(region, outlier_type) %>%
  summarise(n = n()) %>%
  arrange(outlier_type) %>%
  write.table(., "results/outlier_genes_count_by_region_A_ancestry.txt", na = "n.s.",
              sep = "\t", quote = F, col.names = T, row.names = F)

# some summaries for executive report:
# number of genes per type of outlier
bind_rows(high.shared.genes, low.shared.genes, high.CA.genes, high.AR.genes, low.AR.genes) %>%
  group_by(outlier_type) %>%
  summarise(n = n()) %>%
  dplyr::select(n) %>%
  sum()
# length of outlier regions:
summary(with(high.shared.outliers, end - start))
sum(high.shared.outliers$end - high.shared.outliers$start)/max(chr_lengths$chromosome_end)
nrow(low.shared.outliers)
nrow(high.CA.outliers)
nrow(high.AR.outliers)
nrow(low.AR.outliers)
high.shared.outliers$length = high.shared.outliers$end - high.shared.outliers$start
high.shared.outliers
summary(high.shared.outliers$length)


# get 1 top outlier SNP per outlier region (to draw ancestry clines):
get_SNPs_in_outlier_regions <- lapply(outlier_sets, function(x)
  bedr(
    engine = "bedtools", 
    input = list(a = A_AR_CA,
                 b = x), 
    method = "map", 
    params = "-g ../data/honeybee_genome/ordered_scaffolds.lengths -c 4,5 -o collapse",
    check.chr = F
  ) %>%
    data.table::setnames(c(colnames(A_AR_CA), "FDR_region", "region")) %>%
    filter(region != ".") %>%
    mutate(FDR_region = as.numeric(FDR_region)) %>%
    dplyr::select("snp_id", "FDR_region", "region") %>% # otherwise converts everything to strings
    left_join(., A_AR_CA, by = "snp_id") %>%
    mutate(combined = (AR*nrow(AR_pops_included) + CA*nrow(CA_pops_included))/(nrow(AR_pops_included) + nrow(CA_pops_included)))
  )
top_SNP_outliers <- vector("list", 5)
top_SNP_outliers[[1]] <- get_SNPs_in_outlier_regions[[1]] %>%
  group_by(region, scaffold, chr) %>%
  summarize(max_combined = max(combined),
            FDR_region = min(FDR_region),
            top_snp = snp_id[which(combined == max(combined))]) %>%
  arrange(chr, scaffold) %>%
  mutate(outlier_type = outlier_set_names[1])
top_SNP_outliers[[2]] <- get_SNPs_in_outlier_regions[[2]] %>%
  group_by(region, scaffold, chr) %>%
  summarize(min_combined = min(combined),
            FDR_region = min(FDR_region),
            top_snp = snp_id[which(combined == min(combined))]) %>%
  arrange(chr, scaffold) %>%
  mutate(outlier_type = outlier_set_names[2])
top_SNP_outliers[[3]] <- get_SNPs_in_outlier_regions[[3]] %>%
  group_by(region, scaffold, chr) %>%
  summarize(max_AR = max(AR),
            FDR_region = min(FDR_region),
            top_snp = snp_id[which(AR == max(AR))]) %>%
  arrange(chr, scaffold) %>%
  mutate(outlier_type = outlier_set_names[3])
top_SNP_outliers[[4]] <- get_SNPs_in_outlier_regions[[4]] %>%
  group_by(region, scaffold, chr) %>%
  summarize(min_AR = min(AR),
            FDR_region = min(FDR_region),
            top_snp = snp_id[which(AR == min(AR))]) %>%
  arrange(chr, scaffold) %>%
  mutate(outlier_type = outlier_set_names[4])
top_SNP_outliers[[5]] <- get_SNPs_in_outlier_regions[[5]] %>%
  group_by(region, scaffold, chr) %>%
  summarize(max_CA = max(CA),
            FDR_region = min(FDR_region),
            top_snp = snp_id[which(CA == max(CA))]) %>%
  arrange(chr, scaffold) %>%
  mutate(outlier_type = outlier_set_names[5])
top_SNP_outliers_all <- do.call(rbind, top_SNP_outliers)








# plot top outlier regions

outlier_types <- c("shared_high", "AR_high", "CA_high", # no CA_low outliers
                   "shared_low", "AR_low")

# get approximate cumulative positions for all points with ancestry calls
pretty_label_zone = data.frame(zone = c("CA", "AR"),
                               zone_pretty = c("N. America", "S. America"),
                               #zone_pretty = c("California", "Argentina"),
                               stringsAsFactors = F)
# buffer for visibility only:
buffer_4_visibility = 25000*4 # add 50kb for visibility only
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



# zoom in on chr1:
A_AR_CA_cumulative %>%
  filter(., chr_n == 1) %>%
  filter(., pos >= 0.75*10^7 & pos <= 1.25*10^7) %>%
  mutate(FDR = ifelse(zone == "CA", FDR_CA_high, FDR_AR_high)) %>% # just ind. zone FDR's
  mutate(color_by = ifelse(is.na(FDR), ifelse((chr_n %% 2 == 0), # even chromosomes different color
                                              "n.s. - even chr", "n.s. - odd chr"), FDR)) %>%
  #apply(., 2, function(col) sum(is.na(col)))
  #table(.$color_by)
  ggplot(.) +
  geom_point(data = . %>%
               filter(is.na(FDR)), # plot grey points first, rasterized
             aes(x = cum_pos, y = A_ancestry, 
                 color = color_by), size = .05)+
  geom_point(data = . %>%
               filter(!is.na(FDR)), # then plot sig points on top 
             aes(x = cum_pos, y = A_ancestry, 
                 color = color_by), size = .05) +
  xlab("Position (bp)") +
  ylab("mean African ancestry") +
  scale_colour_manual(name = NULL,
                      values = c("0.01"="red", "0.05"="orange", "0.1"="skyblue", 
                                 "n.s. - even chr"="darkgrey", 
                                 "n.s. - odd chr"="grey"),
                      limits = c("0.01", "0.05", "0.1"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.10 FDR")
  ) + 
  geom_segment(data = left_join(rename(outliers_all, scaffold = chr), 
                                #geom_segment(data = left_join(rename(high.shared.outliers, scaffold = chr) %>%
                                #mutate(outlier_type = "high_shared"),
                                chr_lengths, by = "scaffold") %>%
                 filter(chr_n == 1) %>%
                 filter(outlier_type == "high_shared"),
               aes(x = start, xend = end, y = 0.7, yend = 0.7, color = min_FDR),
               lwd = 4) +
  #scale_x_discrete(limits=c("2", "0.5", "1")) +
  #scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  theme_classic() +
  #theme(legend.position = "none") +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1) + 
  ggtitle("5Mb on chr 1 containing shared outlier cluster")
ggsave("plots/A_frequency_plot_AR_CA_FDR_chr1_5mb_outlier_cluster.png",
       height = 5, width = 10, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/A_frequency_plot_AR_CA_FDR_chr1_5mb_outlier_cluster.pdf",
       height = 5, width = 10, units = "in", device = "pdf")

# zoom in more on chr 1
# zoom in on chr1:
A_AR_CA_cumulative %>%
  filter(., chr_n == 1) %>%
  #filter(., pos >= 0.75*10^7 & pos <= 1.25*10^7) %>%
  filter(., pos >= 1*10^7 & pos <= 1.25*10^7) %>%
  mutate(FDR = ifelse(zone == "CA", FDR_CA_high, FDR_AR_high)) %>% # just ind. zone FDR's
  mutate(color_by = ifelse(is.na(FDR), ifelse((chr_n %% 2 == 0), # even chromosomes different color
                                              "n.s. - even chr", "n.s. - odd chr"), FDR)) %>%
  #apply(., 2, function(col) sum(is.na(col)))
  #table(.$color_by)
  ggplot(.) +
  geom_point(data = . %>%
               filter(is.na(FDR)), # plot grey points first, rasterized
             aes(x = cum_pos, y = A_ancestry, 
                 color = color_by), size = .05)+
  geom_point(data = . %>%
               filter(!is.na(FDR)), # then plot sig points on top 
             aes(x = cum_pos, y = A_ancestry, 
                 color = color_by), size = .05) +
  xlab("Position (bp)") +
  ylab("mean African ancestry") +
  scale_colour_manual(name = NULL,
                      values = c("0.01"="red", "0.05"="orange", "0.1"="skyblue", 
                                 "n.s. - even chr"="darkgrey", 
                                 "n.s. - odd chr"="grey"),
                      limits = c("0.01", "0.05", "0.1"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.10 FDR")
  ) + 
  geom_segment(data = left_join(rename(outliers_all, scaffold = chr), 
                                #geom_segment(data = left_join(rename(high.shared.outliers, scaffold = chr) %>%
                                #mutate(outlier_type = "high_shared"),
                                chr_lengths, by = "scaffold") %>%
                 filter(chr_n == 1) %>%
                 filter(outlier_type == "high_shared") %>%
                 filter(start > 1*10^7 & end < 1.25*10^7),
               aes(x = start, xend = end, y = 0.7, yend = 0.7, color = min_FDR),
               lwd = 4) +
  #scale_x_discrete(limits=c("2", "0.5", "1")) +
  #scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  theme_classic() +
  #theme(legend.position = "none") +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1) + 
  ggtitle("2.5Mb around peak of shared outlier cluster on chr1")
ggsave("plots/A_frequency_plot_AR_CA_FDR_chr1_5mb_outlier_cluster_ZOOM_IN.png",
       height = 5, width = 10, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/A_frequency_plot_AR_CA_FDR_chr1_5mb_outlier_cluster_ZOOM_IN.pdf",
       height = 5, width = 10, units = "in", device = "pdf")


# zoom in on chr11:
A_AR_CA_cumulative %>%
  filter(chr_n == 11) %>%
  filter(., pos > 1.25*10^7) %>%
  ggplot(.) + # raster looks pretty terrible -- I could plot instead every 5th or 10th non-sig point or s.t.
  geom_point(data = . %>%
               filter(is.na(FDR)), # plot grey points first
             aes(x = pos, y = A_ancestry, 
                 color = color_by), size = .05) +
  geom_point(data = . %>%
               filter(!is.na(FDR)), # then plot sig points on top 
             aes(x = pos, y = A_ancestry, 
                 color = color_by), size = .05) +
  xlab("Position (bp)") +
  ylab("mean African ancestry") +
  scale_colour_manual(name = NULL,
                      values = c("0.01"="red", "0.05"="orange", "0.10"="skyblue", 
                                 "n.s. - even chr"="darkgrey", 
                                 "n.s. - odd chr"="grey"),
                      limits = c("0.01", "0.05", "0.10"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.10 FDR")
  ) + 
  #geom_segment(data = left_join(rename(outliers_all, scaffold = chr), 
  #                              chr_lengths, by = "scaffold") %>%
  #               filter(outlier_type == "low_AR") %>%
  #               filter(chr == "Group11"),
  #             aes(x = start, xend = end,
  #                 y = 0.7, yend = 0.7, color = min_FDR),
  #             lwd = 4) +
  #scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  #theme(legend.position = "none") +
  theme_classic() +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1) + 
  ggtitle("Low A outlier on chr 11")
ggsave("plots/A_frequency_plot_AR_CA_FDR_chr11_outlier.png",
       height = 5, width = 10, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/A_frequency_plot_AR_CA_FDR_chr11_outlier.pdf",
       height = 5, width = 10, units = "in", device = "pdf")


#********************************************************************************************************************
# all 3 ancestries together
ACM_AR_CA <- A_AR_CA[ , c("scaffold", "start", "end", "chr", "pos", "cum_pos", "snp_id", 
                          "FDR_shared_high", "FDR_AR_high", "FDR_CA_high", "FDR_AR_low")] %>%
  mutate(A_AR = meanA_AR, A_CA = meanA_CA,
         M_AR = meanM_AR, M_CA = meanM_CA,
         C_AR = meanC_AR, C_CA = meanC_CA) %>%
  pivot_longer(data = ., cols = c("A_AR", "A_CA", "C_AR", "C_CA", "M_AR", "M_CA"), 
               values_to = "ancestry_freq", names_to = "ancestry") %>%
  tidyr::separate(., "ancestry", c("ancestry", "zone")) %>%
  left_join(., pretty_label_zone, by = "zone")
ACM_means <- ACM_AR_CA %>%
  group_by(ancestry, zone, zone_pretty) %>%
  summarise(ancestry_sd = sd(ancestry_freq),
            ancestry_freq = mean(ancestry_freq))

# zoom in on chr11:
col_ACM_chr11_outliers = c(col_ACM, overlap_aims_low_M_outliers_11$color)
names(col_ACM_chr11_outliers) = c(names(col_ACM), overlap_aims_low_M_outliers_11$AIM)
ACM_AR_CA %>%
  filter(chr == "Group11") %>%
  filter(., pos > 1.3*10^7 & pos < 1.6*10^7) %>%
  ggplot(.) +
  geom_point(# then plot sig points on top 
             aes(x = pos, y = ancestry_freq, 
                 color = ancestry), size = .05) +
  xlab("Position (bp)") +
  ylab("Mean ancestry frequency") +
  geom_vline(data = overlap_aims_low_M_outliers_11, 
             aes(xintercept = pos), color = "black") +
  #overlap_aims_low_M_outliers_11
  #geom_segment(data = left_join(rename(outliers_all, scaffold = chr), 
  #                              chr_lengths, by = "scaffold") %>%
  #               filter(outlier_type == "low_AR") %>%
  #               filter(chr == 11),
  #             aes(xmin = start, xend = end,
  #                 y = 0.7, yend = 0.7),
  #             lwd = 4) +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr), 
                                chr_lengths, by = "scaffold") %>%
                 filter(outlier_type == "low_AR") %>%
                 filter(chr == 11) %>% 
              mutate(zone = "SA", zone_pretty = "S. America"),
               aes(xmin = start, xmax = end,
                   ymin = -Inf, ymax = Inf),
            alpha = .2) +
  geom_hline(data = ACM_means, aes(yintercept = ancestry_freq, color = ancestry),
             linetype = "solid") + # dashed

  #geom_rect(data = ACM_means,
  #             aes(xmin = 1.3*10^7, xmax = 1.6*10^7,
  #                 ymin = ancestry_freq - 2*ancestry_sd, ymax = ancestry_freq + 2*ancestry_sd,
  #                 fill = ancestry),
  #          alpha = .2) +
  scale_color_manual(values = col_ACM, name = "Ancestry") +
  scale_fill_manual(values = col_ACM, name = "Ancestry") +
  #scale_color_manual(values = col_ACM_chr11_outliers, name = "Ancestry") +
  #scale_fill_manual(values = col_ACM_chr11_outliers, name = "Ancestry") +
  #scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  #theme(legend.position = "none") +
  theme_classic() +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1)
ggsave("plots/ACM_frequency_plot_AR_CA_FDR_chr11_outlier_ACM2.png",
       height = 4, width = 6, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/ACM_frequency_plot_AR_CA_FDR_chr11_outlier2.tiff",
       height = 4, width = 6, units = "in", dpi = 600, device = "tiff")

# chr1
A_AR_CA[which.max(meanA),]
top_pos_shared_high = A_AR_CA[which.max(meanA),]$pos
mean(meanA)
col_ACM_shared_outliers = c(col_ACM, overlap_aims_shared_outliers$color)
names(col_ACM_shared_outliers) = c(names(col_ACM), overlap_aims_shared_outliers$AIM)
left_join(overlap_aims_shared_outliers, A_AR_CA, by = c("chr"="scaffold", "pos"))
sapply(overlap_aims_shared_outliers$pos, function(x) x %in% A_AR_CA[A_AR_CA$chr == "Group1", "pos"])
# none of these markers are in the original ancestry calling set.
ACM_AR_CA %>%
  filter(chr == "Group1") %>%
  filter(., pos > 1*10^7 & pos < 1.25*10^7) %>%
  ggplot(.) +
  geom_vline(data = overlap_aims_shared_outliers, 
             aes(xintercept = pos, color = AIM)) +
  geom_point(# then plot sig points on top 
    aes(x = pos, y = ancestry_freq, 
        color = ancestry), size = .05) +
  xlab("Position (bp)") +
  ylab("Mean ancestry frequency") +
  #geom_segment(data = left_join(rename(outliers_all, scaffold = chr), 
  #                              chr_lengths, by = "scaffold") %>%
  #               filter(outlier_type == "high_shared2") %>%
  #               filter(chr == 1),
  #             aes(xmin = start, xend = end,
  #                 y = 0.7, yend = 0.7),
  #             lwd = 4) +
  geom_rect(data = left_join(rename(outliers_all, scaffold = chr), 
                             chr_lengths, by = "scaffold") %>%
              filter(outlier_type == "high_shared2") %>%
              filter(chr == 1),
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf),
            alpha = .2) +
  #geom_rect(data = left_join(rename(outliers_all, scaffold = chr), 
  #                           chr_lengths, by = "scaffold") %>%
  #            filter(outlier_type == "high_AR") %>%
  #            filter(chr == 1 & start > 1*10^7 & end < 1.25*10^7) %>%
  #            mutate(zone_pretty = "S. America"),
  #          aes(xmin = start, xmax = end,
  #              ymin = -Inf, ymax = Inf),
  #          alpha = .2) +
  #geom_rect(data = left_join(rename(outliers_all, scaffold = chr), 
  #                           chr_lengths, by = "scaffold") %>%
  #            filter(outlier_type == "high_CA") %>%
  #            filter(chr == 1 & start > 1*10^7 & end < 1.25*10^7) %>%
  #            mutate(zone_pretty = "N. America"),
  #          aes(xmin = start, xmax = end,
  #              ymin = -Inf, ymax = Inf),
  #          alpha = .2) +
  geom_hline(data = ACM_means, aes(yintercept = ancestry_freq, color = ancestry),
             linetype = "solid") +
  #geom_rect(data = ACM_means,
  #             aes(xmin = 1.3*10^7, xmax = 1.6*10^7,
  #                 ymin = ancestry_freq - 2*ancestry_sd, ymax = ancestry_freq + 2*ancestry_sd,
  #                 fill = ancestry),
  #          alpha = .2) +
  scale_color_manual(values = col_ACM_shared_outliers, name = "Ancestry") +
  scale_fill_manual(values = col_ACM_shared_outliers, name = "Ancestry") +
  #scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  #theme(legend.position = "none") +
  theme_classic() +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1)
ggsave("plots/ACM_frequency_plot_AR_CA_FDR_chr1_outlier_ACM.png",
       height = 4, width = 6, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/ACM_frequency_plot_AR_CA_FDR_chr1_outlier.tiff",
       height = 4, width = 6, units = "in", dpi = 600, device = "tiff")







# plot with QTLs for varroa on top of 
ggplot() +
  geom_point(data = A_AR_CA_cumulative, aes(x = pos_cum, y = A_ancestry, 
                                            color = color_by), size = .05) +
  xlab("bp position on chromosome") +
  ylab("mean African ancestry") +
  scale_colour_manual(name = NULL,
                      values = c("0.01"="red", "0.05"="orange", "0.1"="skyblue", 
                                 "n.s. - even chr"="darkgrey", 
                                 "n.s. - odd chr"="grey",
                                 "Self-Grooming QTL"="limegreen", 
                                 "Varro-Sensitive Hygiene GWAS SNP"="darkgreen", 
                                 "Hygiene QTL"="olivedrab", 
                                 "Varroa-Sensitive Hygiene putative QTL"="lightseagreen", 
                                 "Varroa-Sensitive Hygiene QTL"="lightgreen"),
                      limits = c("0.01", "0.05", "0.1", 
                                 "n.s. - even chr", 
                                 "n.s. - odd chr",
                                 "Self-Grooming QTL", 
                                 "Varro-Sensitive Hygiene GWAS SNP", 
                                 "Hygiene QTL", 
                                 "Varroa-Sensitive Hygiene putative QTL", 
                                 "Varroa-Sensitive Hygiene QTL"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.1 FDR", 
                                 "n.s. - even chr", 
                                 "n.s. - odd chr",
                                 "Self-Grooming QTL", 
                                 "Varro-Sensitive Hygiene GWAS SNP", 
                                 "Hygiene QTL", 
                                 "Varroa-Sensitive Hygiene putative QTL", 
                                 "Varroa-Sensitive Hygiene QTL")
  ) + 
  scale_x_discrete(limits=c("2", "0.5", "1")) +
  scale_x_continuous(label = chr_lengths$chr, breaks = chr_lengths$chromosome_midpoint) +
  #theme(legend.position = "none") +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1) + 
  geom_segment(data = mutate(varroa_QTL_cumulative, # add padding to GWAS SNPs just to visualize (o/w too small to see)
                             cum_start = ifelse(Outlier_type == "GWAS_SNPs", cum_start - 20000, cum_start),
                             cum_end = ifelse(Outlier_type == "GWAS_SNPs", cum_end + 20000, cum_end)),
               aes(x=cum_start, xend=cum_end, y=0.68, yend=0.68, color = name),
               lwd = 4) +
  theme(legend.position = "bottom") +
  ggtitle("Ancestry outliers across whole genome and co-localization with Varroa defense loci")
ggsave("plots/A_frequency_plot_AR_CA_FDR_whole_genome_wide_VarroaQTL.png",
       height = 5, width = 9, units = "in", device = "png")
# note: the QTL with overlap on chr1 is a putative QTL for removal of varroa-infested brood
# based on an unpublished study in Apis mellifera carnica
# Spotter 2012 "Denser spacing was chosen for nine genomic regions because a preliminary study 
# based on 245 microsatellite loci (M. Brink, M. Solignac, K. Bienefeld, unpublished data) 
# identified QTL for the trait ‘removal of Varroa‐infested brood’ in these regions."

# zoom in just on chromosome 1:
ggplot() +
  geom_point(data = filter(A_AR_CA_cumulative, chr == 1), aes(x = pos_cum, y = A_ancestry, 
                                                              color = color_by), size = .1) +
  xlab("bp position on chromosome") +
  ylab("mean African ancestry") +
  #theme(legend.position = "none") +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1) + 
  geom_segment(data = filter(varroa_QTL_cumulative, chr == 1), 
               aes(x=cum_start, xend=cum_end, y=0.68, yend=0.68, color = name),
               lwd = 4) +
  scale_colour_manual(name = NULL,
                      values = c("0.01"="red", "0.05"="orange", "0.1"="skyblue", 
                                 "n.s. - even chr"="darkgrey", 
                                 "n.s. - odd chr"="grey",
                                 "Varroa-Sensitive Hygiene putative QTL"="lightseagreen", 
                                 "Varroa-Sensitive Hygiene QTL"="lightgreen"),
                      limits = c("0.01", "0.05", "0.1", 
                                 "n.s. - even chr", 
                                 "n.s. - odd chr",
                                 "Varroa-Sensitive Hygiene putative QTL", 
                                 "Varroa-Sensitive Hygiene QTL"),
                      labels = c("0.01 FDR", "0.05 FDR", "0.1 FDR", 
                                 "n.s. - even chr", 
                                 "n.s. - odd chr",
                                 "Varroa-Sensitive Hygiene putative QTL", 
                                 "Varroa-Sensitive Hygiene QTL")
  ) +
  theme(legend.position = "bottom") +
  ggtitle("Ancestry outliers across chromosome 1 and co-localization with Varroa defense loci")
ggsave("plots/A_frequency_plot_AR_CA_FDR_chr1_wide_VarroaQTL.png",
       height = 5, width = 9, units = "in", device = "png")


# how much of the genome do different QTL types cover?
varroa_QTL_cumulative %>%
  group_by(Outlier_type, Reference) %>% 
  summarise(sum = sum(QTL_length))










#################################### OLD ###################################
# load mean ancestry for genes
genes0 <- read.table("results/genes_mean_AR_CA_ancestry.bed", stringsAsFactors = F, 
                     na.strings = c("."), sep = "\t")
colnames(genes0) <- c("scaffold", "start", "end", "gene_list", "gff3_type", "gene_ID", "AR", "CA")

# load FDR for genes -- which genes are outliers?

outlier_genes <- do.call(cbind,
                         lapply(outlier_types, function(x) 
                           read.table(paste0("results/outliers_", x, "_genes_mean_AR_CA_ancestry.bed"),
                                      stringsAsFactors = F, na.strings = c("."))$V9))
colnames(outlier_genes) <- outlier_types
genes <- cbind(genes0, outlier_genes) %>%
  mutate(ID = substr(gene_ID, 4, 100)) %>% # take off the ID= in the gene ID
  mutate(combined = (AR*21 + CA*17)/(17+21)) # combined average across all pops included

# what does the ancestry distribution across genes look like?
ggplot(genes0, aes(x = CA, y = AR)) +
  geom_point() +
  ggtitle("mean ancestry for all genes")
# 168 genes have no ancestry calls -- look at later whether they're candidates
hist(apply(genes0[,c("AR", "CA")], 1, mean))

# how many genes are outliers?
sapply(outlier_types, function(x)
  table(genes[ , x]))


# assess overlap with QTLs. top hits? enrichment?
# first varroa qtls
hygeine_QTL <- read.table("../data/honeybee_genome/Amel_4.5/QTLs/Harpur_2019_social_immunity_hygeine/DatasetS1_PreviousAssociated_hygeine_QTLs.txt",
                          header = T, sep = "\t", stringsAsFactors = F) %>%
  mutate(Scaffold_Start = ifelse(is.na(Scaffold_Start), NA, paste0("Group", Scaffold_Start))) %>%
  # note: typo, should be Oxley 2010 (2008 was a study on QTLs for worker sterility)
  mutate(Scaffold_End = ifelse(is.na(Scaffold_End), NA, paste0("Group", Scaffold_End)))
grooming_QTL <- read.table("Amel4.5_results/grooming_QTL_Arechavaleta-Velasco_2012.txt",
                           header = T, sep = "\t", stringsAsFactors = F) %>%
  mutate(Outlier_type = "QTL")
# Using the conservatively chosen confidence intervals forHyg1,2and3, a  total  of 339  known and  putative  geneswere  identified.  Of  these,  218  had  orthologs  with  GeneOntology  annotations  (The  Gene  Ontology  Consortium2000).
# v3 snp file (with seq. next to it? 250bp..) plus positions of putative QTLs Spoetter 2012: https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1755-0998.2011.03106.x

# there's also uncap1, uncap2, and rem1 from Oxley et al. 2010 and hyg1-3 but I don't easily have markers: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2010.04569.x
# ok hyg1-3 are in the harpur list, only the hyg 1 on chr2 was sig.
# tentative other QTLs from that study: "2 suggested QTLs for uncapping & 1 for removal"
#Component 	Uncapping 	Uncapping 	Removal
#Position 	Chr 9, 218cM 	Chr 16, 0cM 	Chr 10, 98cM
#Nearest marker 	AT128 	K1601 	AC074 (microsats from Solignac 2007; 155 markers, mean marker spacing 27.7cM)
# also first half of chromosome 7 has a sig. QTL for suppression of reproduction of varroa; 
# epistatic effects in Gotland Varroa tolerant honey bees (Behrens 2011 Ecol Evol.)
varroa_QTL0 <- bind_rows(hygeine_QTL, grooming_QTL) %>%
  mutate(chr = as.numeric(Chromosome))

# list of genes and regions associated with hygeine from Harpur 2019
hygeine_genes <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/TableS3_hygeine_associated_gene_list.csv",
                            header = F, sep = "\t", stringsAsFactors = F)$V1
hygeine_regions <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/TableS2_regions_selected_and_associated_with_hygeine.csv",
                              header = T, sep = "\t", stringsAsFactors = F)

varroa_scaffolds_on_chr <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/scaffolds_on_chr.txt",
                                      sep = "\t", header = T, stringsAsFactors = F)

# any overlap with my outlier gene list? none with these 73 genes
genes$hygeine_Harpur2019 <- sapply(genes$ID, function(i) i %in% hygeine_genes)
#genes %>% 
#  filter(hygeine_Harpur2019) %>% 
#  View() # no overlap
table(genes$hygeine_Harpur2019)
genes %>% # very small number, only 73 genes, this could all be noise.
  tidyr::gather(., key = "zone", value = "A_ancestry", c("combined", "CA", "AR")) %>%
  ggplot(., aes(y = A_ancestry, x = zone, fill = hygeine_Harpur2019)) + 
  geom_boxplot() +
  ggtitle("A ancestry in hygeine genes (Harpur 2019) vs. all other genes")

# next step: put all sites with ancestry calls on chromosomes in the same order as QTLs listed
scaff_to_chr_pos <- function(scaffold, pos, mapg = varroa_scaffolds_on_chr){
  pos_chr <- ifelse(mapg[mapg$scaffold == scaffold, "orientation"] == "-",
                    mapg[mapg$scaffold == scaffold, "stop"] - (pos - 1), # reverse -
                    mapg[mapg$scaffold == scaffold, "start"] + (pos - 1)) # forward or unoriented +/0
  return(pos_chr)
}
A_AR_CA$chr_start = sapply(1:nrow(A_AR_CA), function(x) scaff_to_chr_pos(scaffold = A_AR_CA[x, "scaffold"], pos = A_AR_CA[x, "start"]))
A_AR_CA$chr_end = sapply(1:nrow(A_AR_CA), function(x) scaff_to_chr_pos(scaffold = A_AR_CA[x, "scaffold"], pos = A_AR_CA[x, "end"]))
A_AR_CA$chr = as.numeric(sapply(A_AR_CA$scaffold, function(s) substr(strsplit(s, split = "[.]")[[1]][1], 6, 100)))
varroa_QTL0$chr = as.numeric(apply(varroa_QTL0, 1, function(s) ifelse(is.na(s["chr"]), 
                                                                      as.numeric(substr(strsplit(s["Scaffold_Start"], split = "[.]")[[1]][1], 6, 100)),
                                                                      s["chr"])))
varroa_QTL0$chr_start = sapply(1:nrow(varroa_QTL0), function(x) ifelse(is.na(varroa_QTL0[x, "Start"]), 
                                                                       scaff_to_chr_pos(scaffold = varroa_QTL[x, "Scaffold_Start"], 
                                                                                        pos = varroa_QTL0[x, "Scaffold_Pos_Start"]),
                                                                       varroa_QTL0[x, "Start"]))
varroa_QTL0$chr_end = sapply(1:nrow(varroa_QTL0), function(x) ifelse(is.na(varroa_QTL0[x, "End"]), 
                                                                     scaff_to_chr_pos(scaffold = varroa_QTL0[x, "Scaffold_End"], 
                                                                                      pos = varroa_QTL0[x, "Scaffold_Pos_End"]),
                                                                     varroa_QTL0[x, "End"]))

QTL_names <- unique(varroa_QTL0[c("Reference", "Outlier_type", "Name")])
QTL_names$name <- c("Varroa-Sensitive Hygiene putative QTL", "Varro-Sensitive Hygiene GWAS SNP", "Varroa-Sensitive Hygiene QTL", "Hygiene QTL", "Self-Grooming QTL")
QTL_names$source <- c("Spoetter et al. 2012", "Spoetter et al. 2016", "Tsuruda et al. 2012", "Oxley et al. 2008", "Arechavaleta-Velasco et al. 2012")
varroa_QTL <- left_join(varroa_QTL0, QTL_names, by = c("Reference", "Outlier_type", "Name"))

# plot oultiers on their scaffolds with genes under outliers shown. 
# I also want a whole-genome view to make sure the genes 'out of range' near the edges of scaffolds (past ancestry calls) 
# aren't likely ancestry outliers
A_AR_CA %>%
  mutate(CA = -1*CA) %>% # flip CA axis
  #rename(FDR_shared_low = "\n                       FDR_shared_low") %>% # rename one weird column
  arrange(chr) %>% # sort by chromosome order
  tidyr::gather(., "zone", "A_ancestry", c("CA", "AR")) %>%
  mutate(FDR = apply(., 1, function(x) ifelse(x["zone"] == "CA", 
                                              min(x[c("FDR_shared_high", "FDR_shared_low", 
                                                      "FDR_CA_high", "FDR_CA_low")], na.rm = T),
                                              min(x[c("FDR_shared_high", "FDR_shared_low", 
                                                      "FDR_AR_high", "FDR_AR_low")], na.rm = T)))) %>%
  ggplot(aes(x = chr_start, y = A_ancestry, color = FDR)) + # I should get true mid position, not start, but ok for now
  geom_point(size = .2) +
  facet_wrap(~chr, scales = "free_x") +
  ggtitle("African ancestry frequency in Argentina (top) and California (mirrored, bottom)")
ggsave("plots/A_frequency_mirrored_plot_AR_CA_FDR_whole_genome.png",
       height = 10, width = 12, units = "in", device = "png")

# make a plot 'wide' format 'manhattan style':
# manhattan style plot


varroa_QTL_cumulative <- varroa_QTL %>%
  left_join(., chr_lengths, by = "chr") %>%
  mutate(cum_start = chr_start + chromosome_start,
         cum_end = chr_end + chromosome_start) %>% # get cumulative chromosome positions
  mutate(Name = ifelse(Name == "", 
                       ifelse(Outlier_type == "NCBI_QTLs", "Putative Varroa Removal (unpub.)", 
                              Outlier_type), Name)) %>%
  mutate(QTL_length = cum_end - cum_start)


# get QTLs:
QTL = qtl %>%
  left_join(., chr_lengths, by = c("AMEL_CHR"="scaffold")) %>%
  # start is the lower one, end is the higher endpoint
  mutate(start = apply(., 1, function(x) min(as.integer(x[["Start"]]), as.integer(x[["End"]]))),
         end = apply(., 1, function(x) max(as.integer(x[["Start"]]), as.integer(x[["End"]]))))

# simple plot of outliers, whole genome
ggplot() + # raster looks pretty terrible -- I could plot instead every 5th or 10th non-sig point or s.t.
  #ggrastr::geom_point_rast(data = A_AR_CA_cumulative %>%
  #             filter(is.na(FDR)), # plot grey points first
  #           aes(x = cum_pos, y = A_ancestry, 
  #               color = color_by), size = .02,
  #           raster.dpi = 600) +
  geom_point(data = (A_AR_CA_cumulative %>% # plot grey points first; every 10th point only
                       filter(is.na(FDR)))[c(T, rep(F, 9)), ], 
             aes(x = cum_pos, y = A_ancestry, 
                 color = color_by), size = .02) +
  geom_point(data = A_AR_CA_cumulative %>%
               filter(!is.na(FDR)), # then plot sig points on top 
             aes(x = cum_pos, y = A_ancestry, 
                 color = color_by), size = .02) +
  xlab("Position (bp)") +
  ylab("mean African ancestry") +
  #scale_colour_manual(name = NULL,
  #                    values = c("0.01"="red", "0.05"="orange", "0.10"="skyblue", 
  #                               "n.s. - even chr"="darkgrey", 
  #                               "n.s. - odd chr"="grey"),
  #                    limits = c("0.01", "0.05", "0.10"),
  #                    labels = c("0.01 FDR", "0.05 FDR", "0.10 FDR")
  #) + 
  geom_segment(data = left_join(rename(outliers_all, scaffold = chr) %>%
                                  rename(., length_outlier = length), 
                                chr_lengths, by = "scaffold") %>%
                 #mutate(cum_start = start + chr_start,
                 #      cum_end = end + chr_start)
                 filter(outlier_type == "high_shared"),
               aes(x = start + chr_start - buffer_4_visibility, 
                   xend = end + chr_start + buffer_4_visibility, # add 50kb for visibility only
                   y = 0.7, yend = 0.7, color = min_FDR),
               lwd = 4) +
  # add in QTLs:
  geom_segment(data = QTL,
               aes(x = start + chr_start - buffer_4_visibility, # add 50kb for visibility only
                   xend = end + chr_start + buffer_4_visibility, 
                   y = 0.7, yend = 0.7, color = Publication),
               lwd = 4) +
  scale_x_continuous(label = chr_lengths$chr_n, breaks = chr_lengths$chr_mid) +
  #theme(legend.position = "none") +
  theme_classic() +
  facet_wrap(~zone_pretty, nrow = 2, ncol = 1) + 
  ggtitle("Ancestry outliers across whole genome")
ggsave("plots/A_frequency_plot_AR_CA_FDR_whole_genome_wide.png",
       height = 5, width = 10, units = "in", device = "png")
ggsave("../../bee_manuscript/figures/A_frequency_plot_AR_CA_FDR_whole_genome_wide.pdf",
       height = 5, width = 10, units = "in", device = "pdf")





# to do later: better idea is to color the shared outliers and to draw threshold lines for the ind zone outliers
# also later I want to make plots of frequency across latitude for all the major outliers
# AND make C-M outlier plots too (maybe plot first before thinking about FDR)

# mirrored plot example:
#par(mfrow=c(2,1))
#Make the plot
#par(mar=c(0,5,3,3))
#plot(density(x1) , main="" , xlab="", ylim=c(0,1) , xaxt="n", las=1 , col="slateblue1" , lwd=4 )
#par(mar=c(5,5,0,3))
#plot(density(x2) , main="" , xlab="Value of my variable", ylim=c(1,0) , las=1 , col="tomato3" , lwd=4)  

###-------------------------------------------------------###
# read in genes list with beebase IDs:
beebase_outlier_genes <- read.table("results/DAVID_results_genes_0.1FDR_combined_with_beebase.csv",
                                    sep = ",", header = T)
#"([A-Z]+|[a-z]+)[=]"
remove_id <- function(x, link) stringr::str_replace(x, paste0("(.+)[", link, "]"), "") # turns ID=3253 into just 3253
get_id <- function(x, link) stringr::str_replace(stringr::str_extract(x, paste0("(.+)[", link, "]")), paste0("[", link, "]"), "") # returns "ID"
get_gene_id_info <- function(x, split = ";", link = "="){ # takes in dictionary and returns named vector
  v = strsplit(x, split = split)[[1]] # turns a dictionary string e.g. "ID=89;Gene=20" into a vector, split by ;
  y = remove_id(v, link = link) # returns just the values: 89, 20
  names(y) = get_id(v, link = link) # returns just the IDs: ID, Gene
  return(y)
}

#gene_id_cols <- c("ID", "Dbxref", "Name", "gbkey", "gene", "gene_biotype")

get_gene_id_info(x = "BEEBASE:2,GeneID:LOC213", split = ",", link = ":")
genes_combined2 = do.call(bind_rows,
                          lapply(genes_combined$gene_info,
                                 function(x) get_gene_id_info(x, split = ";", link = "="))) %>%
  cbind(dplyr::select(genes_combined, -c(gene_info, gene)), .) %>%
  cbind(., #dplyr::select(., -Dbxref)
        do.call(bind_rows,
                lapply(.$Dbxref, function(x)
                  get_gene_id_info(x, split = ",", link = ":")))) %>%
  left_join(., beebase_outlier_genes[ , c("BEEBASE_ID", "Gene.Name", "Related.Genes", "Species")], 
            by = c("BEEBASE"="BEEBASE_ID")) %>%
  rename(DAVID_gene_name = Gene.Name, 
         DAVID_related_genes = Related.Genes, 
         DAVID_species = Species) %>%
  dplyr::arrange(desc(outlier_type), scaffold, start) %>%
  mutate(FDR_AR_high = as.numeric(FDR_AR_high),
         FDR_CA_high = as.numeric(FDR_CA_high),
         FDR_AR_low = as.numeric(FDR_AR_low))
table(!is.na(genes_combined2$BEEBASE))
dim(beebase_outlier_genes)
beebase_outlier_genes$BEEBASE_ID[!(beebase_outlier_genes$BEEBASE_ID %in% genes_combined2$BEEBASE)]
genes_combined2$BEEBASE[!(genes_combined2$BEEBASE %in% beebase_outlier_genes$BEEBASE_ID) & !is.na(genes_combined2$BEEBASE)]

# write genes and DAVID functional information to file
write.table(dplyr::select(genes_combined2, -Dbxref), 
            "results/genes_0.1FDR_combined_with_DAVID_functions.txt",
            sep = "\t",
            col.names = T, row.names = F)  
