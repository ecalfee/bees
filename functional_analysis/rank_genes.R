
# this script ranks genes & plots them by evidence for selection
# load mean ancestry for genes
genes0 <- read.table("results/genes_mean_AR_CA_ancestry.bed", stringsAsFactors = F, 
                       na.strings = c("."), sep = "\t")
colnames(genes0) <- c("scaffold", "start", "end", "gene_list", "gff3_type", "gene_ID", "AR", "CA")

# load FDR for genes -- which genes are outliers?
outlier_types <- c("shared_high", "AR_high", "CA_high", # no CA_low outliers
                   "shared_low", "AR_low")
outlier_genes <- do.call(cbind,
                         lapply(outlier_types, function(x) 
                            read.table(paste0("results/outliers_", x, "_genes_mean_AR_CA_ancestry.bed"),
                                       stringsAsFactors = F, na.strings = c("."))$V9))
colnames(outlier_genes) <- outlier_types
genes <- cbind(genes0, outlier_genes) %>%
  mutate(ID = substr(gene_ID, 4, 100)) %>% # take off the ID= in the gene ID
  mutate(combined = (AR*21 + CA*17)/(17+21)) # combined average across all pops included

# load mean ancestry all snps with ancestry calls genomewide
A_AR_CA <- read.table("../local_ancestry/results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/anc/mean_ancestry_AR_CA_included.bed", 
              sep = "\t", stringsAsFactors = F, na.strings = c("."))
colnames(A_AR_CA) <- c("scaffold", "start", "end", "snp_id", "AR", "CA", 
                       "FDR_shared_high", "FDR_AR_high", "FDR_CA_high", "
                       FDR_shared_low", "FDR_AR_low", "FDR_CA_low")


# load false-discovery-rate cutoffs based on MVN simulations
FDRs <- read.table("../local_ancestry/results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", 
                   stringsAsFactors = F, header = T, sep = "\t")

# what does the ancestry distribution across genes look like?
ggplot(genes0, aes(x = CA, y = AR)) +
  geom_point() +
  ggtitle("mean ancestry for all genes")
# 168 genes have no ancestry calls -- look at later whether they're candidates
hist(apply(genes0[,c("AR", "CA")], 1, mean))

# how many genes are outliers?
sapply(outlier_types, function(x)
       table(genes[ , x]))

# print list of candidate genes
# keep clustered by scaffold, but sort priority within scaffold based on mean_AR_CA ancestry
# high
genes %>%
  filter(!is.na(shared_high) | !is.na(AR_high) | !is.na(CA_high)) %>%
  arrange(desc(combined)) %>%
  dplyr::select(scaffold, start, end, gene_list, gene_ID, AR, CA, combined, shared_high, AR_high, CA_high) %>%
  write.table(., "results/outlier_genes_list_high_A_ancestry.txt", na = "n.s.",
              sep = "\t", quote = F, col.names = T, row.names = F)
# count how many genes per scaffold (i.e. together)?
genes %>%
  filter(!is.na(shared_high) | !is.na(AR_high) | !is.na(CA_high)) %>%
  group_by(scaffold) %>%
  summarise(n = n()) %>%
  write.table(., "results/outlier_genes_count_per_scaffold_high_A_ancestry.txt",
              sep ="\t", quote = F, col.names = T, row.names = F)

# low
genes %>%
  filter(!is.na(shared_low) | !is.na(AR_low)) %>%
  arrange(combined) %>%
  mutate(CA_low = NA) %>% # no CA low outliers
  dplyr::select(scaffold, start, end, gene_list, gene_ID, AR, CA, combined, shared_low, AR_low, CA_low) %>%
  write.table(., "results/outlier_genes_list_low_A_ancestry.txt", na = "n.s.",
              sep = "\t", quote = F, col.names = T, row.names = F)
genes %>%
  filter(!is.na(shared_low) | !is.na(AR_low)) %>%
  group_by(scaffold) %>%
  summarise(n = n()) %>%
  write.table(., "results/outlier_genes_count_per_scaffold_low_A_ancestry.txt",
              sep ="\t", quote = F, col.names = T, row.names = F)

# assess overlap with QTLs. top hits? enrichment?
# first varroa qtls
varroa_QTL <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/DatasetS1_PreviousAssociated_hygeine_QTLs.txt",
                         header = T, sep = "\t", stringsAsFactors = F)

# list of genes and regions associated with hygeine from Harpur 2019
hygeine_genes <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/TableS3_hygeine_associated_gene_list.csv",
                            header = F, sep = "\t", stringsAsFactors = F)$V1
hygeine_regions <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/TableS2_regions_selected_and_associated_with_hygeine.csv",
                              header = T, sep = "\t", stringsAsFactors = F)

varroa_scaffolds_on_chr <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/scaffolds_on_chr.txt",
                                      sep = "\t", header = T, stringsAsFactors = F)

# any overlap with my outlier gene list? none with these 73 genes
genes$hygeine_Harpur2019 <- sapply(genes$ID, function(i) i %in% hygeine_genes)
genes %>% filter(hygeine_Harpur2019) %>% View() # no overlap

# next step: put all sites with ancestry calls on chromosomes in the same order as QTLs listed
scaff_to_chr_pos <- function(scaffold, pos, mapg = varroa_scaffolds_on_chr){
  pos_chr <- ifelse(mapg[mapg$scaffold == scaffold, "orientation"] == "-",
                    mapg[mapg$scaffold == scaffold, "stop"] - (pos - 1), # reverse -
                    mapg[mapg$scaffold == scaffold, "start"] + (pos - 1)) # forward or unoriented +/0
  return(pos_chr)
}
A_AR_CA$chr_start = apply(A_AR_CA, 1, function(x) scaff_to_chr_pos(scaffold = x["scaffold"], pos = x["start"]))
A_AR_CA$chr_end = apply(A_AR_CA, 1, function(x) scaff_to_chr_pos(scaffold = x["scaffold"], pos = x["end"]))
A_AR_CA$chr = sapply(A_AR_CA$scaffold, function(s) substr(strsplit(s, split = "[.]")[[1]][1], 6, 100))
# ***************** oops! something is not working!
scaff_to_chr_pos(scaffold = "Group16.8", pos = 2)
sapply(unique(A_AR_CA$scaffold), function(x) x %in% unique(varroa_scaffolds_on_chr$scaffold))
sum(is.na(A_AR_CA$pos))



# plot oultiers on their scaffolds with genes under outliers shown. 
# I also want a whole-genome view to make sure the genes 'out of range' near the edges of scaffolds (past ancestry calls) 
# aren't likely ancestry outliers

str(A_AR_CA)


