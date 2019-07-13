library(dplyr)
library(ggplot2)

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
hygeine_QTL <- read.table("../data/honeybee_genome/QTLs/Harpur_2019_social_immunity_hygeine/DatasetS1_PreviousAssociated_hygeine_QTLs.txt",
                         header = T, sep = "\t", stringsAsFactors = F) %>%
  mutate(Scaffold_Start = ifelse(is.na(Scaffold_Start), NA, paste0("Group", Scaffold_Start))) %>%
  mutate(Scaffold_End = ifelse(is.na(Scaffold_End), NA, paste0("Group", Scaffold_End)))
grooming_QTL <- read.table("results/grooming_QTL_Arechavaleta-Velasco_2012.txt",
                        header = T, sep = "\t", stringsAsFactors = F)
# there's also uncap1, uncap2, and rem1 from Oxley et al. 2010 and hyg1-3 but I don't easily have markers: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2010.04569.x
# ok hyg1-3 are in the harpur list, only the hyg 1 on chr2 was sig.
# tentative other QTLs from that study: "2 suggested QTLs for uncapping & 1 for removal"
#Component 	Uncapping 	Uncapping 	Removal
#Position 	Chr 9, 218cM 	Chr 16, 0cM 	Chr 10, 98cM
#Nearest marker 	AT128 	K1601 	AC074 (microsats from Solignac 2007; 155 markers, mean marker spacing 27.7cM)
# also first half of chromosome 7 has a sig. QTL for suppression of reproduction of varroa; 
# epistatic effects in Gotland Varroa tolerant honey bees (Behrens 2011 Ecol Evol.)
varroa_QTL <- bind_rows(hygeine_QTL, grooming_QTL) %>%
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
varroa_QTL$chr = as.numeric(apply(varroa_QTL, 1, function(s) ifelse(is.na(s["chr"]), 
                                                                      as.numeric(substr(strsplit(s["Scaffold_Start"], split = "[.]")[[1]][1], 6, 100)),
                                                         s["chr"])))
varroa_QTL$chr_start = sapply(1:nrow(varroa_QTL), function(x) ifelse(is.na(varroa_QTL[x, "Start"]), 
                                                                      scaff_to_chr_pos(scaffold = varroa_QTL[x, "Scaffold_Start"], 
                                                                                       pos = varroa_QTL[x, "Scaffold_Pos_Start"]),
                                                                      varroa_QTL[x, "Start"]))
varroa_QTL$chr_end = sapply(1:nrow(varroa_QTL), function(x) ifelse(is.na(varroa_QTL[x, "End"]), 
                                                                     scaff_to_chr_pos(scaffold = varroa_QTL[x, "Scaffold_End"], 
                                                                                      pos = varroa_QTL[x, "Scaffold_Pos_End"]),
                                                                     varroa_QTL[x, "End"]))

# plot oultiers on their scaffolds with genes under outliers shown. 
# I also want a whole-genome view to make sure the genes 'out of range' near the edges of scaffolds (past ancestry calls) 
# aren't likely ancestry outliers
A_AR_CA %>%
  mutate(CA = -1*CA) %>% # flip CA axis
  rename(FDR_shared_low = "\n                       FDR_shared_low") %>% # rename one weird column
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
# Compute chromosome size

chr_lengths <- varroa_scaffolds_on_chr %>%
  mutate(chr = as.numeric(substr(chr, 4, 100))) %>%
  group_by(chr) %>%
  arrange(chr) %>%
  summarise(length = max(stop)) %>%
  mutate(chromosome_start = cumsum(length) - length,
         chromosome_end = cumsum(length)) %>%
  mutate(chromosome_midpoint = (chromosome_start + chromosome_end)/2)

varroa_QTL_cumulative <- varroa_QTL %>%
  left_join(., chr_lengths, by = "chr") %>%
  mutate(cum_start = chr_start + chromosome_start,
         cum_end = chr_end + chromosome_start) %>% # get cumulative chromosome positions
  mutate(Name = ifelse(Name == "", 
                       ifelse(Outlier_type == "NCBI_QTLs", "Putative Varroa Removal (unpub.)", 
                              Outlier_type), Name)) %>%
  mutate(cum_start = ifelse(Outlier_type == "GWAS_SNPs", cum_start - 20000, cum_start),
         cum_end = ifelse(Outlier_type == "GWAS_SNPs", cum_end + 20000, cum_end)) # just to visualize

# get approximate cumulative positions for all points with ancestry calls
A_AR_CA_cumulative <- A_AR_CA %>%
  mutate(pos = (chr_start + chr_end)/2) %>% # proxy ok for plotting; regions are short
  left_join(., chr_lengths, by = "chr") %>%
  mutate(pos_cum = pos + chromosome_start) %>%
  #mutate(CA = -1*CA) %>% # flip CA axis
  rename(FDR_shared_low = "\n                       FDR_shared_low") %>% # rename one weird column
  arrange(chr) %>% # sort by chromosome order
  tidyr::gather(., "zone", "A_ancestry", c("CA", "AR")) %>%
  mutate(FDR = apply(., 1, function(x) ifelse(x["zone"] == "CA", 
                                              min(x[c("FDR_shared_high", "FDR_shared_low", 
                                                      "FDR_CA_high", "FDR_CA_low")], na.rm = T),
                                              min(x[c("FDR_shared_high", "FDR_shared_low", 
                                                      "FDR_AR_high", "FDR_AR_low")], na.rm = T)))) %>%
  mutate(color_by = ifelse(FDR == "NA", ifelse((chr %% 2 == 0), # even chromosomes different color
                                               "n.s. - even chr", "n.s. - odd chr"), FDR))

# plot with QTLs for varroa on top of 
ggplot() +
  geom_point(data = A_AR_CA_cumulative, aes(x = pos_cum, y = A_ancestry, 
                                            color = color_by), size = .1) +
  xlab("bp position on chromosome") +
  ylab("mean African ancestry") +
  scale_colour_manual(values = c("red", "orange", "skyblue", 
                                 "limegreen", "darkgreen", "olivedrab", 
                                 "darkgrey", "grey",
                                 "lightseagreen", "lightgreen")) + 
  scale_x_continuous(label = chr_lengths$chr, breaks = chr_lengths$chromosome_midpoint) +
  #theme(legend.position = "none") +
  facet_wrap(~zone, nrow = 2, ncol = 1) + 
  geom_segment(data = varroa_QTL_cumulative, 
               aes(x=cum_start, xend=cum_end, y=0.68, yend=0.68, color = Name),
               lwd = 4) +
  ggtitle("A ancestry in both hybrid zones, showing varroa-resistance QTLs in green")
ggsave("plots/A_frequency_plot_AR_CA_FDR_whole_genome_wide_VarroaQTL.png",
       height = 5, width = 10, units = "in", device = "png")
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
  scale_colour_manual(values = c("red", "orange", "skyblue", 
                                 "darkgrey",
                                 "lightseagreen", "lightgreen")) + 
  #theme(legend.position = "none") +
  facet_wrap(~zone, nrow = 2, ncol = 1) + 
  geom_segment(data = filter(varroa_QTL_cumulative, chr == 1), 
               aes(x=cum_start, xend=cum_end, y=0.68, yend=0.68, color = Name),
               lwd = 4) +
  ggtitle("A ancestry in both hybrid zones CHR1 showing varroa-resistance QTLs in green")
ggsave("plots/A_frequency_plot_AR_CA_FDR_chr1_wide_VarroaQTL.png",
       height = 5, width = 10, units = "in", device = "png")


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


