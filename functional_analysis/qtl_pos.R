library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Get QTL positions re-mapped to new genome HAv3.1
# Spoetter 2012
# chr data for HAv3.1
chr_lengths <- read.table("../data/honeybee_genome/chr.list", sep = "\t", stringsAsFactors = F)
colnames(chr_lengths) <- c("scaffold", "length", "chr_group", "chr_lg")

# maximum mapping evalue to keep a blastn hit:
max_eval = .01
# spoetter qtls on Amel3 coordinates
spoetter <- read.table("../data/Spoetter_2012_QTLs/spoetter_2012_Amel3_start_end.csv",
                       sep = ",", header = T, stringsAsFactors = F) %>%
  left_join(., rename(chr_lengths, chr = chr_group), by = "chr") %>%
  mutate(QTL = 1:nrow(.)) %>%
  dplyr::select(scaffold, start, end, QTL, chr_lg) %>%
  #dplyr::select(chr, scaffold, start, end, QTL, chr_lg) %>%
  rename(chr = scaffold, LG = chr_lg)
  #rename(LG = chr_lg)
# 44k snp chip on amel3 coordinates
snps_44k <- read.table("../data/Spoetter_2012_QTLs/SNP_Info_HoneyBeeSNPAssay.csv",
                       sep = ",", header = T, stringsAsFactors = F)
head(snps_44k)
snps_44k %>%
  arrange(linkage.group, Cooridnate) %>%
  write.table("../data/Spoetter_2012_QTLs/SNPs_44k.txt",
              quote = F, sep = "\t", col.names = T, row.names = F)
# separately I ran a local blastn for short seq. alignment
# to get positions for left.flanking seq
# load mapped coordinates (using blastn) to HAv3.1
Hav3.1_coord_44k <- do.call(bind_rows, 
                            lapply(1:16, function(i)
                              read.table(paste0("results/SNPs_44k_blast_pos/Group", i, "_blast_results.txt"),
                                         header = F, stringsAsFactors = F) %>%
                                data.table::setnames(c("lsid", "scaffold", "start", "end",
                                                       "perc_identical", "evalue", "bitscore", "score")))) %>%
  dplyr::arrange(., lsid, evalue) %>%
  filter(!duplicated(lsid)) # get rid of all blast matches but highest score (lowest evalue)
snps_44k_coord <- snps_44k %>%
  left_join(., chr_lengths, by = c("linkage.group"="chr_group")) %>%
  rename(pos = Cooridnate) %>%
  rename(Amel3_scaffold = scaffold) %>%
  dplyr::select(lsid, Amel3_scaffold, linkage.group, pos) %>%
  left_join(., Hav3.1_coord_44k, by = "lsid") %>%
  filter(!is.na(start) & evalue <= max_eval) %>% # filter out snps with no match or low match
  rename(chr = Amel3_scaffold,
         start_HAv3.1 = start,
         end_HAv3.1 = end,
         scaffold_HAv3.1 = scaffold) #%>%
  #mutate(start = pos - 1, 
        # end = pos) %>%
  #dplyr::select(chr, start, end, 
   #             scaffold_HAv3.1, start_HAv3.1, end_HAv3.1,
    #            perc_identical, lsid)
# now which snps fall within qtls?
qtl_spoetter_snps <- do.call(rbind,
                        lapply(1:nrow(spoetter), function(i){
  filter(snps_44k_coord, chr == spoetter$chr[i] &
           pos >= spoetter$start[i] &
           pos <= spoetter$end[i]) %>%
                            mutate(qtl = spoetter$QTL[i])
}))
#qtl_spoetter <- bedr(
#  engine = "bedtools", 
 # input = list(a = snps_44k_coord,
  #             b = spoetter), 
#  method = "map", 
# params = "-g ../data/honeybee_genome/chr.lengths -c 4 -o distinct",
 # check.chr = F
#) %>%
 # data.table::setnames(c(colnames(snps_44k_coord), "qtl"))

table(paste0(qtl_spoetter_snps$qtl, qtl_spoetter_snps$scaffold_HAv3.1))
qtl_spoetter_snps %>%
  ggplot(., aes(x = pos/10^6, y = end_HAv3.1/10^6, color = scaffold_HAv3.1)) +
  geom_point() +
  facet_wrap(~paste("QTL", qtl, "from", chr)) +
  xlab("Old marker position (Mbp) on Amel3") +
  ylab("New mapped position (Mbp) HAv3.1") +
  labs(color = "New scaffold HAv3.1")
ggsave("plots/remap_spoetter_2012_QTLs_to_HAv3.1.png",
       height = 6, width = 7.5, dpi = 600, device = "png")
# a few small segments map to the wrong chromosome, but otherwise no big gaps or mapping outside target region,
# so I'll just get approximate endpoints with min and max of the end (end not start because mapped left flanking region)
qtl_spoetter <- qtl_spoetter_snps %>%
  filter(scaffold_HAv3.1 == chr) %>%
  group_by(qtl, chr) %>%
  summarise(start = min(end_HAv3.1),
            end = max(end_HAv3.1)) %>%
  mutate(citation = "Spoetter et al. 2012",
         positions_from = "SNPs within QTL mapping interval from 44K SNP chip published in Spoetter et al. 2012 (Supplement, left flanking sequences) mapped to HAv3.1 with BLASTn",
         phenotype = "Varroa Sensitive Hygiene") %>%
  left_join(., spoetter %>%
              mutate(comments = paste0("Amel3 mapping interval (Table 1 Spoetter et al. 2012) ", LG,
                                       ":", start, "-", end)) %>%
              dplyr::select(chr, QTL, comments), by = c("chr", "qtl"="QTL")) %>%
  rename(scaffold = chr) %>%
  left_join(., chr_lengths[ , c("scaffold", "chr_group")], by = "scaffold") %>%
  mutate(QTL_name = paste0("unnamed_hygeine", qtl))
  
  

# Oxley
oxley <- read.table("results/blast_results_oxley_2010_markers.txt",
                    header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("marker_id", "scaffold", "start", "end",
                         "perc_identical", "evalue", "bitscore", "score")) %>%
  dplyr::arrange(., marker_id, evalue) %>%
  filter(!duplicated(marker_id)) %>% # get rid of all blast matches but highest score (lowest evalue)
  left_join(., chr_lengths[ , c("scaffold", "chr_group")], by = "scaffold") %>%
  tidyr::separate(., col = "marker_id", into = c("QTL_name", "marker_id"), sep = "_") %>%
  tidyr::separate(., col = "marker_id", into = c("marker", "orig_chr", "orig_pos", "primer"), sep = "-") %>%
  tidyr::pivot_longer(., col = c("start", "end"), names_to = "pos_type", values_to = "pos") %>%
  filter(evalue <= max_eval) %>% # meets mapping threshold
  group_by(QTL_name, marker, scaffold) %>%
  summarise(start = min(pos), end = max(pos)) %>%
  left_join(., chr_lengths[ , c("scaffold", "chr_group")], by = "scaffold") %>%
  mutate(phenotype = "Hygienic behaviour",
         comments = "Nearest marker",
         citation = "Oxley et al. 2010",
  positions_from = "P1 and P2 probe sequences published in Solignac et al. 2007 (Table S2) mapped to HAv3.1 with BLASTn")

# A-Vel
av <- read.table("results/blast_results_arechavaleta-velasco_2012_markers.txt",
                 header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("marker", "scaffold", "start", "end",
                         "perc_identical", "evalue", "bitscore", "score")) %>%
  dplyr::arrange(., marker, evalue) %>%
  filter(!duplicated(marker)) %>% # get rid of all blast matches but highest score (lowest evalue)
  filter(evalue <= max_eval) %>% # meets mapping threshold
  left_join(., chr_lengths[ , c("scaffold", "chr_group")], by = "scaffold") %>%
  mutate(phenotype = "Grooming",
         QTL_name = "unnamed_grooming1",
         citation = "Arechavaleta-Velasco et al. 2012",
         positions_from = "Probe sequences from Arechavaleta-Velasco et al. 2012 (Table S1) mapped to HAv3.1 with BLASTn",
         comments = ifelse(marker == "5_1644271", 
                           "LOD-1.5 support interval; marker closest to peak", 
                           "LOD-1.5 support interval"))

tsuruda <- read.table("results/blast_results_tsuruda_2012.txt",
                   header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("marker", "scaffold", "start", "end",
                         "perc_identical", "evalue", "bitscore", "score")) %>%
  dplyr::arrange(., marker, evalue) %>%
  filter(!duplicated(marker)) %>% # get rid of all blast matches but highest score (lowest evalue)
  filter(evalue <= max_eval) %>% # meets mapping threshold
  left_join(., chr_lengths[ , c("scaffold", "chr_group")], by = "scaffold") %>%
  mutate(phenotype = "Varroa sensitive hygiene",
         QTL_name = ifelse(chr_group == "Group1", "unnamed_VSH1",
                           ifelse(chr_group == "Group9", "unnamed_VSH2", NA)),
         citation = "Tsuruda et al. 2012",
         positions_from = "Probe sequences from Tsuruda et al. 2012 (Table S3) mapped to HAv3.1 with BLASTn",
         comments = ifelse(marker == "9_9224292", 
                           "LOD-1.5 support interval; marker closest to peak", 
                           "LOD-1.5 support interval"))
sting <- read.table("../data/QTL_markers_Hunt_lab/sting1-3_markers.tsv",
                    sep = "\t", header = T, stringsAsFactors = F) %>%
  rename(comments = note)
hunt0 <- read.table("results/blast_results_hunt.txt",
                      header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("marker", "scaffold", "start", "end",
                         "perc_identical", "evalue", "bitscore", "score")) %>%
  dplyr::arrange(., marker, evalue) %>%
  filter(!duplicated(marker)) %>% # get rid of all blast matches but highest score (lowest evalue)
  filter(evalue <= .1) %>% # meets mapping threshold
  rename(marker_id = marker) %>%
  tidyr::separate(., col = "marker_id", into = c("marker_type", "marker"), remove = F, sep = "_") %>%
  tidyr::pivot_longer(., cols = c("start", "end"), names_to = "start_or_end", values_to = "pos") %>%
  group_by(marker, scaffold) %>%
  summarise(start = min(pos),
            end = max(pos)) %>%
  left_join(., chr_lengths[ , c("scaffold", "chr_group")], by = "scaffold") %>%
  left_join(sting, ., by = "marker") %>%
  mutate(phenotype = "Defense response",
         seq = ifelse(seq == "", "-", seq),
         primer1 = ifelse(primer1 == "", "-", primer1),
         primer2 = ifelse(primer2 == "", "-", primer2),
         sequence = ifelse(seq_source != "Not found", 
                           paste0("seq:", seq, "; p1:", primer1, "; p2:", primer2),
                           ""),
         citation = "Hunt et al. 1998; Hunt et al. 2007",
         positions_from = ifelse(seq_source == "Not found",
                                 "No sequences found.",
                                 paste0("Sequences from ",
                                 ifelse(seq_source=="", 
                                        primer_seq_source,
                                        ifelse(primer_seq_source=="", 
                                               seq_source,
                                               paste(primer_seq_source, "and", seq_source))),
                                 " mapped to HAv3.1 with BLASTn")),
         positions_from = str_replace(positions_from, "Harpur", "Harpur 2020 (Dryad)"),
         positions_from = str_replace(positions_from, "Emore 2003", "Emore 2003 (main text)"),
         positions_from = str_replace(positions_from, "1998", "1998 (main text)")) %>%
  rename(QTL_name = qtl) %>%
  mutate(alternative_notation = ifelse(alternative_notation == "GRP7.17-3'", "GRP7.17-3", alternative_notation)) %>% # prime messes up reading in values to table
  dplyr::select(phenotype, QTL_name, marker, comments, chr_group, scaffold, start, end, citation, positions_from, alternative_notation, sequence)
common_chr = hunt0 %>% 
  group_by(QTL_name) %>% 
  slice(which.max(table(chr_group))) %>% dplyr::select(QTL_name, chr_group) %>%
  rename(common_chr = chr_group)
hunt <- hunt0 %>%
  left_join(., common_chr, by = "QTL_name") %>%
  mutate(comments_extra = ifelse(!is.na(chr_group) & chr_group != common_chr,
                          "; (! different chromosome)",
                           ""),
         comments = paste0(comments, comments_extra)) %>%
  dplyr::select(-c(common_chr, comments_extra))

qtl_all <- bind_rows(hunt, tsuruda, qtl_spoetter, oxley, av) %>%
  dplyr::select(colnames(hunt)) %>%
  dplyr::select(., -sequence) %>%
  rename(LG = "chr_group")

qtl_all[qtl_all == ""] <- NA
#View(qtl_all)
# write out results
write.table(qtl_all, "results/Approximate_QTL_positions_HAv3.1.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")
write.table(qtl_all, "../../bee_manuscript/files_supp/S5_Table_Approximate_QTL_positions_HAv3.1.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")
