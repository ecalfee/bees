# chr data
chr_lengths <- read.table("../data/honeybee_genome/chr.list", sep = "\t", stringsAsFactors = F)
colnames(chr_lengths) <- c("scaffold", "length", "chr_group", "chr_lg")

# spoetter qtls
spoetter <- read.table("../data/Spoetter_2012_QTLs/spoetter_2012_Amel3_start_end.csv",
                       sep = ",", header = T, stringsAsFactors = F) %>%
  left_join(., rename(chr_lengths, chr = chr_group), by = "chr") %>%
  mutate(QTL = 1:nrow(.)) %>%
  dplyr::select(scaffold, start, end, QTL, chr_lg) %>%
  rename(chr = scaffold, LG = chr_lg)
# 44k snp chip
snps_44k <- read.table("../data/Spoetter_2012_QTLs/SNP_Info_HoneyBeeSNPAssay.csv",
                       sep = ",", header = T, stringsAsFactors = F)
head(snps_44k)
spoetter$nearest_SNP <- unlist(lapply(1:nrow(spoetter), function(i) filter(snps_44k, 
                                 linkage.group == spoetter$chr[i]) %>%
              summarise(best = Cooridnate[which.min(abs(spoetter$pos[i] - Cooridnate))])))
spoetter_snps <- left_join(spoetter, snps_44k[ , c("linkage.group", "Cooridnate", "left.flanking.sequence", "right.flanking.sequence")],
                           by = c("chr"="linkage.group", "nearest_SNP"="Cooridnate"))
spoetter_snps %>%
  dplyr::select(QTL_study, QTL_name, marker_type, chr, pos, source, nearest_SNP, left.flanking.sequence, right.flanking.sequence, chr_HAv3.1, pos_HAv3.1_start, pos_HAv3.1_end) %>%
  write.table(., "results/Spoetter_2012_QTL_coordinates.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)
snps_44k %>%
  arrange(linkage.group, Cooridnate) %>%
  write.table("../data/Spoetter_2012_QTLs/SNPs_44k.txt",
              quote = F, sep = "\t", col.names = T, row.names = F)

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
  filter(!is.na(start)) %>% # filter out snps with no match
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
qtl_spoetter <- do.call(rbind,
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

table(paste0(qtl_spoetter$qtl, qtl_spoetter$scaffold_HAv3.1))

