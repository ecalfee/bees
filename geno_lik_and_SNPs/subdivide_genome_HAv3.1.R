# take in fasta.fai file and information file from NCBI
# print out subsets of the genome with names as .lengths files
library(dplyr)
library(bedr)

fai <- read.table("../data/honeybee_genome/Amel_HAv3.1.fasta.fai",
                  header = F, sep = "\t")[ , 1:2] %>%
  data.table::setnames(c("scaffold_id", "length"))
ncbi <- read.table("../data/honeybee_genome/Amel_HAv3.1_NCBI.info",
                   header = F, skip = 33, sep = "\t")[ , c(1,3,7)] %>%
  data.table::setnames(c("name", "LG", "scaffold_id"))
genome <- left_join(fai, ncbi, by = "scaffold_id")
tail(genome)
chromosomes = paste0("Group", 1:16)

# list of spanned gaps in genome
gaps <- read.table("../data/honeybee_genome/allcontig.agp.gz",
                   skip = 3, header = F) %>%
  filter(., V5 == "N") %>% # only look at gaps, not scaffolds
  data.table::setnames(c("scaffold_id", "start_1", "end", "n_on_scaffold", "N",
                         "gap_length", "gap_type", "linkage", "evidence")) %>%
  mutate(start = start_1 - 1) # put on bed 0-based coordinates


# write files
# known gaps
write.table(gaps[ , c("scaffold_id", "start", "end")], 
            "../data/honeybee_genome/gaps.bed",
            col.names = F, row.names = F, quote = F, sep = "\t")

# just scaffolds placed on chromosomes
filter(genome, name %in% chromosomes) %>%
  dplyr::select(scaffold_id, length) %>%
  write.table("../data/honeybee_genome/chr.lengths", 
              col.names = F, row.names = F, quote = F, sep = "\t")
filter(genome, name %in% chromosomes) %>%
  write.table("../data/honeybee_genome/chr.list", 
              col.names = F, row.names = F, quote = F, sep = "\t")
# just mitochondria
filter(genome, name == "MT") %>%
  dplyr::select(scaffold_id, length) %>%
  write.table("../data/honeybee_genome/mtDNA.lengths", 
              col.names = F, row.names = F, quote = F, sep = "\t")
filter(genome, name == "MT") %>%
  write.table("../data/honeybee_genome/mtDNA.list", 
              col.names = F, row.names = F, quote = F, sep = "\t")
# all nuclear scaffolds
filter(genome, name != "MT") %>%
  dplyr::select(scaffold_id, length) %>%
  write.table("../data/honeybee_genome/scaffolds.lengths", 
              col.names = F, row.names = F, quote = F, sep = "\t")
filter(genome, name != "MT") %>%
  write.table("../data/honeybee_genome/scaffolds.list", 
              col.names = F, row.names = F, quote = F, sep = "\t")
# all nuclear scaffolds plus mitochondria
genome %>%
  dplyr::select(scaffold_id, length) %>%
  write.table("../data/honeybee_genome/scaffolds_with_mtDNA.lengths", 
              col.names = F, row.names = F, quote = F, sep = "\t")
genome %>%
  write.table("../data/honeybee_genome/scaffolds_with_mtDNA.list", 
              col.names = F, row.names = F, quote = F, sep = "\t")

# coverage
depth <- read.table("results/combined_sept19/coverage/chr.random_pos_1000.txt",
                    header = F, stringsAsFactors = F)
colnames(depth) <- c("scaffold", "start", "end", ids$Bee_ID)
ind_mean <- apply(depth[, ids$Bee_ID], 2, mean)
site_total <- apply(depth[ , ids$Bee_ID], 1, sum)
n_zeros <- apply(depth[ , ids$Bee_ID], 1, function(x) sum(x == 0))
ids <- read.table("../bee_samples_listed/combined_sept19.list",
                      stringsAsFactors = F) %>%
  data.table::setnames("Bee_ID")
meta <- read.table("../bee_samples_listed/all.meta",
                   stringsAsFactors = F,
                   sep = "\t",
                   header = T)
bees <- left_join(ids, meta, by = "Bee_ID")
summary(site_total)
hist(site_total)
# outliers high coverage
table(site_total > mean(site_total)*2)/length(site_total)
# outliers low coverage
table(site_total < mean(site_total)/2)/length(site_total)
table(site_total < mean(site_total)/2)/length(site_total)
# approximate mean depth is 2748 ~ 2750 * 2 = 5500
# half the individuals is 348/2 = 174 which I'll use as a lower cutoff
# write out estimated individual coverage per bee
write.table(data.frame(Bee_ID = ids$Bee_ID,
                       est_coverage = apply(depth[ , ids$Bee_ID], 2, mean)),
            "results/combined_sept19/coverage/mean_ind_coverage.chr.random_pos_1000.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")

# what is an ok coverage cutoff for the higher coverage bees?
hist(depth[ , "SRR5270368"])
table(depth[ , "SRR5270368"] < 6)
colMeans(depth[ , ids$Bee_ID] < 10)
colMeans(depth[ , ids$Bee_ID] < 7)

