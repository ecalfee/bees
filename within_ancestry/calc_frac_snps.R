# what is the snp density genomewide?
chr_length_tot <- sum(read.table("../data/honeybee_genome/chr.lengths")$V2)
gaps_tot <- read.table("../data/honeybee_genome/gaps.bed") %>%
  mutate(length = V3 - V2) %>%
  summarise(tot = sum(length)) %>%
  unlist(.)
n_snps <- 3510834 # from wc -l chr.var.sites
# genome_size <- chr_length_tot - gaps_tot #236*10^6 # genome size
# actually, use all sites passing quality filters for the denominator of pi:
all_sites <- data.frame(
  bin_name = read.table("../geno_lik_and_SNPs/results/1cM_bins.names", header = F, stringsAsFactors = F)$V1,
  n_sites = read.table("../geno_lik_and_SNPs/results/1cM_bins_total_counts.txt",
                       header = F, sep = "\t")$V2)
genome_size <- sum(all_sites$n_sites)
frac_snps <- n_snps/genome_size # multiplier for any snp-based heterozygosity estimate
# fraction of the genome that is variable, i.e. snps