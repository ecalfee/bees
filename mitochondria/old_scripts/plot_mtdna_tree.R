# script for creating a haplotype tree and plotting it for my mtDNA haplotypes
library(ape)

# read in sample IDs:
ids <- read.table("ARCA2018_ref.list", header = F, stringsAsFactors = F)$V1

# read in metadata:
meta <- left_join(data.frame(Bee_ID = ids, stringsAsFactors = F), 
                  read.table("../bee_samples_listed/all.meta", sep = "\t", header = T, stringsAsFactors = F),
                  by = "Bee_ID")

# read in haplotypes from ANGSD
haplos <- read.table("results/haplotypes/ARCA2018_ref.haplo.gz", stringsAsFactors = F, header = T, na.strings = "N")

# look at data and filter for quality of SNPs - how many are tri-allelic sites?
sum(apply(haplos, 1, function(x) sum(!is.na(unique(x)))) >4)/nrow(haplos) # about 6%, filter out:
haplos_no3 <- haplos[apply(haplos, 1, function(x) sum(!is.na(unique(x)))) <= 4, ] %>%
  dplyr::select(., - c("chr", "pos", "major")) # keep only genotype columns
colnames(haplos_no3) <- ids
haplos_no3.snps <- haplos[apply(haplos, 1, function(x) sum(!is.na(unique(x)))) <= 4, ] %>%
  dplyr::select(., c("chr", "pos", "major")) # keep only snp position columns
# how much information do I have about haplotypes?
summary(apply(haplos_no3, 2, function(x) sum(!is.na(x)))) # not that much missing data


# convert to ape recognized format:
a <- as.DNAbin(t(as.matrix(haplos_no3)), fill.with.gaps = T)
print(a) 

# alternatively, create a binary matrix where rows are loci and columns are individuals. 
# 1 means major allele 0 means minor allele
a.01 <- haplos_no3
for (i in 1:nrow(haplos_no3)){
  for (j in 1:ncol(haplos_no3)){
    a.01[i, j] <- as.integer(haplos_no3[i, j] != haplos_no3.snps[i, "major"]) # 0 for major allele, 1 for minor, NA for missing data
  }
}

# use NJ tree from ape package
# get distance under a particular model of molecular evolution
# Q: is it ok that I only include variant sites? I may need the whole fasta sequence or my tree is just wonky because of low coverage data
a.dist <- dist.dna(a, model = "K80", variance = FALSE,
         gamma = FALSE, pairwise.deletion = FALSE,
         base.freq = NULL, as.matrix = FALSE)
# plot neighbor-joining tree
a.nj <- nj(a.dist)
plot(a.nj, tip.color = rainbow(6)[as.factor(meta$group)])
print(a.nj)


