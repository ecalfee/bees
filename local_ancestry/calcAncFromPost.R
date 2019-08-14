# helper script calculates ancestry from posterior
#!/usr/bin/env Rscript

# this script takes in an admixed population #
# and a directory where to find the ID.posterior files from ancestry_hmm
# and outputs 3 types files into a new 'anc' subdirectory:
# a ID.A.anc ID.C.anc ID.M.anc files with ancestry for all individuals in a population individually
# a POP.A.anc files with mean ancestry for all individuals in a population POP

library(dplyr)

# to run:
# Rscript calcAncFromPost.R AR01 results/ancestry_hmm/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
POP = args[1]
# paths to input and output directories
dir_input = args[2]
# individuals file with IDs of samples with local ancestry calls:
ids_file = paste0("../bee_samples_listed/byPop/", POP, ".list")
dir_output = paste0(dir_input, "/anc")
if (file.exists(dir_input) & !file.exists(dir_output)){ # make sure output directory exists first (or create)!
  dir.create(file.path(dir_output), recursive = T)
}

# find included individuals from current pop
pop_ids <- read.table(ids_file, stringsAsFactors = F,
                      header = F, sep = "\t")$V1

# genotypes
genotypes = c("CC", "CM", "AC", "MM", "AM", "AA") # used in calc_anc_from_post
CMA <- c("C", "M", "A")

# get ancestry for all individuals
# rows = SNPs; columns = individuals
# helper function reads posterior files
read_post <- function(id, dir = dir_post){
  post <- read.table(paste0(dir, "/", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., CC = X2.0.0) %>%
    rename(., CM = X1.1.0) %>%
    rename(., AC = X1.0.1) %>%
    rename(., MM = X0.2.0) %>%
    rename(., AM = X0.1.1) %>%
    rename(., AA = X0.0.2) %>%
    mutate(., Bee_ID = id)
}

# ancestry proportions
calc_anc_from_post <- function(post, genos = genotypes){
  anc1 <- c(1, .5, .5, 0, 0, 0)
  anc2 <- c(0, .5, 0, 1, .5, 0)
  anc3 <- c(0, 0, .5, 0, .5, 1)
  anc_matrix <- cbind(anc1, anc2, anc3)
  anc <- as.data.frame(t(apply(post[ , genos], 1, function(row) row %*% anc_matrix)), stringsAsFactors = F)
  colnames(anc) <- c("C", "M", "A")
  return(anc)
}

# get ancestry proportions across the genome
anc <- lapply(pop_ids, function(id)
  calc_anc_from_post(post = read_post(id = id, dir = dir_input)))

# helper function to write individual ancestry files
write_anc <- function(anc, label){
  lapply(CMA, function(x) 
    write.table(anc[ , x],
                paste0(dir_output, "/", label, ".", x, ".anc"),
                col.names = F, row.names = F, quote = F, sep = "\t"))
}

# write individual ancestries
lapply(1:length(pop_ids), function(i) write_anc(anc = anc[[i]], label = pop_ids[i]))

# calculate and write out mean individual genomewide means for C M and A ancestries, alpha
a <- data.frame(ID = pop_ids, 
                      do.call(rbind,
             lapply(anc, function(x) apply(x, 2, mean))))
write.table(a, paste0(dir_output, "/", POP, ".alpha.anc"), quote = F, col.names = T, row.names = F)

# calculate mean population ancestry across individuals w/in a pop
pop_anc <- anc[[1]]
for (i in 2:length(anc)){
  pop_anc <- pop_anc + anc[[i]]
}
pop_anc <- pop_anc/length(anc)
#pop_anc <- Reduce("+", anc)/length(pop_ids) # didn't work

# write population ancestries
write_anc(anc = pop_anc, label = POP)
