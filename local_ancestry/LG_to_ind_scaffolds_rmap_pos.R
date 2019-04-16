# because scaffold orientation changes and some scaffolds
# are not oriented, I want to bring positions back to local scaffold positions
# and just run ancestry_hmm on individual scaffolds
# this script takes LG level positions and rmap and converts to scaffold-level
# eventually I'll need to fix how I convert recombination rate positions in the first place due to backwards oriented scaffolds..
# but this is a quick fix for now:
dir <- "results/SNPs/thin1kb_common3"
sites <- read.table(file.path(dir, "included.var.sites"), stringsAsFactors = F, header = F)
colnames(sites) <- c("scaffold", "pos", "major", "minor")
# get recombination distances between each SNP
# assumes 50kb between scaffolds, but contiguous chromosomes
rdist <- read.table(file.path(dir, "included.rdist"), header = F, stringsAsFactors = F)
rdist_scaffolds <- rdist$V1
# reset the distance at the start of each scaffold to 1
rdist_scaffolds[(cumsum(rle(sites$scaffold)$lengths) + 1)] <- 1 

# write results to file
write.table(sites[, c("scaffold", "pos")], file = file.path(dir, "included_scaffolds.pos"),
            quote = F, row.names = F, col.names = F, sep = "\t")

write.table(rdist_scaffolds, file = file.path(dir, "included_scaffolds.rdist"),
            quote = F, row.names = F, col.names = F, sep = "\t")

