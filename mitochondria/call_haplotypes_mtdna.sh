#!/bin/bash

# this script uses ANGSD to determine mtDNA haplotypes at SNPs

# command line arguments:
# note: all paths relative to bees/mitochrondria/
BAMS_LIST="../bee_samples_listed/combined_sept19.bams" # list of paths to mtDNA bee bams
MAX_MISSING=174 # maximum number of samples with missing data to include a SNP = half of sample
FILE_OUT="results/haplotypes_mtDNA" # name of file for results
SITES="../geno_lik_and_SNPs/results/combined_sept19/variant_sites/mtDNA_highDepthOK.var.sites"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

echo "calling variants and GL using ANGSD on BAMS for mtDNA"

angsd -out "$FILE_OUT" \
-sites "$SITES" \
-r NC_001566.1: \
-ref "../data/honeybee_genome/Amel_HAv3.1.fasta" \
-bam "$BAMS_LIST" \
-doMajorMinor 3 \
-dohaplocall 2 \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doCounts 1 \
-minMinor 2 \
-maxMis "$MAX_MISSING"

echo "all done!"

# settings:
# -remove_bads removes reads with flags like duplicates
# -bam list of bams to include
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# dohaplocall 2 uses the majority base (random base for ties)
# minMinor 2 only includes minor alleles found in >= 2 samples
# maxMis sets the maximum number of individuals with no reads to include a SNP
