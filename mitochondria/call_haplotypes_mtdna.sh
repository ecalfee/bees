#!/bin/bash

# this script uses ANGSD to determine mtDNA haplotypes

# command line arguments:
# note: all paths relative to bees/mitochrondria/
BAMS_LIST=$1 # list of paths to mtDNA bee bams
MAX_MISSING=$2 # maximum number of samples with missing data to include a SNP
DIR_OUT="results/haplotypes"
NAME_OUT=$3 # name of file for results
REF="../data/honeybee_genome/mitochondria/KY926884.fasta"

# also requires honeybee reference genome (indexed by samtools faidx)

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

echo "calling variants and GL using ANGSD on BAMS for mtDNA"

angsd -out "$DIR_OUT/$NAME_OUT" \
-ref "$REF" \
-bam "$BAMS_LIST" \
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
