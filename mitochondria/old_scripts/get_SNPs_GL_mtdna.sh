#!/bin/bash

# this script uses ANGSD to calculate genotype likelihoods (GL)
# at all putative SNPs for a region of the genome

# command line arguments:
# note: all paths relative to bees/geno_lik_and_SNPs/
BAMS_LIST=$1 # list of paths to bee bams
SCAFFOLD=$2 # region to run GL on = a scaffold
MIN_DEPTH=$3 # minimum depth for total sample to keep a site
MAX_DEPTH=$4 # maximum depth for total sample to keep a site
DIR_OUT=$5 # output file goes in this directory, labelled by  scaffold.

# also requires honeybee reference genome (indexed by samtools faidx)

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

echo "calling variants and GL using ANGSD on BAMS for hilo genomic region: "${REGION}

angsd -out ${DIR_OUT}"/scaffold_"${SCAFFOLD} \
-r ${SCAFFOLD}: \
-ref "../data/honeybee_genome/Amel_4.5_scaffolds.fa" \
-bam ${BAMS_LIST} \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 2 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-GL 1 -doGlf 2 \
-P 1 \
-setMaxDepth ${MAX_DEPTH} \
-setMinDepth ${MIN_DEPTH}

echo "all done!"

# settings:
# -r specifies a genomic region to work on, e.g. scaffoldName: specifies the whole scaffold
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 2: infer major and minor from allele counts
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# we use AlleleCounts method for MAF (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x minor allele freq.
# -minInd N: only keep sites with information (at least one read) from N individuals (skipped for now)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
# -setMaxDepth -setMinDepth filters out sites where total depth is below or exceeds some threshold
