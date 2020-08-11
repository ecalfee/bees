#!/bin/bash

# this script uses ANGSD to calculate genotype likelihoods (GL)
# at all putative SNPs in the mitochondria
# snp calling for mtDNA is a bit different because
# each individual is haploid and we expect more variable coverage

# command line arguments:
# note: all paths relative to bees/geno_lik_and_SNPs/
BAMS_LIST=$1 # list of paths to bee bams
SCAFFOLD_NAME="mtDNA"
SCAFFOLD_REGION="NC_001566.1" # region to run GL on = a scaffold
MIN_IND=$2 # minimum number of individuals with data for total sample to keep a site
DIR_OUT=$3 # output file goes in this directory, labelled by  scaffold.
INBRED_FILE=$4

# also requires honeybee reference genome (indexed by samtools faidx)

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

echo "calling variants and GL using ANGSD on BAMS for hilo genomic region: "${REGION}

angsd -out ${DIR_OUT}"/"${SCAFFOLD_NAME} \
-r ${SCAFFOLD_REGION}: \
-ref "../data/honeybee_genome/Amel_HAv3.1.fasta" \
-bam ${BAMS_LIST} \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 4 \
-minMaf 0.05 \
-doMaf 2 \
-GL 1 \
-doGlf 2 \
-P 1 \
-minInd ${MIN_IND} \
-indFname ${INBRED_FILE}
#-setMaxDepth ${MAX_DEPTH}

echo "all done!"

# settings:
# -r specifies a genomic region to work on, e.g. scaffoldName: specifies the whole scaffold
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 4: set reference allele as major allele; infers minor from genotype likelihoods
# -doMaf 1: use the inferred major/minor allele to calculate MAF based on genotype likelihoods
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x minor allele freq.
# -minInd N: only keep sites with information (at least one read) from N individuals (skipped for now)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
# -setMaxDepth -setMinDepth filters out sites where total depth is below or exceeds some threshold
