#!/bin/bash -l

# this script calculates genotype likelihoods for all reference samples +
# bams for CA/AR populations filtered to include only regions of high confidence for focal ancestry
# to see how similar populations are within A/C/M ancestry


# to run: ./calc_GL_within_ancestry_4PCA.sh C combined_sept19 Group1 prunedBy250 ../global_ancestry/results/input/combined_sept19_chr_prunedBy250.var.sites
ANCESTRY="${1}"
PREFIX="${2}"
SCAFFOLD="${3}"
OUT="${4}"
SNP_FILE="${5}" # snps
BAM_LIST="results/combined_sept19/${ANCESTRY}/allele_freq/${PREFIX}.bams"
REGIONS_FILE="../geno_lik_and_SNPs/results/${PREFIX}/variant_sites/${SCAFFOLD}.regions"
DIR_OUT="results/${PREFIX}/${ANCESTRY}/GL/${OUT}"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_HAv3.1.fasta"

echo "Calculating genotype likelihoods within ancestry "$ANCESTRY" for scaffold $SCAFFOLD"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

echo "finding within-ancestry genotype likelihoods"
angsd -out "${DIR_OUT}/${SCAFFOLD}" \
-bam "$BAM_LIST" \
-ref "$REF" \
-rf "$REGIONS_FILE" \
-remove_bads 1 \
-minMapQ 30 \
-doMajorMinor 3 \
-minQ 20 \
-GL 1 \
-doGlf 2 \
-sites "${SNP_FILE}" \
-P 1

echo "done getting allele frequencies!"


# options
# basic quality filtering for reads
# -rf limits angsd to walking through indexed regions in the regions file, not the whole genome
# -GL 1 uses samtools GL method
# -doGlf
# -doMajorMinor 3: takes major & minor allele from sites file
# ANGSD calculates freq of minor allele (ignoring all other non maj-min alleles seen)
# -minMapQ 30: filter out sites with low mapping quality.
# -minQ 20 filters for minimum base quality of 20
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -remove_bads removes reads that don't pass quality filters (here these are actually already removed by deduplicating)
