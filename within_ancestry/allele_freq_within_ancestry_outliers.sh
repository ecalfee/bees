#!/bin/bash -l

# this script calculates allele frequencies within ancestry for a population and list of bams


# to run: ./allele_freq_within_ancestry.sh outlier1_chr1 NA.C NC_037638.1:5000000-15000000 combined_sept19 Group1

OUTLIER_NAME="${1}"
POP_ANCESTRY="${2}"
REGION="${3}"
ANCESTRY_HMM_PREFIX="${4}"
SNP_PREFIX="${5}"

BAM_LIST="../bee_samples_listed/byPop/${ANCESTRY_HMM_PREFIX}_pops.${POP_ANCESTRY}.bams"
DIR_OUT="results/combined_sept19/outliers/${OUTLIER_NAME}"
FILE_OUT="${DIR_OUT}/${POP_ANCESTRY}"
SNP_FILE="../geno_lik_and_SNPs/results/${ANCESTRY_HMM_PREFIX}/variant_sites/${SNP_PREFIX}.var.sites"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_HAv3.1.fasta"

echo "Calculating allele frequencies for OUTLIER: $OUTLIER_NAME $POP_ANCESTRY $REGION"

# make directory to store output (if doesn't yet exist)

mkdir -p "$DIR_OUT"

# bams to include
echo "bam list:" $BAM_LIST

echo "finding site allele frequencies"
angsd -out "${FILE_OUT}" \
-bam "$BAM_LIST" \
-ref "$REF" \
-underFlowProtect 1 \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-r "$REGION" \
-sites "${SNP_FILE}" \
-doCounts 1 \
-doMaf 8 \
-P 4

echo "done getting allele frequencies!"

# options
# basic quality filtering for reads
# -doMajorMinor 3: takes major & minor allele from sites file
# ANGSD calculates freq of minor allele (ignoring all other non maj-min alleles seen)
# -doMaf 8 is an unbiased estimator of pop freq using straight counts and no HWE assumption; -doMaf 8 with -doCounts 1
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)

# underFlowProtect is necessary for large #s of bams
