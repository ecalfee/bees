#!/bin/bash -l

# this script calculates allele frequencies for a population and list of bams

# to run: ./allele_freq.sh AR01 Group1 combined_sept19

POP="${1}"
SCAFFOLD="${2}"
PREFIX="${3}"
BEE_ID_FILE="../bee_samples_listed/byPop/${POP}.list"
SNP_FILE="../geno_lik_and_SNPs/results/${PREFIX}/variant_sites/${SCAFFOLD}.var.sites"
REGIONS_FILE="../geno_lik_and_SNPs/results/${PREFIX}/variant_sites/${SCAFFOLD}.regions"
DIR_OUT="results/${PREFIX}/allele_freq/combined/${SCAFFOLD}"
BAM_LIST="../bee_samples_listed/byPop/${POP}.bams"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_HAv3.1.fasta"

echo "Calculating allele frequencies for POP: "$POP" across all ancestries"

# make directory to store output (if doesn't yet exist)

mkdir -p "$DIR_OUT"

echo "bam list:" $BAM_LIST

echo "finding site allele frequencies"
angsd -out "${DIR_OUT}/${POP}" \
-bam "$BAM_LIST" \
-ref "$REF" \
-rf "$REGIONS_FILE" \
-underFlowProtect 1 \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites "${SNP_FILE}" \
-doCounts 1 \
-doMaf 8 \
-P 1

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
