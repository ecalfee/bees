#!/bin/bash -l

# this script calculates allele frequencies within ancestry for a population and list of bams


# to run: ./allele_freq_within_ancestry.sh AR01 AA thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot ../bee_samples_listed/byPop/AR01.list ../data/bees_new_positions/ALL.var.sites

POP="${1}"
ANCESTRY="${2}"
ANCESTRY_HMM_PREFIX="${3}"
BEE_ID_FILE="${4}"
SNP_FILE="${5}"
DIR_BAMS=results/"${ANCESTRY_HMM_PREFIX}"/"${ANCESTRY}"/bams
DIR_OUT=results/"$ANCESTRY_HMM_PREFIX"/"$ANCESTRY"/allele_freq
BAM_LIST="$DIR_OUT"/"$POP".bams

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_4.5_scaffolds.fa"

echo "Calculating allele frequencies for POP: "$POP" and ancestry "$ANCESTRY

# make directory to store output (if doesn't yet exist)

mkdir -p "$DIR_OUT"

# make list of bams to include
for i in $(cat "$BEE_ID_FILE")
	do echo $DIR_BAMS/$i.sort.dedup.baq.bam
done > $BAM_LIST

echo "bam list:" $BAM_LIST

echo "finding site allele frequencies"
angsd -out "${DIR_OUT}/${POP}" \
-bam "$BAM_LIST" \
-ref "$REF" \
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