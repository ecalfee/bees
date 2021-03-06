#!/bin/bash


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables

# this script takes in a bed file for the ancestry caller's tracts
# and outputs a bed file for an individual sample where that sample has 
# homozygous ancestry calls above some threshold (0.8) for ancestry ANCESTRY, e.g. A
# to run: ./getHighPosteriorTracts.sh CA1207 3 C combined_sept19 results/SNPs/combined_sept19/chr.var.sites.bed

BEE_ID="${1}"
POSTERIOR_COLUMN="${2}"
ANCESTRY="${3}"
ANCESTRY_HMM_PREFIX="${4}"
SNP_BED="${5}"
DIR_OUT=results/tracts/"${ANCESTRY_HMM_PREFIX}"/"${ANCESTRY}"
# sort to genome order for Amel_HAv3.1
GENOME_INDEX="../data/honeybee_genome/Amel_HAv3.1.fasta.fai"

# make output and scratch directories
mkdir -p "$DIR_OUT"

echo getting "$ANCESTRY" tracts for bee "$BEE_ID"

# take posterior calls over threshold for ancestry specified
# then merge adjacent tracts using bedtools
pr -mt -s$'\t' "$SNP_BED"  <(tail -n +2 results/ancestry_hmm/"$ANCESTRY_HMM_PREFIX"/posterior/"$BEE_ID".posterior | cut -f"$POSTERIOR_COLUMN") | awk '$5 > 0.8 {print $0}' | bedtools merge -sorted -d 1 | bedtools sort -faidx "$GENOME_INDEX" -i stdin > "$DIR_OUT"/"$BEE_ID".bed
