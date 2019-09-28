#!/bin/bash


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables

# this script takes in a bed file for the ancestry caller's tracts
# and outputs a bed file for an individual sample where that sample has 
# homozygous ancestry calls above some threshold (0.8) for ancestry ANCESTRY, e.g. C
# to run: ./getBamsForTracts.sh CA1207 3 C combined_sept19

BEE_ID="${1}"
ANCESTRY="${2}"
ANCESTRY_HMM_PREFIX="${3}"
TRACTS_FILE=../local_ancestry/results/tracts/"${ANCESTRY_HMM_PREFIX}"/"${ANCESTRY}"/"$BEE_ID".bed
DIR_OUT=results/"$ANCESTRY_HMM_PREFIX"/"$ANCESTRY"/bams

DIR_BAMS="../filtered_bams/results"
BAM=$DIR_BAMS/$BEE_ID.sort.dedup.baq.bam
# genome order for HAv3.1. must list all scaffolds in original bam file.
GENOME_ORDER="../data/honeybee_genome/scaffolds_with_mtDNA.lengths"

# make output and scratch directories
mkdir -p "$DIR_OUT"

echo getting "$ANCESTRY" tracts for bee "$BEE_ID"

# filter bams using bedtools and tracts
echo "filtering bam"
bedtools intersect -sorted -a "$BAM" \
-b "$TRACTS_FILE" \
-g "$GENOME_ORDER" \
| samtools view -b - > $DIR_OUT/"$BEE_ID".sort.dedup.baq.bam

echo "now indexing new bam!"
sleep 5s # because index needs to have a later timestamp
samtools index $DIR_OUT/"$BEE_ID".sort.dedup.baq.bam

echo "all done!"

# note: I do not impose any additional mapping or other quality filters on these bams, just high conf. ancestry calls.
