#!/bin/bash

echo "argument 1 is: $1"
echo "argument 2 is: $2"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables

# this script takes in a sites file and outputs an angsd counts file with
# columns totA totC totG totT as well as a .pos.gz file with all positions with data
# to run: ./countReadsACGT.sh CA0906 TEST Group1_highLD_CMA

BEE_ID="${1}"
DIR_OUT="results/SNPs/${2}"
SITES_FILE="results/SNPs/${2}/${3}.var.sites"
DIR_IN="../filtered_bams/results/"
BAM="${DIR_IN}/${BEE_ID}.sort.dedup.baq.bam"
DIR_SINGLE_BAM_LIST="../bee_samples_listed/list_indiv_bams"
REF="../data/honeybee_genome/Amel_HAv3.1.fasta"

# make output and scratch directories
mkdir -p "${DIR_OUT}/countsACGT"
mkdir -p ${DIR_SINGLE_BAM_LIST}

# make file with single focus bam listed
echo "${BAM}" > ${DIR_SINGLE_BAM_LIST}/${BEE_ID}.only.list

echo "counting reads supporting ACGT alleles for bee ${BEE_ID}"
# counts reads for sites in sites file that pass basic quality filtering
angsd -out "${DIR_OUT}/countsACGT/${BEE_ID}" \
-ref "$REF" \
-bam "${DIR_SINGLE_BAM_LIST}/${BEE_ID}.only.list" \
-minQ 20 -minMapQ 30 \
-doCounts 1 -dumpCounts 3 \
-remove_bads 1 \
-sites "${SITES_FILE}"

echo "all done!"
