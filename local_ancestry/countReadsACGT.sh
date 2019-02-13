#!/bin/bash -l

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# this script takes in a sites file and outputs an angsd counts file with
# columns totA totC totG totT as well as a .pos.gz file with all positions with data
# to run: ./countReadsACGT CA0906

ID=$1
DIR_OUT="results/SNPs/thin1kb_common3"
DIR_IN="../filtered_bams/results/"
BAM="${DIR_IN}/${ID}.sort.dedup.baq.bam"
DIR_SINGLE_BAM_LIST="../bee_samples_listed/list_indiv_bams"
REF="../data/honeybee_genome/honeybee_Amel_4.5"

# make output and scratch directories
mkdir -p "${DIR_OUT}/countsACGT"
mkdir -p ${DIR_SINGLE_BAM_LIST}

# make file with single focus bam listed
echo ${BAM} > ${DIR_SINGLE_BAM_LIST}/${ID}.only.list

echo "counting reads supporting ACGT alleles for bee ${ID}"
# counts reads for sites in sites file that pass basic quality filtering
angsd -out "${DIR_OUT}/countsACGT/${ID}" \
-ref "$REF" \
-bam "${DIR_SINGLE_BAM_LIST}/${ID}.only.list" \
-minQ 20 -minMapQ 30 \
-doCounts 1 -dumpCounts 3 \
-remove_bads 1 \
-sites "${DIR_OUT}/included.var.sites"

echo "all done!"
