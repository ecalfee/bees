#!/bin/bash

# A script to sort and merge bams into one file per sample,
# then filter out PCR duplicates in PICARD

# general bash script settings to make sure if any errors in the pipeline fail
# the it is a 'fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# command line arguments:
# note: all paths relative to bees/filtered_bams/
ID=$1 # bee ID
DIR_TMP="tmp/"${ID} # for memory overflow

# make results directory (if necessary)
mkdir -p ${DIR_TMP}
mkdir -p merged
echo "working directory:"${PWD} # print current working directory

echo "sorting reads with samtools"

for RG in L006 L007 L008
do samtools sort -m 2G -T ${DIR_TMP} \
-o ${RG}/${ID}.sort.bam ${RG}/${ID}.bam
done

echo "all sorted! merging sorted bams"
ls L*/${ID}.sort.bam > merged/${ID}.list # list of bams to merge
# merging
samtools merge -b merged/${ID}.list \
-O merged/${ID}.sort.bam

echo "all done!"
