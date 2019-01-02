#!/bin/bash

# A script to add readgroups to bams using picardtools
# this is only relevant for bams that I did not map myself with bowtie2
# because read groups can also be added when mapping.

# general bash script settings to make sure if any errors in the pipeline fail
# the it is a 'fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# directory to PICARD tools:
PICARD="/home/erin/Software/picard-2.18.2"

# command line arguments:
# note: all paths relative to bees/filtered_bams/
ID=$1 # bee ID
RG_ID=$2 # readgroup ID / library name
DIR_IN=$3 # directory for input
DIR_OUT=${DIR_IN}"/RG_added" # directory for output
BAM_IN=${DIR_IN}"/"${ID}".sort.dedup.baq.bam" # full path to starting bam
BAM_OUT=${DIR_OUT}"/"${ID}".sort.dedup.baq.bam"

# make results directory (if necessary)
mkdir -p ${DIR_OUT}

echo "bee: "${ID}
echo "input bam: "${BAM_IN}
echo "working directory:"${PWD} # print current working directory
echo "picard path:" ${PICARD}

echo "adding readgroup IDs with PICARD"
java -jar ~/Software/picard/picard.jar AddOrReplaceReadGroups INPUT=${BAM_IN} \
OUTPUT=${BAM_OUT} SORT_ORDER=coordinate RGID=${RG_ID} RGLB=${RG_ID} RGPL=illumina \
RGPU=NA RGSM=${ID}

echo "indexing with samtools"
sleep 5s # because index needs to have a later timestamp
samtools index ${BAM_OUT}

echo "all done!"
