#!/bin/bash

# A script to sort reads, before calculating base adjusted quality scores (BAQ)
# and filter out PCR duplicates (while saving a stats log for % duplication)

# general bash script settings to make sure if any errors in the pipeline fail
# the it's a 'fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# command line arguments:
# note: all paths relative to bees/filtered_bams
ID=$1 # bee ID
DIR_IN=$2 # directory where bam is stored (pre-filtering)..local path from bees/filtered_bams
BAM_IN=${DIR_IN}"/"${ID}".bam" # full path to starting bam
DIR_OUT=${DIR_IN} # results directory
BAM_SORTED=${DIR_OUT}"/"${ID}".sort.bam"
DIR_TMP=${DIR_OUT}"/tmp/"${ID}

# make results and temporary directory (if necessary)
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_TMP}
echo ${PWD} # print current working directory

echo "sorting reads with samtools"
# (1) SAMTOOLS sort reads by coordinate position
samtools sort -m 6G -T ${DIR_TMP} \
-o ${BAM_SORTED} \
${BAM_IN}

echo "all sorted!"
