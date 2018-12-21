#!/bin/bash

# A script to filter out PCR duplicates in PICARD
# (while saving a stats log file in filtered_bams/metrics/ for % duplication)
# and re-calibrate base quality scores using SAMTOOLS BAQ
# note: input bams must be sorted
# and $PICARD must be a globally stored variable with path to picard.jar file
# (e.g. add PICARD=global_path_to_picard) to your bash profile ~/.profile)

# general bash script settings to make sure if any errors in the pipeline fail
# the it's a 'fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# move from scripts to data directory
cd ../data

# command line arguments:
# note: all paths relative to ../data/
ID=$1 # bee ID
DIR_IN=$2 # directory for input
BAM_SORTED=${DIR_IN}"/"${ID}".sort.bam" # full path to starting bam
DIR_OUT="filtered_bams/results" # results directory
DIR_METRICS="filtered_bams/metrics" # metrics directory
BAM_OUT=${DIR_OUT}"/"${ID}".sort.dedup.baq.bam"
DIR_TMP=${DIR_OUT}"/tmp/"${ID} # for memory overflow
REF="honeybee_genome/Amel_4.5_scaffolds.fa" # honeybee reference genome (indexed by samtools faidx)

# make results directory (if necessary)
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_TMP}
mkdir -p ${DIR_METRICS}
echo "working directory:"${PWD} # print current working directory

echo "marking duplicates with PICARD and calculating BAQ with SAMTOOLS"
java -Xmx8g -jar ${PICARD}/picard.jar MarkDuplicates \
INPUT=${DIR_OUT}/${ID}.sort.bam OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true \
TMP_DIR=${DIR_TMP} \
METRICS_FILE=${DIR_METRICS}/${ID}.metrics.txt | \
samtools calmd -SbArE --reference ${REF} - > ${BAM_OUT}

echo "all done removing duplicates and calculating BAQ, now indexing!"
sleep 5s # because index needs to have a later timestamp
samtools index ${BAM_OUT}
echo "done indexing"


# options:
# -Xmx8g will not spill to tmp until 8G of memory are used
# samtools calmd calculates adjusted base quality score (BAQ)
# -S specifies that input is sam, not binary bam
# - sets stdin as input (stdout as output is default)
# -Ar compute BAQ and cap base quality with BAQ
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default
