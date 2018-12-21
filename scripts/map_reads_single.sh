#!/bin/bash

# this script maps reads that are single reads (no pairs), e.g. reference bees downloaded from NCBI

# A script to do raw read mapping/alignment to the honeybee reference genome
# e.g. to run one sample: scripts$ ./map_reads_single.sh SRR957082 Harpur_2014_NCBI
# or using nohup and gnu parallel to run many samples:
#scripts$ nohup parallel --noswap --joblog ../logs/map_reads_Harpur_2014_NCBI.log --jobs 4 \
#./map_reads_single.sh {1} Harpur_2014_NCBI :::: ../data/Harpur_2014_NCBI/samples.list \
# &> ../logs/map_reads_Harpur_2014_NCBI.out &

# general bash script settings to make sure if any errors in the pipeline fail
# the it's a 'fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# command line arguments:
ID=$1 # bee ID
SEQ_RUN=$2 # identifier for sequenced group
FASTQ_PREFIX=${SEQ_RUN}"/fastq_files/"${1}
DIR_OUT=${SEQ_RUN}"/bam_files"

# move from scripts to data directory
cd ../data

mkdir -p ${DIR_OUT}

echo "working directory:"${PWD} # print current working directory
echo "output directory:"${DIR_OUT}

echo "mapping reads with bowtie2"

bowtie2 --seed 2014 --very-sensitive-local --local \
-x honeybee_genome/honeybee_Amel_4.5 \
-U ${FASTQ_PREFIX}*_1.fq.gz \
--rg-id ${SEQ_RUN} --rg SM:${ID} | \
samtools view -b - > ${DIR_OUT}/${ID}.bam

echo "alignment done!"

# options:
# choosing very sensitive local alignment (i.e. reads may be soft-clipped to map)
# random seed is used to choose a mapping location when there are ties
# -x is for indexed reference genome
# -U is for unpaired reads
# converting sam directly to bam using samtools
# --rg-id and --rg create read group headers required by GATK downstream analysis
