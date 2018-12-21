#!/bin/bash

# this script generalizes input for fastq_prefix so it works with the downloaded ncbi data too

# A script to do raw read mapping/alignment to the honeybee reference genome
# e.g. to run one sample: scripts$ ./map_reads.sh AR0302 C202SC18101772
# or using nohup and gnu parallel to run many samples:
#scripts$ nohup parallel --noswap --joblog ../logs/map_reads_Harpur_2014_NCBI.log --jobs 4 \
#./map_reads_ncbi.sh {1} Harpur_2014_NCBI :::: ../data/Harpur_2014_NCBI/samples.list \
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
#FASTQ_PREFIX="novo_seq/"${2}"/raw_data/"${ID}"/"${ID}"_"
FASTQ_PREFIX=$3
OUTPUT_PREFIX=$4

# move from scripts to data directory
cd ../data

mkdir -p novo_seq/bam_files
echo ${PWD} # print current working directory

echo "mapping reads with bowtie2"

bowtie2 --seed 2014 --very-sensitive-local --local \
-x honeybee_genome/honeybee_Amel_4.5 \
-1 ${FASTQ_PREFIX}*_1.fq.gz \
-2 ${FASTQ_PREFIX}*_2.fq.gz \
--rg-id ${SEQ_RUN} --rg SM:${ID} | \
samtools view -b - > novo_seq/bam_files/${ID}.bam

echo "alignment done!"

# options:
# choosing very sensitive local alignment (i.e. reads may be soft-clipped to map)
# random seed is used to choose a mapping location when there are ties
# -x is for indexed reference genome
# -1 and -2 are for paired reads; -U would be used for unpaired reads
# converting sam directly to bam using samtools
# --rg-id and --rg create read group headers required by GATK downstream analysis
