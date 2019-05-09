#!/bin/bash

# this script maps reads that are single reads (no pairs), e.g. reference bees downloaded from NCBI
# to a mtDNA reference genome

# A script to do raw read mapping/alignment to the honeybee reference genome
# e.g. to run one sample: scripts$ ./map_reads_single.sh SRR957082 Harpur_2014_NCBI
# or using nohup and gnu parallel to run many samples:
#mitochondria$ nohup parallel --noswap --joblog logs/map_mtdna_Harpur_2014_NCBI.log --jobs 4 \
#./map_reads_single.sh {1} Harpur_2014_NCBI :::: ../data/Harpur_2014_NCBI/samples.list \
# &> logs/map_mtdna_Harpur_2014_NCBI.out &

# note: the fastq files from NCBI end with fastq.gz not fq.gz

# general bash script settings to make sure if any errors in the pipeline fail
# the it is a "fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# command line arguments:
ID=$1 # bee ID
SEQ_RUN=$2 # identifier for sequenced group
FASTQ_PREFIX=../data/${SEQ_RUN}"/fastq_files/"${1}
DIR_OUT="results/bam_files"
REF="../data/honeybee_genome/mitochondria/KY926884"

mkdir -p ${DIR_OUT}


echo "mapping reads with bowtie2"

bowtie2 --seed 2014 --very-sensitive-local --local \
-x honeybee_genome/honeybee_Amel_4.5 \
-U ${FASTQ_PREFIX}_1.fastq.gz \
--rg-id ${SEQ_RUN} --rg SM:${ID} | \
samtools view -q 1 -b - > ${DIR_OUT}/${ID}.bam

echo "alignment done!"

# options:
# choosing very sensitive local alignment (i.e. reads may be soft-clipped to map)
# random seed is used to choose a mapping location when there are ties
# -x is for indexed reference genome
# -U is for unpaired reads
# converting sam directly to bam using samtools
# --rg-id and --rg create read group headers required by GATK downstream analysis
# samtools -q 1 skips reads with 0 mapping quality (i.e. did not map)
