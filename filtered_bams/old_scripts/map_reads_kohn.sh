#!/bin/bash

# A script to do raw read mapping/alignment to the honeybee reference genome
# e.g. to run one sample: Mexico_Honeybee_001_S1_L008_R1_001.fastq.gz
# scripts$ ./map_reads_kohn.sh SanDiego001 L007 San_Diego_Honeybee_001

# general bash script settings to make sure if any errors in the pipeline fail
# the it's a 'fail' and it passes all errors to exit and allows no unset variables
# -x is for verbose output of commands as they are run
set -o pipefail
set -o errexit
set -o nounset
#set -x

# command line arguments:
ID=$1 # bee ID, e.g. SanDiego001
SEQ_RUN=$2 # e.g. L007
PREFIX=$3 # e.g. San_Diego_Honeybee_001

# make output directory in filtered_bams, one per lane of sequencing
mkdir -p $SEQ_RUN

echo ${PWD} # print current working directory

echo "mapping reads with bowtie2"
bowtie2 --seed 2014 --very-sensitive-local --local \
-x ../data/honeybee_genome/honeybee_Amel_4.5 \
-1 ../data/kohn_data/${PREFIX}*${SEQ_RUN}_R1_001.fastq.gz \
-2 ../data/kohn_data/${PREFIX}*${SEQ_RUN}_R2_001.fastq.gz \
--rg-id ${SEQ_RUN} --rg SM:${ID} | \
samtools view -b - > $SEQ_RUN/${ID}.bam

echo "alignment done!"

# options:
# choosing very sensitive local alignment (i.e. reads may be soft-clipped to map)
# random seed is used to choose a mapping location when there are ties
# -x is for indexed reference genome
# -1 and -2 are for paired reads; -U would be used for unpaired reads
# converting sam directly to bam using samtools
# --rg-id and --rg create read group headers required by GATK downstream analysis
