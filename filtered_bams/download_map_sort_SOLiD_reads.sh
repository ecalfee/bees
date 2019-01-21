#!/bin/bash

# A script to download SOLiD data from NCBI, map using SHRiMP, sort by coordinate,
# and output a sorted bam file
# to run:
#filtered_bams$ ./download_map_sort_SOLiD_reads.sh SRR1151485 SRS549155

# general bash script settings to make sure if any errors in the pipeline fail
# the it's a 'fail' and it passes all errors to exit and allows no unset variables
set -o pipefail
set -o errexit
set -o nounset

# command line arguments:
# note: all paths relative to bees/filtered_bams
SRA=$1 # unique NCBI SRA accession number (could be many per sample due to diff. lanes)
ID=$2 # unique bee sample ID
DIR_OUT="SOLiD/"${ID} # results directory
DIR_TMP="/tmp/"${SRA}
REF="../data/honeybee_genome/Amel_4.5_scaffolds.fa"

# make results and temporary directory (if necessary)
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_TMP}
echo ${PWD} # print current working directory

# download data from NCBI using sra-toolkit's fastq-dump
echo "downloading data from NCBI, SRA: "${SRA}
fastq-dump --dumpcs --skip-technical --split-files -O ${DIR_OUT} "${SRA}"

# reformat using cutadapt so colorspace fastq can be read by the mapping program SHRiMP
echo "re-formating fastq file with cutadapt"
cutadapt -c --format=sra-fastq ${DIR_OUT}/${SRA}_1.fastq > ${DIR_OUT}/${SRA}.fastq

echo "delete intermediate file here: "${DIR_OUT}"/"${SRA}"_1.fastq"

echo "map reads with SHRiMP"
/Users/ecalfee/Software/SHRiMP_2_2_2/bin/gmapper-cs --sam --fastq --single-best-mapping \
${DIR_OUT}/${SRA}.fastq --read-group ${SRA},${ID} ${REF} | \
samtools view -b - > ${DIR_OUT}/${SRA}.bam

echo "delete intermediate file here: "${DIR_OUT}"/"${SRA}".fastq"

# SAMTOOLS sort reads by coordinate position
echo "sorting reads with samtools"
samtools sort -m 6G -T ${DIR_TMP} \
-o ${DIR_OUT}/${SRA}.sort.bam \
${DIR_OUT}/${SRA}.bam

echo "delete intermediate file here: "${DIR_OUT}"/"${SRA}".bam"

echo "all done!"
