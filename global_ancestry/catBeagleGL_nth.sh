#!/bin/bash

# concatenates beagle genotype likelihood files, skipping headers and only taking 1 in N positions
# for set of scaffolds specified in SCAFFOLDS_FILE
N=$1
PREFIX=$2
DIR_GL_IN="../geno_lik_and_SNPs/results/${PREFIX}/GL_by_scaffold"
WHICH_SCAFFOLDS=$3 # with names of scaffolds to include, e.g. chr or GroupUN or scaffolds
SCAFFOLDS_FILE="../data/honeybee_genome/${WHICH_SCAFFOLDS}.names"
DIR_OUT="results/input"
FILE_OUT="${DIR_OUT}/${PREFIX}_${WHICH_SCAFFOLDS}_prunedBy${N}.beagle.gz"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
mkdir -p $DIR_OUT

# concatenate header from first file with every Nth SNP from each subsequent file
zcat "$DIR_GL_IN/"$(head -n 1 $SCAFFOLDS_FILE).beagle.gz | head -n 1 | gzip > "$FILE_OUT"; \ # header
for i in $(cat $SCAFFOLDS_FILE); do zcat "$DIR_GL_IN/$i".beagle.gz | tail -n +2; done | awk -v N=$N 'NR % N == 0' | gzip >> "$FILE_OUT"

echo "done running with cat beagle.gz every $N th position for scaffolds in $SCAFFOLDS_FILE, outputting to $FILE_OUT"
