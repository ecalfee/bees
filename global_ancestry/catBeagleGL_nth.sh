#!/bin/bash

# concatenates beagle genotype likelihood files, skipping headers and only taking 1 in N positions
# can change which subset is selected by choosing a different order of scaffolds/regions in input file SCAFFOLDS_LIST
N=$1
GL_FILE_LIST=$2
DIR_GL_IN=$3
OUTPUT_NAME=$4
DIR_OUT="results/input"
FILE_OUT="${DIR_OUT}/${OUTPUT_NAME}_prunedBy${N}.beagle.gz"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
mkdir -p $DIR_OUT

HEADER=$(zcat "$DIR_GL_IN/"$(head -n 1 $GL_FILE_LIST) | head -n 1)

# concatenate header from first file with every Nth SNP from each subsequent file
zcat $DIR_GL_IN"/"$(head -n 1 $GL_FILE_LIST) | head -n 1 | gzip > "$FILE_OUT"; \
for i in $(cat $GL_FILE_LIST); do zcat "$DIR_GL_IN/$i" | tail -n +2; done | awk -v N=$N 'NR % N == 0' | gzip >> $FILE_OUT

echo "done running with cat beagle.gz every $N th position for scaffolds in $GL_FILE_LIST, outputting to $FILE_OUT"
