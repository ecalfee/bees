#!/bin/bash

# concatenates var.sites files to match nth GL file, only taking 1 in N positions
# for set of scaffolds specified in SCAFFOLDS_FILE
N=$1
PREFIX=$2
DIR_IN="../geno_lik_and_SNPs/results/${PREFIX}/variant_sites"
WHICH_SCAFFOLDS=$3 # with names of scaffolds to include, e.g. chr or GroupUN or scaffolds
SCAFFOLDS_FILE="../data/honeybee_genome/${WHICH_SCAFFOLDS}.names"
DIR_OUT="results/input"
FILE_OUT="${DIR_OUT}/${PREFIX}_${WHICH_SCAFFOLDS}_prunedBy${N}.var.sites"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
mkdir -p $DIR_OUT

for i in $(cat $SCAFFOLDS_FILE); do cat "$DIR_IN/$i".var.sites; done | awk -v N=$N 'NR % N == 0' > "$FILE_OUT"

echo "done"
