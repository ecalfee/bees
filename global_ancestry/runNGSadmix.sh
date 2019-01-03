#!/bin/bash

# slurm array task id sets number of genetic clusters, e.g.
# set an --array=2 for K = 2 or --array=2-4 to test K = 2, 3, 4 etc.

# set VARIABLES
# The first argument K is the number of clusters
# and the second argument is the prefix for the genotype likelihood input file (also used to name output file)
k=$1
GL_FILE="results/input/$2.beagle.gz"
DIR_OUT="results/NGSAdmix"
FILE_OUT="${DIR_OUT}/K${k}_$2"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

echo "running NGSadmix"
NGSadmix -likes "${GL_FILE}" \
-K "${k}" -P 1 \
-o ${FILE_OUT}

# settings:
# -likes beagle genotype likelihood file
# -K 2 for number of subpopulations/clusters to consider in admixture model
# -P n splits the analysis job over n nodes (but does not distribute I/O)
# -o output
