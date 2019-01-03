#!/bin/bash

# input one variable -- the prefix for the beagle genotype likelihood input file
GL_FILE="results/input/$1.beagle.gz"
DIR_OUT="results/PCA"
FILE_OUT="${DIR_OUT}/$1"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# run PCAngsd
mkdir -p ${DIR_OUT}

# assumes GL data is already filtered for a minimum MAF; doesn't re-filter
python2 ~/Software/pcangsd/pcangsd.py \
-beagle "${GL_FILE}" \
-threads 2 -iter 100 \
-minMaf 0 -admix \
-o "${FILE_OUT}"

# -admix option calculates admixture proportions in addition to genotype covariance matrix for PCA
# -iter specifies number of EM steps

echo "done running PCAngsd"
