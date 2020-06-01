#!/bin/bash

# to run: ./runPCAngsd_within_ancestry.sh M_CA_AR M

# input the prefix for the beagle genotype likelihood input file and the focal ancestry (to find input/output directories)
BEAGLE_PREFIX="$1"
ANCESTRY="$2"
GL_FILE="results/combined_sept19/$ANCESTRY/$GL/${BEAGLE_PREFIX}.beagle.gz"
DIR_OUT="results/combined_sept19/$ANCESTRY/PCA"
FILE_OUT="${DIR_OUT}/${BEAGLE_PREFIX}"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# run PCAngsd
mkdir -p ${DIR_OUT}

python2 ~/Software/pcangsd/pcangsd.py \
-beagle "${GL_FILE}" \
-threads 2 -iter 100 \
-minMaf 0.5 -admix \
-o "${FILE_OUT}"

# -admix option calculates admixture proportions in addition to genotype covariance matrix for PCA
# -iter specifies number of EM steps
# use default

echo "done running PCAngsd within ancestry $ANCESTRY for file $GL_FILE"
