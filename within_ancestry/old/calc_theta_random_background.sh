#!/bin/bash -l

# this script calculates pi for a population and list of bams, e.g. A ref pop
# over a set of random regions specified in the regions file (bed file converted to +1 start position so it's 1 not zero indexed for angsd)
# bed is ok even though it's not 


# to run: ./calc_theta_random_background.sh A ../bee_samples_listed/byPop/A.list 1Mb_of_non-outlier_100bp combined ../filtered_bams/results_Amel4.5 

POP="${1}"
BEE_ID_FILE="${2}"
REGIONS_PREFIX="${3}"
REGIONS_FILE="results/$REGIONS_PREFIX.regions"
ANCESTRY="${4}"
DIR_BAMS="${5}"

DIR_OUT="results/non-outlier_regions/${REGION_PREFIX}/${ANCESTRY}"
BAM_LIST="${DIR_OUT}/${POP}.bams"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_4.5_scaffolds.fa"

echo Calculating pi for POP: "$POP"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

# make list of bams to include
for i in $(cat "$BEE_ID_FILE")
	do echo $DIR_BAMS/$i.sort.dedup.baq.bam
done > $BAM_LIST

echo "finding site allele frequencies"
# steps:
# get saf approximate
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-fold 0 \
-underFlowProtect 1 \
-rf "$REGIONS_FILE" \
-bam "$BAM_LIST" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-doSaf 1 \
-P 2

# make SFS file
realSFS "$DIR_OUT/$POP".saf.idx -fold 0 -P 2 > "$DIR_OUT/$POP".sfs

# make thetas
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-doThetas 1 \
-fold 0 \
-doSaf 1 \
-pest "$DIR_OUT/$POP".sfs \
-underFlowProtect 1 \
-rf "$REGIONS_FILE" \
-bam "$BAM_LIST" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-P 2


# calculate overall theta and (no windows!):
thetaStat do_stat "$DIR_OUT/$POP".thetas.idx -outnames "$DIR_OUT/$POP".thetas.all



echo "all done!"

# options
# basic quality filtering for reads
# -rf to specify a set of regions from a file
# anc polarizes SFS by reference genome .. but that's not a true ancestral polarization, 
# we need the unfolded SFS to get a 2D SFS estimate (-fold 0)
# but just know some stats more complicated than pi/Taj D require correct polarization
# underFlowProtect is necessary for large #s of bams
