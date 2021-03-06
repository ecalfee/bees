#!/bin/bash -l

# this script calculates pi for a population and list of bams, e.g. A ref pop
# over a small region, e.g. an outlier


# to run: ./calc_theta_regions.sh A ../bee_samples_listed/byPop/A.list Group1.23:57322-186870 5

POP="${1}"
BEE_ID_FILE="${2}"
REGION="${3}"
REGION_N="${4}"
DIR_BAMS="../filtered_bams/results_Amel4.5"
DIR_OUT="results/outlier_regions/region_${REGION_N}/combined"
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
-r "$REGION" \
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
-r "$REGION" \
-bam "$BAM_LIST" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-P 2


# calculate overall theta and windows:
thetaStat do_stat "$DIR_OUT/$POP".thetas.idx -outnames "$DIR_OUT/$POP".thetas.all
thetaStat do_stat "$DIR_OUT/$POP".thetas.idx -win 5000 -step 1000 -outnames "$DIR_OUT/$POP".thetas.windows



echo "all done!"

# options
# basic quality filtering for reads
# -r to specify just a small region
# anc polarizes SFS by reference genome ..
# but that's not a true ancestral polarization, so we cant interpret some stats like FaysH
# we have to use the unfolded sfs for Fst calcs later though (-fold 0)
# underFlowProtect is necessary for large #s of bams
