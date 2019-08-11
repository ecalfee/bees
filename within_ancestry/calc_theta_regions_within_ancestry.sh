#!/bin/bash -l

# this script calculates pi for a population and list of bams, within an ancestry type, e.g. AA
# over a small region, e.g. an outlier


# to run: ./calc_theta_regions_within_ancestry.sh NA_3 ../bee_samples_listed/byPop/NA_3.list Group1.23:57322-186870 5 AA results/thin1kb_common3/byPop/output_byPop_CMA_ne670000_scaffolds_Amel4.5_noBoot/AA/bams

POP="${1}"
BEE_ID_FILE="${2}"
REGION="${3}"
REGION_N="${4}"
ANCESTRY="${5}"
DIR_BAMS="${6}"
DIR_OUT="results/outlier_regions/region_${REGION_N}/${ANCESTRY}"
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
-fold 1 \
-underFlowProtect 1 \
-r "$REGION" \
-bam "$BAM_LIST" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-doSaf 1 \
-P 2

# make folded SFS file
realSFS "$DIR_OUT/$POP".saf.idx -fold 1 -P 2 > "$DIR_OUT/$POP".folded.sfs

# make thetas
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-doThetas 1 \
-fold 1 \
-doSaf 1 \
-pest "$DIR_OUT/$POP".folded.sfs \
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
# anc polarizes SFS by reference genome .. but that's not a true ancestral polarization, so we use the folded sfs instead (-fold 1)
# underFlowProtect is necessary for large #s of bams
