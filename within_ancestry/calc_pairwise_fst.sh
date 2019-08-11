#!/bin/bash -l

# this script calculates fst for a pair of populations
# individual population sfs and saf files (e.g. made in calc_theta_regions.sh)
# must be prepared ahead of time and found in the directory $DIR


# to run: ./calc_pairwise_fst.sh A SA_3 results/outlier_regions/region_5/combined 

POP1="${1}"
POP2="${2}"
DIR="${3}"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

echo "Calculating fst for POP1 $POP1 and POP2 $POP2"
echo "Directory: $DIR"

echo "getting 2D SFS"

realSFS "$DIR"/"$POP1".saf.idx "$DIR"/"$POP2".saf.idx -fold 1 -P 2 > "$DIR"/"$POP1"-"$POP2".folded.sfs

echo "now making Fst"
realSFS fst index "$DIR"/"$POP1".saf.idx "$DIR"/"$POP2".saf.idx \
-sfs "$DIR"/"$POP1"-"$POP2".folded.sfs \
-whichFst 1 \
-fold 1 \
-fstout "$DIR"/"$POP1"-"$POP2"

# print outputs for fst:
realSFS fst stats "$DIR"/"$POP1"-"$POP2".fst.idx > "$DIR"/"$POP1"-"$POP2".fst.stats
realSFS fst stats2 "$DIR"/"$POP1"-"$POP2".fst.idx > t"$DIR"/"$POP1"-"$POP2".fst.all
realSFS fst stats2 "$DIR"/"$POP1"-"$POP2".fst.idx -win 5000 -step 1000 > "$DIR"/"$POP1"-"$POP2".fst.windows

echo "all done!"

# options
# all setting are in the making of the sfs and saf.idx files
# this script uses the folded site frequency spectrum to calculate pairwise population Fst
