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

if [[ "$POP1" =~ ^(A|C|M)$ ]]; then
    echo "pop1 is A/C/M go to combined ancestry directory for saf.idx"
    DIR1="$DIR"/../combined
else
    DIR1="$DIR"
fi
if [[ "$POP2" =~ ^(A|C|M)$ ]]; then
    echo "pop2 is A/C/M go to combined ancestry directory for saf.idx"
    DIR2="$DIR"/../combined
else
    DIR2="$DIR"
fi


echo "getting 2D SFS"

realSFS "$DIR1"/"$POP1".saf.idx "$DIR2"/"$POP2".saf.idx -fold 0 -P 2 > "$DIR"/"$POP1"-"$POP2".sfs

echo "now making Fst"
realSFS fst index "$DIR1"/"$POP1".saf.idx "$DIR2"/"$POP2".saf.idx \
-sfs "$DIR"/"$POP1"-"$POP2".sfs \
-whichFst 1 \
-fold 0 \
-fstout "$DIR"/"$POP1"-"$POP2"

# print outputs for fst:
realSFS fst stats "$DIR"/"$POP1"-"$POP2".fst.idx > "$DIR"/"$POP1"-"$POP2".fst.stats
realSFS fst stats2 "$DIR"/"$POP1"-"$POP2".fst.idx > "$DIR"/"$POP1"-"$POP2".fst.all
realSFS fst stats2 "$DIR"/"$POP1"-"$POP2".fst.idx -win 5000 -step 1000 > "$DIR"/"$POP1"-"$POP2".fst.windows

echo "all done!"

# options
# all setting are in the making of the sfs and saf.idx files
# this script uses the unfolded site frequency spectrum to calculate pairwise population Fst
# the SFS is arbitrarily polarized by the reference genome, which makes not diff. for pi and fst but some stats can't be used
