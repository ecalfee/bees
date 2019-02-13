#!/bin/bash

# this script creates input for MOSAIC based on fastPHASE haplotypes
PREFIX=$1
CHR=$2
DIR_OUT="results/MOSAIC/$PREFIX"
fastPHASE_FILE="results/fastPHASE/${PREFIX}_01.${CHR}_hapguess_switch.out"

# to run: ./transpose_fastPHASE_to_MOSAIC_input.sh ACM_Riv2014

# make output directory
mkdir -p $DIR_OUT 

# first transpose the data from haplotype rows X SNP columns
# to SNP rows X haplotypes columns
colN=$(head -n 1 "${fastPHASE_FILE}" | wc -m)
echo "number of columns = SNPs = $colN"

for n in $(seq 1 ${colN})
do cut -c"$n" "$fastPHASE_FILE" | \
cat -v | sed "s:\^@:?:g" | \
paste -s -d"\0"
done > "${DIR_OUT}/allgenotypes.$CHR"

# then split this transposed file into multiple files, one per group of haplotypes (population)
echo "splitting genotypes into one file per population"
cp "results/geno/${PREFIX}_var_sites.fam" "${DIR_OUT}/sample.names"


for j in $(cut -d" " -f1 "${DIR_OUT}/sample.names" | sort | uniq)
do echo $j
# creates comma separated list of haplotypes to take:
cut -c"$(for i in $(cut -d" " -f1 "${DIR_OUT}/sample.names" | grep -n "$j" | cut -d":" -f1)
	do echo $(( ($i * 2) - 1 )); echo $(( $i * 2 ))
	done | paste -s -d",")" "${DIR_OUT}/allgenotypes.$CHR" > "${DIR_OUT}/${j}genofile.$CHR" 
done


# options
# cut -c cuts out column character i
# cat -v makes special characters visible, then I use
# sed to replace special character ^@ in output with ? for missing genotype
# where \ is the leading escape character for ^ in ^@
# paste -s pastes character column i together on one line 
# -d"\0" specifies NULL deliminator separating SNPs
