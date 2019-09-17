#!/bin/bash -l

# this script calls genotypes for reference bees A, C, M


# to run: ./call_genotypes.sh A combined_sept19 Group1

POP="${1}"
BAM_LIST="../bee_samples_listed/byPop/${POP}.bams"
PREFIX="${2}"
SCAFFOLD="${3}"
SNP_FILE="results/${PREFIX}/variant_sites/${SCAFFOLD}.var.sites"
REGIONS_FILE="results/${PREFIX}/variant_sites/${SCAFFOLD}.regions"
DIR_OUT="results/${PREFIX}/genotypes/${POP}"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REF="../data/honeybee_genome/Amel_HAv3.1.fasta"

echo "Calculating allele frequencies for POP: "$POP" scaffold $SCAFFOLD"

# make directory to store output (if doesn't yet exist)

mkdir -p "$DIR_OUT"

echo "finding site allele frequencies"
angsd -out "${DIR_OUT}/${SCAFFOLD}" \
-bam "$BAM_LIST" \
-ref "$REF" \
-rf "$REGIONS_FILE" \
-underFlowProtect 1 \
-remove_bads 1 \
-doGeno 2 \
-doPost 1 \
-geno_minDepth 6 \
-geno_minMM 0.6 \
-minMapQ 30 \
-doMajorMinor 3 \
-doMaf 1 \
-GL 1 \
-sites "${SNP_FILE}" \
-doCounts 1 \
-P 1

echo "done getting allele frequencies!"


# options
# basic quality filtering for reads
# -rf limits angsd to walking through indexed regions in the regions file, not the whole genome
# -doGeno 2 writes genotypes as 0, 1, 2 copies of the minor allele (or -1 for no genotype called)
# -doPost 1 estimates genotype posteriors using allele freq and HW as prior
# -geno_minDepth=6 requires a minimum depth of 6 reads to call a genotype.
# -geno_minMM=0.6 requires at minimum 60% of reads to match the major or minor allele.
# -GL 1 uses samtools GL method
# -doMajorMinor 3: takes major & minor allele from sites file
# ANGSD calculates freq of minor allele (ignoring all other non maj-min alleles seen)
# -doMaf 1 uses known major and minor alleles + gets a ML estimate of the allele freq based on EM and the genotype likelihoods (note: assumes HW). http://www.popgen.dk/angsd/index.php/Allele_Frequencies
# -minMapQ 30: filter out sites with low mapping quality. No filter for base/BAQ quality (error rates are incorporated into the model)
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)

# underFlowProtect is necessary for large #s of bams
