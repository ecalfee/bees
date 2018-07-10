#!/bin/bash
# to run: bees/sims_downsample$ type=admix bim=../data/bees_new_positions/ALL pop=pop1 out=../data/sims_downsample/multimix/input ./make_multimix_input_geno.sh
echo "making multimix input genotype file."
echo "population: "$pop
echo "bim file: "$bim
echo "type of output: "$type
echo "to output directory: "$out
# check for input errors
if ! [ -r $bim.bim ]
then
    echo "error: bim file doesn't exist: "$bim".bim" >&2
    exit 1
fi
if [ $type != "admix" ] && [ $type != "ref" ]
then
    echo "error: type should be ref or admix" >&2
    exit 1
fi
if [ $type == "admix" ] && ! [ -r $out/input/haplotypes/legend_files/chr1.legend ]
then
    echo "error: you need to create legends file first before admix files (use bp_scaffold_to_bp_LG_4multimix.R)"
    exit 1
fi

# make temporary directory
mkdir -p $out/temp/$pop
# get list of all individual IDs to process
# first make an individuals file if it doesn't yet exits:
# with population ID, individual ID, then individual ID again in 3rd column
if ! [ -r $bim.ind ]
then
    awk '{print $1, '\t', $2, '\t', $2}' $bim.fam > $bim.ind
else
    echo "file already exists, won't rewrite: "$bim.ind
fi

IDs=$(awk -v pop=$pop '$1 == $pop {print $2}' $bim.ind)
CHROMOSOMES={1..16}

# get individual counts for all SNPs in bim file; write to temp directory $id.frq.counts
for id in $IDs; \
do plink --bfile $bim --freq counts --keep-allele-order --filter individuals.txt $ind \
--out $out/temp/$pop/$id; \
done



# make output & output directories
if [ $type == "ref" ]
then
    # make output directory
    mkdir -p $out/input/haplotypes/$pop
    # make each individual's output
    for chr in $CHROMOSOMES;
    do for id in $IDs;
    # counts of A1 haplotype (column 5)
    do awk -v chr=$chr '$1 == $chr { if ($7 == 0) print $5; else print 9 }' > $out/temp/$pop/$id.chr$chr.strp;
    echo $id >> $out/input/haplotypes/$pop/ref.ind; # save IDs in order
    done
    # paste it all together
    paste -d" " $(for id in $IDs; do ls $out/temp/$pop/$id.chr$chr.strp; done) > $out/input/haplotypes/$pop/chr$chr.genos
fi

if [ $type == "admix" ]
then
    # make output directory admixed samples
    mkdir -p $out/input/Samples/genos/$pop
    # make each individual's output
    for chr in $CHROMOSOMES;
    do for id in $IDs;
    # counts of A1 haplotype (column 5)
    do awk -v chr=$chr '$1 == $chr { if ($7 > 0) print 9 9 9; elif ($5 == 0) print 0 0 1; elif ($5 == 1) print 0 1 0; else print 0 0 1}' > $out/temp/$pop/$id.chr$chr.strp;
    echo $id >> $out/input/haplotypes/$pop/ref.ind; # save IDs in order
    done

    # make SNP ID output
    # get this part from bp_scaffold_to_bp_LG_4multimix.R for absolute bp positions w/ 50kb gaps
    # paste it all together
    paste -d" " $out/input/haplotypes/legend_files/chr$chr.legend $(for id in $IDs; do ls $out/temp/$pop/$id.chr$chr.strp; done) > $out/input/Samples/genos/$pop/genos_chr$chr.genos
fi

# remove temporary directory
rm -r $out/temp/$pop

