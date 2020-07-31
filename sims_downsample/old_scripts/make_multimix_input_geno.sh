#!/bin/bash
# to run: bees/sims_downsample$ type=admix bim=../data/bees_new_positions/ALL chr=16 pop=pop1 out=../data/sims_downsample/multimix/input ./make_multimix_input_geno.sh
# note: this script counts the A1 allele, not A2 (A2 appears to be the major allele in the plink files)
echo "making multimix input genotype file."
echo "population: "$pop
echo "bim file: "$bim
echo "chromosome: "$chr
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
mkdir -p $out/temp_$pop"_"chr$chr/
# get list of all individual IDs to process
# first make an individuals file if it doesn't yet exits:
# with population ID, individual ID, then individual ID again in 3rd column
if ! [ -r $bim.ind ]
then
    awk '{print $1, '\t', $2, '\t', $2}' $bim.fam > $bim.ind
else
    echo "file already exists, won't rewrite: "$bim.ind
fi

IDs=$(awk -v pop=$pop '$1 == pop {print $2}' $bim.ind)
echo "IDs to include: " $IDs


# get individual counts for all SNPs in bim file; write to temp directory $id.frq.counts
for id in $IDs; \
do plink --bfile $bim --freq counts --keep-allele-order --filter $bim.ind $id \
--out $out/temp_$pop"_"chr$chr/$id; \
done



# make output & output directories
if [ $type == "ref" ]
then
    # make output directory
    mkdir -p $out/input/haplotypes/$pop
    # save IDs in order
    echo $IDs > $out/input/haplotypes/$pop/ref.ind;
    # make each individual's output
    for id in $IDs;
    # counts of A1 haplotype (column 5)
    do awk '$1 != "CHR" { if ($7 == 0) print $5; else print 9 }' $out/temp_$pop"_"chr$chr/$id.frq.counts > $out/temp_$pop"_"chr$chr/$id.strp;
    done
    # paste it all together
    paste -d" " $(for id in $IDs; do ls $out/temp_$pop"_"chr$chr/$id.strp; done) > $out/input/haplotypes/$pop/chr$chr.genos
fi

if [ $type == "admix" ]
then
    # make output directory admixed samples
    mkdir -p $out/input/Samples/genos/$pop
    # save IDs in order
    echo $IDs > $out/input/Samples/genos/$pop/admix.ind;
    # make each individual's output
    for id in $IDs;
    # counts of A1 haplotype (column 5)
    do awk '$1 != "CHR" { if ($7 > 0) print "9 9 9"; else if ($5 == 0) print "0 0 1"; else if ($5 == 1) print "0 1 0"; else print "1 0 0"}' \
$out/temp_$pop"_"chr$chr/$id.frq.counts > $out/temp_$pop"_"chr$chr/$id.strp;
    done

    # make SNP ID output
    # get this part from bp_scaffold_to_bp_LG_4multimix.R for absolute bp positions w/ 50kb gaps
    # paste it all together
    paste -d" " $out/input/haplotypes/legend_files/chr$chr.legend $(for id in $IDs; do ls $out/temp_$pop"_"chr$chr/$id.strp; done) > $out/input/Samples/genos/$pop/genos_chr$chr
fi

# remove temporary directory
rm -r $out/temp_$pop"_"chr$chr

