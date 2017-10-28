# bees
Genomic project analyzing admixture in Africanized honey bees

# data
data/bees_new_positions
For initial ancestry calling, I use a set of SNPs produced for me by Julie using her quality controls. 
These SNPs should include all high-quality sites polymorphic in the combined data for CA bees and European/African panels from Harpur.
Physical positions (bp) are from alignment to v4.5 of the honey bee genome, 
while "the map file is calculated based on the mean recombination rate per chromosome from the Beye 2006 paper.  So clearly it is an approximation."

# scripts
bee_filter.txt and make_bee_filter.R go with the CA/Harpur data and use the very approximate recombination distance noted above.
bee_elai has some command line scripts initially run on the data for threeway admixture between C M and A in CA bees.
plink_commands.txt are for general commands in plink

# downloading NCBI data - Harpur 2014 not sure they're actually paired reads
while read p; do fastq-dump --split-files -clip --gzip --skip-technical -O ../bees/data/Harpur_2014_NCBI/ "${p}"; done < ../bees/SRR_Acc_List_Harpur_2014.txt;
