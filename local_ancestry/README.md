# Local Ancestry
## List of Scripts
#### Main script for plotting local ancestry results, e.g. Fig 5
- plotLocalAncestryTracts.R  
#### Main script for analyzing local ancestry results and running simulations
- local_ancestry_results_and_simulations.R
#### Ancestry variance-covariance calculations (K matrix)
- k_matrix.R
#### K matrix correlations by recombination rate quintile and chromosome
- K_by_r_and_chr.R
#### Calculating false-discovery rates for ancestry outliers
- calc_FDRs.R
#### Scripts to produce input files for ancestry\_hmm
- countReadsACGT.sh
- getMajorMinorReadCountsFromACGT.R
- thin_SNPs_add_genos_ancestry_hmm.R
#### Ancestry\_hmm
- see pipeline below for ancestry\_hmm command run
- coverageLocalAncSites.R # calculates coverage across sites used for ancestry\_hmm
#### Post-processing ancestry\_hmm results
- calcAncFromPost.R
- ancestry_sites_to_tracts.R # converts ancestry calls at individual SNPs to a bed windows file
#### Identifies high confidence homozygous ancestry tracts (> 0.8 posterior)
-getHighPosteriorTracts.sh
#### Get highest posterior ancestry state (MAP value for ancestry mapping)
- getHighestPostAncState.R

## Record of scripts run
- commands.txt # includes preliminary analyses not in the paper

## Additional notes for running ancestry\_hmm pipeline
Note: Directories and file names may need to be changed in several scripts.  
#### 2 possible pipelines
I have two different pipelines. The one I've run on the new genome HAv3.1 uses low-coverage data from the very beginning in calling SNPs, o everything is run out of ANGSD (no vcf).  
If using VCFs, there is an older pipeline that used SNPs and genotypes from high-coverage data only:  
- old_scripts/old_plink_commands.txt # note the filtering etc. has been updated in the newer pipeline. This earlier SNP set I used for Amel4.5 comes out of the pipeline from Cridland 2017

#### recombination map used:
I used the 10kb map from Wallberg (obtained via email correspondence):
apis_ALL.rates.bp1.w_headers.csv.windows.10000.csv
The map does not reach the ends of the chromosome. I extended the recombination rate from the most most distal recombination rate bin to the end of the chromosome:
map_rates_extended_10kb.bed

#### steps to creating ancestry_hmm input files and running local ancestry calling:
1. Get ACM genotypes using samtools GL for all variant sites (see geno_lik_and_SNPs/):
- call SNPs in ANGSD which creates a mafs.gz file (with SNPs) and GL file (for calculating global ancestry proportions in NGSAdmix)
calc_GL_all_SNPs.sh
- from SNPs, create a .var.sites file for input to ANGSD -sites. Index sites file and create a regions file for ANGSD too.  
```
geno_lik_and_SNPs$ for i in $(cat ../data/honeybee_genome/chr.names); do zcat results/combined_sept19/GL_by_scaffold/${i}.mafs.gz | \
tail -n +2 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > results/combined_sept19/variant_sites/${i}.var.sites; sleep 2s; \
angsd sites index results/combined_sept19/variant_sites/${i}.var.sites; done
```
- also make a regions file for each sites file (so angsd doesn't walk through the whole genome):
```
geno_lik_and_SNPs$ for i in $(cat ../data/honeybee_genome/chr.names); do cut -f1 results/combined_sept19/variant_sites/${i}.var.sites | sort | uniq > results/combined_sept19/variant_sites/${i}.regions; done
```
- call genotypes for A/C/M reference bees using ANGSD for all of these sites. Creates a .genos file:
call_genotypes.sh

2. Get map position for each SNP in sites files by running Rscript bp_to_r_WallbergHAv3.1.R (creates .rpos files), then
thin/filter SNPs running script thin_SNPs_add_genos_ancestry_hmm.R (I just ran interactively once):
- MAF in at least one group (A/C/M) >= .3 and
- all A/C/M ancestry groups having >=6 individuals with called genotypes.
- minimum .005 cM spacing between closest SNPs, based on recomb. map. This creates output files per chromosome, e.g. local_ancestry/results/SNPs/combined_sept19/Group1.var.sites and Group1.CMA.genos. In total, I keep about 1/7th of all SNPs, 542K SNPs

3. combine sites files into one file for all chromosomes, and indexed that sites file with angsd:
```
local_ancestry/results/SNPs/combined_sept19$ (for i in {1..16}; do cat Group$i.var.sites; done) > chr.var.sites
local_ancestry/results/SNPs/combined_sept19$ angsd sites index chr.var.sites
```
combined CMA allele counts into one file for all chromosome:
```
local_ancestry/results/SNPs/combined_sept19$ (for i in {1..16}; do cat Group$i.CMA.genos; done) > chr.CMA.genos
```
4. get read counts for every individual for all SNPs:
```
local_ancestry$ nohup parallel --joblog logs/ACGTcounts_MajMinCounts_combined_sept19.log --noswap --jobs 4 './countReadsACGT.sh {1} combined_sept19 chr; Rscript ./getMajorMinorReadCountsFromACGT.R {1} combined_sept19 chr' :::: .my_populations.list &> logs/ACGTcounts_MajMinCounts_combined_sept19.out &
```
- countReadsACGT.sh
- getMajorMinorReadCountsFromACGT.R

5. Make ploidy files for ancestry_hmm

6. combined site information and allele counts from A/C/M with individual read counts to make input files (one input file per population):
To do this, I paste together stripped files into correct order for ancestry_hmm input file  
Note: order is CMA and I reset recombination distance at each scaffold (not LG)  
Columns: chromosome, position_bp, allele counts A1 in C, allele counts A2 in C, allele counts A1 in M, allele counts A2 in M, allele counts A1 in A, allele counts A2 in A, distance in Morgans between previous marker and this one, read counts A1 in sample1, read counts A2 in sample1, read counts A1 in sample2 etc.  
```
local_ancestry$ mkdir -p results/ancestry_hmm/combined_sept19; for i in $(cat ../bee_samples_listed/byPop/combined_sept19_pops.list); do paste results/SNPs/combined_sept19/chr.CMA.genos $(for j in $(cut -f1 results/ancestry_hmm/"$i".ploidy); do echo results/SNPs/combined_sept19/countsMajMin/"$j".counts.txt; done) > results/ancestry_hmm/combined_sept19/"$i".counts; done
```

7. run ancestry_hmm for each population:
(estimates time but no bootstrap)
```
local_ancestry$ nohup parallel --noswap --jobs 4 --joblog logs/ancestry_hmm_combined_sept19.log \
'A=$(grep '{1}' {2} | cut -f2); \
C=$(grep '{1}' {2} | cut -f3); \
M=$(grep '{1}' {2} | cut -f4); \
echo pop: "{1}" A:"$A" C:"$C" M:"$M"; \
cd results/ancestry_hmm/combined_sept19/posterior; \
ancestry_hmm_v0.94 -e 3e-3 -a 3 "$C" "$M" "$A" \
-p 0 100000 "$C" -p 1 -100 "$M" -p 2 -60 "$A" \
--tmax 150 --tmin 2 --ne 670000 \
-i ../{1}.counts -s ../../{1}.ploidy' ::: populations
where {1} is the population name, and $A $C and $M are the populations admixture proportions (e.g. from STRUCTURE or in my case NGSAdmix)
```
NOTE: When you start a job, check if ancestry_hmm is actually running by seeing if it is using memory & processing power. There are no error messages when the input is wrong .. it will just spin in an infinite loop without using up memory or cpu. Otherwise it’s pretty quick!
