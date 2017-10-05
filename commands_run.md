# list of commands run
# Starting with called genotype SNP data curated by Julie and in data/new_bee_positions/All.bed

# Thu Oct  5 12:24:55 PDT 2017
Reduced SNPs and data to polymorphic sites in a list of the families to include in my analysis. CA2014 bees and A,C,M
Option mac 2 keeps only SNPs with minimum count of 2 in new subset
```
plink1.9 --bfile ../bees_new_positions/ALL --mac 2 --keep-fam CA2014_AMC.txt --make-bed --out CA2014_AMC
```

Get subset of SNPs in low LD within A family; about 10% of SNPs
--indep window_size (in SNPs) window_shift (in SNPs) VIF
Where VIF is recommended for smaller samples sizes 1.5-2
and VIF = 1/(1-R^2) where R^2 is the multivariate correlation across a whole window of SNPs whereas r^2 is pairwise
Because the data is unphased, it's really a genotype correlation coefficient
```
plink1.9 --bfile CA2014_AMC --keep-fam A.txt --indep 50 5 2 --out low_LD_A2
```
Also ran the pruning with the C and M groups for similar results in # SNPs kept.
And with the whole sample which allowed more SNPs to be kept (likely because of higher sample sizes?)
And with an admixed population, Pacerita_2014, which I was surprised wasn't excluding many more SNPs for LD.
These results are easy to regenerate and in extra/ . Setting VIF to 1 left very few SNPs.
The reason for SNP pruning is NGSadmix does not handle within ancestry LD or admixture LD.

Convert to beagle format to use in NGSadmix. One file with chromosomes combined.
```
plink1.9 --bfile CA2014_AMC --extract low_LD_A2.prune.in --recode beagle-nomap --out low_LD_A2_CA2014_AMC
```
This produced file data/CA2014_AMC/low_LD_A2_CA2014_AMC.beagle.dat
...So this is a beagle genotype data when NGSadmix needs a beagle genotype likelihood file
https://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf
From a BAM file I could create this using ANGST, but from plink I'm not sure
How is missing data encoded anyways?
Mega2 file type conversion may be useful down the road.
I have decided to just start with the aligned read data (BAM format) instead of called SNPs.

