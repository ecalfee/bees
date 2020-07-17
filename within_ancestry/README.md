# Within ancestry diversity (and overall pi)
## List of Scripts

#### Main R script to make plots of diversity genomewide (e.g. Fig 4)
- plot_pi_fst_from_allele_freqs.R

#### Get individual bams for high-confidence homozygous ancestry regions
Using local_ancestry scripts I called high confidence AA/CC/MM HOMOZYGOUS ANCESTRY tracts for all bees using local_ancestry/getHighPosteriorTracts.sh, e.g. local_ancestry/results/tracts/combined_sept19/A/CA1207.bed. The script below filters bams to only save reads that overlap these tracts.
- getBamsForTracts.sh

#### Estimate allele frequencies at all SNPs on a chromosome
Note this uses a sites file for all SNPs created in *** and has no output if there are no reads covering a SNP.
##### For high-confidence homozygous ancestry tracts (within ancestry)
- allele_freq_within_ancestry.sh
##### For all regions (uses full original bams)
- allele_freq.sh

#### Create one file per population that combines allele frequencies from all sites and all chromosomes
- combine_pop_allele_freqs.R

#### Estimate bootstrap uncertainty around pi estimates
- block_bootstrap_pi_within.R
- bootstrap_pi_within.R

#### Pi and Fst within outlier regions:


## Record of scripts run
- commands.txt # includes preliminary analyses not in the paper


## Other scripts
-allele_freq.sh
-allele_freq_within_ancestry.sh
-allele_freq_within_ancestry_outliers.sh
-calc_GL_within_ancestry_4PCA.sh
-calc_pairwise_fst.sh
-calc_theta_random_background.sh
-calc_theta_regions.sh
-calc_tehta_regions_within_ancestry.sh
-combine_pop_allele_freqs.R
extepected_drift_ancestry.R
pi.sh
pi_within_ancestry.sh
plotPCA_within_ancestry.R
plot_pi_fst_outliers.R
runPCAngsd_within_ancestry.sh
