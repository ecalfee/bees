# Within ancestry diversity analyses
## List of Scripts
#### Main script to plot diversity genomewide
- plot_pi_genomewide.R # Fig 4 and S25
#### Run block bootstrap of pi within ancestry
- block_bootstrap_pi_within.R # note: geno_lik_and_SNPs/bp_to_r_WallbergHAv3.sh returns 1cM bins across genome
#### Script to plot Fst across outlier regions (supporting figs)
- plot_fst_outliers.R
#### Script to create new bams that only include reads overlapping high confidence AA, CC, MM ancestry tracts
- getBamsForTracts.sh
#### Scripts to calculate population allele frequencies
- allele_freq_within_ancestry_outliers.sh
- allele_freq_within_ancestry.sh
- allele_freq.sh # ignores ancestry state ('Combined')
- combine_pop_allele_freqs.sh # ensures every site has an allele freq (can be NA for no data)
#### Run and plot PCA analysis (PCAngsd) within ancestry
- calc_GL_within_ancestry_4PCA.sh
- runPCAngsd_within_ancestry.sh
- plotPCA_within_ancestry.R
#### Functions for Fst and pi (small sample size correction)
- het_fst_functions.R
#### Calculates total mappable sites in the genome (used in denominator for pi)
- calc_frac_snps.R

## Record of scripts run
- commands.txt # includes preliminary analyses not in the paper
