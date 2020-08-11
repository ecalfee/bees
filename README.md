# Bees
This github repository accompanies the 2020 manuscript "Selection and hybridization shaped the Africanized honey bee invasion of the Americas" by Erin Calfee, Marcelo Nicolas Agra, Maria Alejandra Palacio, Santiago R. Ramirez and Graham Coop.

## Citation
https://www.biorxiv.org/content/10.1101/2020.03.17.994632v2

## Where to find scripts for main analyses
#### List of samples, by population
- bee_samples_listed/ # see also supplemental files from the manuscript
#### Bioinformatics pipeline fastq -> bams
- filtered_bams/
#### Identifying SNPs and calculating genotype likelihoods
- geno_lik_and_SNPs/
#### Global ancestry inference (NGSAdmix and PCA)
- global_ancestry/ # includes Fig 1 map
#### Local ancestry inference (ancestry\_hmm) and ancestry variance-covariances
- local_ancestry/ # includes Fig 3 and K matrix calculation
#### Genome-wide and individual clines
- clines/ # includes Fig 2 comparing ancestry clines to wing length clines
#### Wing lengths and admixture mapping
- wing_analysis/
#### Ancestry outlier regions and overlap with genes and QTLs
- functional_analysis/
#### Within ancestry diversity (pi, FST, PCA)
- within_ancestry/ # includes Fig 4
#### Mitochondria SNPs and clines
- clines/
#### Ancestry informative markers (AIMs)
- clines/
#### Other
- maps/ # additional scripts for mapping samples
- labwork/ # notes taken during labwork
- sims_downsample/ # preliminary analysis for depth of sequencing needed for local ancestry inference
- mitochrondria/ # attempted tree analysis for mtDNA (insufficient coverage to recover complete mtDNA haplotypes)
- colors.R # color palettes used in plots
