library(dplyr)
library(tidyr)
library(ggplot2)
  
# plot clines in allele frequencies at outlier loci
# to see snps underlying ancestry skews and confirm ancestry calls are consistent wiht allelic clines
source("../colors.R") # get color palette
load("../local_ancestry/results/pops_by_lat.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 

# first load SNPs with ancestry calls:
load("../local_ancestry/results/sites_r.RData") # loads sites dataframe
# and find the outlier loci.

# load allele freq data for each population at all AIMs

# then for chr11 find nearest M AIMs for outlier region & plot the allele freq clines NA and SA

# for chr1 find nearest A AIMs for outlier regions & plot the allele freq clines NA and SA
