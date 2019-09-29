# script uses beta regression to calculate ML cline for individual loci
#!/usr/bin/env Rscript
library(dplyr, warn.conflicts = F)
library(betareg)
library(tidyr)

# to run first 100 snps from observed A ancestry data, creates output results/ind_snp_clines/A_1_100.txt:
# Rscript betareg_clines_ind_loci.R A 1 100 A_1_100

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
DATA_NAME = args[1]
INDEX_START = as.integer(args[2])
INDEX_END = as.integer(args[3])
NAME_OUT = args[4]

# output to clines/results
# make sure output directory exists first (or create)!
dir_output = paste0("results/ind_snp_clines")
if (!file.exists(dir_output)){ 
  dir.create(file.path(dir_output), recursive = T)
}

# load population metadata
load("../local_ancestry/results/pops_by_lat.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 

# load mean ancestry by population
# .RData file contains objects 'sites' which has columns chr and pos and 'DATA_NAME' which has ancestry where sites are rows and pops are columns
load(paste0("../local_ancestry/results/", DATA_NAME, ".RData")) 

# subset data to pops of interest, only SA
A_snp_all <- get(DATA_NAME)[ , meta.AR.order.by.lat$population]
# deal with 0 and 1 values. I don't think zero-1 inflation should be modeled separately, as a diff. process
# in particular since this is almost entirely an issue for the simulated data only, not the observed data.
# so instead I will add/subtract a small value to get everything on (0, 1)
# an alternative would be to fit a zero-inflated mixed model
# e.g. see zoib() An R Package for Bayesian Inference for Beta Regression an dZero/One Inflated Beta Regression by Fang Liu and Yunchuan Kong
# zoib: https://journal.r-project.org/archive/2015/RJ-2015-019/RJ-2015-019.pdf
# or to transform all the data [0, 1] -> (0,1), e.g. (y(n−1) +0.5)/n, where n is the sample size (Smithson and Verkuilen 2006)
A_snp_all[A_snp_all == 0] <- 10^-5
A_snp_all[A_snp_all == 1] <- 1 - 10^-5
# in the observed data, there is only 1 value with A = 1 and none with A = 0
# but in the simulations, many more values hit the bounds.

# complain if index is out of range
if (INDEX_END > nrow(A_snp_all)) print(paste("! Warning end index", INDEX_END, 
                                             "exceeds rows in file:", nrow(A_snp_all), 
                                             ". This script will do calculations to the end of the rows"))
INDEX_END <- min(INDEX_END, nrow(A_snp_all))

# make function to repeatedly run betareg on many snps:
fit_betareg <- function(d){
  m <- betareg(A ~ abs_lat_SA_c, 
               data = d, 
               link = "logit")
  output <- c(m$coefficients$mean, m$coefficients$precision, loglik = m$loglik)
  names(output) <- c("mu", "b_lat", "phi", "loglik")
  return(output)}

# run ML beta regression model on each snp
fits <- do.call(rbind,
                lapply(INDEX_START:INDEX_END, # also works 
                       function(i) A_snp_all[i, ] %>%
                         tidyr::gather(., "population", "A") %>%
                         left_join(meta.AR.order.by.lat, ., by = "population") %>%
                         fit_betareg(.)))

# save txt file with snp information and beta regression output -- actually forget sites for now
#write.table(cbind(sites[INDEX_START:INDEX_END, c("chr", "pos")],
#                  fits),
write.table(cbind(snp_index = INDEX_START:INDEX_END, fits), 
           paste0(dir_output, "/", NAME_OUT, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")
