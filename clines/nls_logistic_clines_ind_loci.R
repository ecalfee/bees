# script uses nonlinear least squares to calculate
# best fitting logistic cline for individual loci
# using population A ancestry frequencies in AR

#!/usr/bin/env Rscript
library(dplyr, warn.conflicts = F)
library(tidyr)

# to run first 100 snps from observed A ancestry data, creates output results/ind_snp_clines/A_1_100.txt:
# Rscript nls_logistic_clines_ind_loci.R A 1 100 A_1_100

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
DATA_NAME = args[1]
INDEX_START = as.integer(args[2])
INDEX_END = as.integer(args[3])
NAME_OUT = args[4]

# output to clines/results
# make sure output directory exists first (or create)!
dir_output = paste0("results/ind_snp_nls_clines")
if (!file.exists(dir_output)){ 
  dir.create(file.path(dir_output), recursive = T)
}

# load population metadata
load("../local_ancestry/results/meta.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 

# load mean ancestry by population
# .RData file contains objects 'sites' which has columns chr and pos and 'DATA_NAME' which has ancestry where sites are rows and pops are columns
load(paste0("../local_ancestry/results/", DATA_NAME, ".RData")) 

# subset data to pops of interest, only SA
A_snp_all <- get(DATA_NAME)[ , meta.AR.order.by.lat$population]

# complain if index is out of range
if (INDEX_END > nrow(A_snp_all)) print(paste("! Warning end index", INDEX_END, 
                                             "exceeds rows in file:", nrow(A_snp_all), 
                                             ". This script will do calculations to the end of the rows"))
INDEX_END <- min(INDEX_END, nrow(A_snp_all))

# define logistic curve
logistic3 <- function(x, mu, b){
  1/(1 + exp(-b*(x - mu)))
}

# make function to repeatedly run nls on many snps:
fit_cline0 <- function(snp){
  start0 <- getInitial(A ~ SSlogis(lat, Asym, 
                                   xmid, scal), 
                       d = snp)
  fit0 <- nls(A ~ logistic3(x = lat, b = b, mu = mu),
              start = list(b = 1/unname(start0["scal"]),
                           mu = unname(start0["xmid"])),
              data = snp,
              trace = F)
  sum0 <- summary(fit0)
  fit <- c(converged = sum0$convInfo$isConv,
           mu = sum0$coefficients["mu", "Estimate"],
           b = sum0$coefficients["b", "Estimate"],
           residual_error = sum0$sigma)
  return(fit)
}
# wrapper to return an error message and empty vector if fit fails
fit_cline <- function(snp){ 
  tryCatch(fit_cline0(snp),
           error=function(e) {
             print(e)
             return(c(converged = NA,
                      mu = NA,
                      b = NA,
                      residual_error = NA))
           }) 
}

# run nls logistic model fit on each snp
fits <- do.call(rbind,
                lapply(INDEX_START:INDEX_END, # also works 
                       function(i) A_snp_all[i, ] %>%
                         tidyr::gather(., "population", "A") %>%
                         left_join(meta.AR.order.by.lat, ., by = "population") %>%
                         fit_cline(.)))


# save txt file output
write.table(cbind(snp_index = INDEX_START:INDEX_END, fits), 
           paste0(dir_output, "/", NAME_OUT, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")
