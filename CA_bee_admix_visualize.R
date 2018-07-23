library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(data.table)
source("../zea/erin_scripts/elai_output_functions.R")
source("../covAncestry/forqs_sim/k_matrix.R") # get functions for k matrix and covariance in ancestry analysis

#test1 = get_ELAI_ancestry("out_chr15_g30_c9_s1", "data/ELAI_input/output/", n_anc = 3)

# use fread from data.table because much faster than read.table, then convert to data.frame
#d1 = read.table("data/ELAI_input/output/out_chr15_g30_c9_s1.ps21.txt", header = F, sep = " ")
n_anc = 3
chr = 5
g = 60
#g=30 differences in generations doesn't appear to make a difference, though g30 is possibly more conservative in population ancestry
d1 = fread(paste0("data/ELAI_input/output/out_chr", chr, "_g", g, "_c9_s1.ps21.txt"))
d2 = as.data.frame(t(d1))[c(T, F, F), ]
# plot(d2[1:50,1], d2[1:50,2])
# test it's the same as using perl to get every 3rd entry (A ancestry contribution)
test = as.data.frame(t(fread(paste0("data/ELAI_input/output/A_chr", chr, "_g", g, "_c9_s1.all.txt"))))


d3 = d2[, c(2,3,4,9,12,13,1,5,6,7,8,10,11,14)]
apply(d3, 2, mean)/2 # mean A ancestry in this region
corrplot(cov2cor(make_K_matrix(t(d3))$K), method="shade")
# 6 are placerita and 8 are riverside (all 2014, see recomb_map.R for list/order)
#dp = apply(d2[,c(2,3,4,9,12,13)], 1, mean)
#dr = apply(d2[,c(1,5,6,7,8,10,11,14)], 1, mean)
dp = apply(d3[,1:6], 1, mean)
dr = apply(d3[,7:14], 1, mean)
#dp = apply(d2[,c(1,3,5,7,11,13)], 1, mean)
#dr = apply(d2[,c(2,4,6,8,10,12,14)], 1, mean) # slight + cov when divided randomly (good check)
dpop = rbind(dp, dr)
corrplot(cov2cor(make_K_matrix(dpop)$K), method = "shade")

d1_pos = fread(paste0("data/ELAI_input/output/out_chr", chr, "_g", g, "_c9_s1.snpinfo.txt"))$pos
d2$pos = d1_pos
plot(d2$pos, d2$V1, type = "l")
plot(d2$pos, dp, type = "l")
lines(d2$pos, dr, type = "l", col = "blue")
plot(d2$pos, apply(dpop, 2, mean), type = "l", col = "purple")
#lines(d2$pos, apply(dpop, 2, mean), type = "l", col = "blue")



# get recombination rates along chromosome and test for correlation with African ancestry
# get mean rate per chrom to translate pos back to bp value
avg_r_per_chr = read.csv(file = "data/recomb_map/Liu_2015/avg_r_per_chr_Liu_2015.csv", header = T, stringsAsFactors = F)
r = read.csv(file = "data/Wallberg_honeybee_recombination_rates_2015/recombination_rates/A.rates.1000.201.low_penalty.csv.cM_Mb.windows.100000.csv",
             sep = "\t", header = F, stringsAsFactors = F)
r_avg = avg_r_per_chr[avg_r_per_chr$chr_n == chr, "avg_r_cMperMb"] # get avg. rate for current chromosome
d2$pos_bp = round(d2$pos/r_avg)
d2$r_group = floor(d2$pos_bp/100000)
colnames(r) = c("chr_group", "start_r_group", "r_group", "r")
d4 = left_join(d2, r[r$chr_group == paste0("Group", chr), c("r_group", "r")], by = c("r_group"))
Abyr = tapply(apply(dpop, 2, mean), INDEX = floor(d4$r), mean)
plot(names(Abyr), Abyr)
abline(lm(Abyr ~ as.numeric(names(Abyr))), col = "blue")
# not summarized
plot(apply(dpop, 2, mean) ~ floor(d4$r))
abline(lm(apply(dpop, 2, mean) ~ floor(d4$r)), col = "blue")

