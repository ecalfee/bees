sfs <- data.frame(Riverside_2014 = t(read.table("test/small_Riverside_2014.sfs", header = F)),
                       CA04 = t(read.table("test/small_CA04.sfs", header = F)),
                  derived_count = 0:16)
sfs$freq = sfs$derived_count/16
sfs$het = sfs$freq*(1-sfs$freq)
# pi for Riverside_2014:
sum(sfs$freq*sfs$Riverside_2014)/sum(sfs$Riverside_2014)
# pi for CA04:
sum(sfs$freq*sfs$CA04)/sum(sfs$CA04)
sfs
png("plots/bias_angsd_sfs_example_low_vs_high_coverage.png", height = 4, width = 5.2, units = "in", res = 600)
with(sfs, plot(derived_count, Riverside_2014/sum(Riverside_2014), col = "orange", pch = 1,
               xlab = "Non-reference allele count", ylab = "Frequency", main = "SFS estimated by ANGSD"))
with(sfs, points(derived_count, CA04/sum(CA04), col = "blue", pch = 3))
legend("topright", legend = c(paste0("Riverside_2014 (high coverage), pi=", round(sum(sfs$freq*sfs$Riverside_2014)/sum(sfs$Riverside_2014), 3)), 
                              paste0("CA04 (low coverage), pi=", round(sum(sfs$freq*sfs$CA04)/sum(sfs$CA04), 3))), 
       col = c("orange", "blue"), pch = c(1, 3))
dev.off()

png("plots/bias_angsd_sfs_example_low_vs_high_coverage_zoom_in_tail.png", height = 4, width = 5.2, units = "in", res = 600)
with(sfs, plot(derived_count[4:17], Riverside_2014[4:17]/sum(Riverside_2014), col = "orange", pch = 1,
               xlab = "Non-reference allele count", ylab = "Frequency", main = "SFS estimated by ANGSD, zoom in 3+ count"))
with(sfs, points(derived_count, CA04/sum(CA04), col = "blue", pch = 3))
legend("topright", legend = c(paste0("Riverside_2014 (high coverage), pi=", round(sum(sfs$freq*sfs$Riverside_2014)/sum(sfs$Riverside_2014), 3)), 
                              paste0("CA04 (low coverage), pi=", round(sum(sfs$freq*sfs$CA04)/sum(sfs$CA04), 3))), 
       col = c("orange", "blue"), pch = c(1, 3))
dev.off()


# conclusion: ANGSD gives biased SFS for low coverage populations, calling more SNPs than in similar
# high coverage populations, and therefore inflating diversity (pi) estimates for low coverage pops