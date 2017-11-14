setwd("/media/erin/3TB/Documents/gitErin/bees")
# load in mapQ
mapQ = read.csv("data/Harpur_2014_NCBI/alignment_stats/Harpur_SRR957060.sorted.bam.mapq")
mapQ = read.csv("data/Kenya_Sheppard_NCBI/alignment_stats/Kenya_SRR5270368.sorted.bam.mapq")
mapQ = read.csv("data/CA_Bee/alignment_stats/SRCD10B.merged.bam.mapq") # peaks >40 and smaller peak at 10
mapQ = read.csv("data/CA_Bee/alignment_stats/SRCD49C.merged.bam.mapq") # peaks at 44 and 10
mapQ = read.csv("data/CA_Bee/alignment_stats/ap41.sequence.bam.mapq") # peaks at 0: much lower depth and mapping quality
mapQ = read.csv("data/CA_Bee/alignment_stats/ap42.sequence.bam.mapq") # peaks at 0: again, poorer dist. of mapping quality
# plot
names(mapQ) = "mq"
hist(mapQ$mq, breaks = 44, freq = F)
density(mapQ$mq)
hist(mapQ$mq[(mapQ$mq<=42) & (mapQ$mq>=2)], breaks = 41, freq = F)
table(mapQ)
# bimodal when you zoom out, but not clear where a cutoff should be -- 10, 20?
# other option: https://academic.oup.com/bioinformatics/article/28/18/i349/249968
ps = c(.1, .01, .05, .001)
plot(ps, -10*log(ps, base = 10))
