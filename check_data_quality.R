setwd("/media/erin/3TB/Documents/gitErin/bees")
# load in mapQ
mapQ = read.csv("data/Harpur_2014_NCBI/alignment_stats/Harpur_SRR957060.sorted.bam.mapq")
mapQ = read.csv("data/Kenya_Sheppard_NCBI/alignment_stats/Kenya_SRR5270368.sorted.bam.mapq")
# plot
names(mapQ) = "mq"
hist(mapQ$mq, breaks = 44, freq = F)
hist(mapQ$mq[(mapQ$mq<=42) & (mapQ$mq>=2)], breaks = 41, freq = F)
# bimodal when you zoom out, but not clear where a cutoff should be -- 10, 20?
# other option: https://academic.oup.com/bioinformatics/article/28/18/i349/249968