library(dplyr)
setwd("/media/erin/3TB/Documents/gitErin/bees/")
# See the results of NGSadmix below:
ngs_admix = t(as.matrix(read.table("data/out_all_NGSadmix.qopt"))) # get NGSadmix results
barplot(ngs_admix,col=1:3,space=0,border=NA,xlab="Individuals",ylab="admixture")
id<-read.table("Bee_bam_include_NGSadmix.txt",as.is=T) # read in sample IDs
id1 = lapply(pop$V1, function(x) strsplit(x, split="/")[[1]][[4]])
id2 = sapply(pop1, function(x) strsplit(x, split="[.]")[[1]][[1]])
# put sample IDs together with admixture results
admix = as.data.frame(list(id = pop2, A = ngs_admix["V1",], M = ngs_admix["V2",], C = ngs_admix["V3",]))
# add population labels
# CA metadata
metaCA = read.csv("data/CA_Bee/bam_files/Population_coverage", sep = "\t", header = F, stringsAsFactors = F)
colnames(metaCA) = c("id", "coverage", "pop", "year")
# Harpur metadata
metaHarpur = read.csv("data/Harpur_2014_NCBI/Harpur_SraRunTable.txt", sep = "\t", header = T, stringsAsFactors = F)
admix = left_join(admix, metaCA, by = "id")
admix = metaHarpur %>%
  mutate(., location = geographic_location_s) %>%
  mutate(., id = paste0("Harpur_", Run_s)) %>%
  mutate(., study = "Harpur") %>%
  select(., one_of(c("label_s", "location", "id", "study"))) %>%
  left_join(admix, ., by = "id")
admix[!is.na(admix$year), "study"] = "CA"
admix[is.na(admix$study), "study"] = "Kenya"
admix[, "pop"] = ifelse(!is.na(admix$pop), admix$pop, 
                        ifelse(!is.na(admix$label_s), 
                               substring(admix$label_s, 1, 1), 
                                        "A"))
admix = arrange(admix, year) %>% # organize by population and within population, by year
  arrange(., pop)
barplot(t(as.matrix(admix[,c("A", "C", "M")])),col=c("red", "green", "blue"),space=0,border=NA,xlab="Individuals",ylab="admixture")
text(tapply(1:length(admix$pop),admix$pop,mean),-0.05,unique(admix$pop),xpd=T)
# I accidentally included a population from Yemen in the Harpur data (O lineage - in #7 in 8th row)
admix[,c("pop")]
# I can repeat this with just 2 populations to confirm A and C/M stand out (check sensitivity)
# I can also repeat with some thinning to see if the results are the same