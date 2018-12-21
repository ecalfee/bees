# calculation/plotting global ancestry for Ramirez CA bees and Harpur A, C, M bees
# for getting Riverside 2014 population alpha for input to ancestry_hmm
admix<-as.matrix(read.table("data/sims_downsample/NGSadmix/CA_Bee_AMC.qopt"))
ids <- read.table("data/bees_new_positions/thin5KB_ids.list", stringsAsFactors = F)
d = cbind(ids, admix)

colnames(d) = c("pop", "id", "kM", "kC", "kA")
d$pop = factor(d$pop, levels = c("A", "C", "M", "Riverside_2014", 
                                     "Avalon_2014", "Placerita_2014",
                                     "Stanislaus_2014", "Stebbins_2014",
                                     "Davis_2014", "Humboldt_2014"))
d = d[order(d$pop), ]
barplot(t(d[, 3:ncol(d)]),col=c("yellow", "blue", "red"),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
legend("topleft", legend = c("M", "C", "A"), col = c("yellow", "blue", "red"), pch = 19)

tapply(d$kA, d$pop, mean)
tapply(d$kC, d$pop, mean)
tapply(d$kM, d$pop, mean)
# Riverside_2014: 42.6% A, 37.6% C, 19.8% M