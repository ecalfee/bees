library(tidyr)
library(dplyr)
library(ggplot2)
library(MASS) # for linear discriminant analysis
library(geomorph)
library(abind) # for combining tps files

# corn says that I need the absolute scaling to be applied before procrustes analysis

# for linux instructions to get rgl dependency:
# https://stackoverflow.com/questions/31820865/error-in-installing-rgl-package
#For mac you need xquartz I think.

# resource: http://www.indiana.edu/~g562/PBDB2013/Day%202B%20-%20Geometric%20Morphometrics%20in%20R.pdf

wings <- readland.tps("wings.tps", specID = "imageID")
# images labeled lineage-subspecies-colony-individual,
# e.g. A-scu-0160-01.dw.png
imageIDs <- dimnames(wings)[[3]]
meta <- data.frame(imageID = tools::file_path_sans_ext(tools::file_path_sans_ext(imageIDs))) %>%
  tidyr::separate(., imageID, c("lineage", "subspecies", "colony", "individual"))

plotAllSpecimens(wings)

wings.procrustes <- gpagen(wings)
plot(wings.procrustes)
# make a tidy data frame with the procrustes results and metadata for each set of points
d_tidy <- data.frame(meta, wings.procrustes$data) %>%
  tidyr::gather(., "label", "coord", 5:42) %>%
  tidyr::separate(., label, c("ignore", "landmark", "axis")) %>%
  dplyr::select(., -ignore) %>%
  spread(., axis, coord)

# plot all landmarks for all individuals
d_tidy %>% 
  ggplot(., aes(x = X, y = Y, color = lineage)) +
  geom_point()

# plot means per colony to visualize group separation (subtle)
d_tidy %>% 
  group_by(lineage, colony, landmark) %>%
  dplyr::summarize(mean_X = mean(X), mean_Y = mean(Y)) %>%
  ggplot(., aes(x = mean_X, y = mean_Y, color = lineage)) +
  geom_point(alpha = .5) +
  ggtitle("A/C/M ref. bees, colony mean wing landmarks")

# example plot to look at an individual landmark
d_tidy %>% 
  filter(landmark == 1) %>%
  ggplot(., aes(x = X, y = Y, color = colony, shape = lineage)) +
  geom_point(alpha = .5) +
  ggtitle("A/C/M ref. bees, ind bees landmark 1")

# plot a pca
pca.wings <- plotTangentSpace(wings.procrustes$coords, 
                              axis1 = 1, axis2 = 2,
                              groups = as.factor(meta$lineage))
unique(meta$lineage) # A is black, C is red, M is green
pca.wings$pc.summary
# hmm, PC1 separates C from A/M, but none of the top PCs appear to separate out A
# The identifly paper used Canonical variate analysis, otherwise known as
# discriminant function analysis
# instead of a PCA

# linear discriminant analysis for reference populations:
d_wide <- data.frame(meta, wings.procrustes$data)
lda.wings <- lda(meta$lineage ~ .,
                 prior = c(1/3, 1/3, 1/3), # uniform prior on the 3 groups
                 data = dplyr::select(wings.procrustes$data, -Csize))
#Warning: In lda.default(x, grouping, ...) : variables are collinear
lda.wings

# make arrows - 
# function taken from DataCamp tutorial: https://campus.datacamp.com/courses/helsinki-open-data-science/clustering-and-classification?ex=7 
lda.arrows <- function(x, myscale = 1, arrow_heads = 0.1, color = "black", 
                         tex = 0.75, choices = c(1,2)){
  heads <- coef(x)
  arrows(x0 = 0, y0 = 0, 
         x1 = myscale * heads[,choices[1]], 
         y1 = myscale * heads[,choices[2]], col=color, length = arrow_heads)
  text(myscale * heads[,choices], labels = row.names(heads), 
       cex = tex, col=color, pos=3)
}

# plot
png("plots/lda.png")
plot(lda.wings, main = "Linear discriminant analysis of A/M/C groups", 
     col = rainbow(3)[as.factor(meta$lineage)])
lda.arrows(lda.wings, myscale = .02)
dev.off()

plot(lda.wings, main = "Linear discriminant analysis of A/M/C groups", 
     col = rainbow(56)[as.factor(meta$colony)])

# how well does the function discriminate individual bees?
# fit 2 with jackknife leave one out
lda.wings.pred <- lda(meta$lineage ~ .,
                 CV=TRUE,
                 prior = c(1/3, 1/3, 1/3), # uniform prior on the 3 groups
                 data = dplyr::select(wings.procrustes$data, -Csize))
# classification table -- very high percent correct
classified.wings <- table(meta$lineage, lda.wings.pred$class)
diag(prop.table(classified.wings, 1))

dimnames(wings.procrustes$coords)[[3]]<-d_wide$lineage

mean.ACM <- lapply(c("A", "C", "M"), function(l) 
  mshape(wings.procrustes$coords[ , , d_wide$lineage == l]))
par(mfrow=c(3,1))
lapply(mean.ACM, plot)
par(mfrow=c(1,1))
# lookup: redundancy analysis
# Chris Klingenberg fly wing morphometrics
# plot ref to target compares two shapes to each other
# maybe split into 2 groups - A vs. non-A

# test
# 5 test hives - 2 ind's per hive
test <- readland.tps("5Hives.tps", specID = "imageID")
# images labeled lineage-subspecies-colony-individual,
# e.g. A-scu-0160-01.dw.png
testIDs <- dimnames(test)[[3]]
test_meta <- data.frame(imageID = substr(testIDs, 9, 100)) %>%
  tidyr::separate(., imageID, c("colony", "individual", "subspeciesN"))
plotAllSpecimens(test)

test.procrustes <- gpagen(test)
par(mfrow=c(1,1))
plot(test.procrustes) # right wing orientation
plot(wings.procrustes) # left wing orientation
# ok, it's easy to flip the coordinates
test.procrustes2 <- gpagen(rotate.coords(test, type = c("flipX")))
plot(test.procrustes2)
# flip the original wings data and combine with the test data
combined_data <- abind(rotate.coords(wings, type = c("flipX")), test, along = 3)
dimnames(combined_data)[[3]] <- c(dimnames(wings)[[3]], dimnames(test)[[3]])

# put metadata together
test_meta2 <- test_meta %>%
  mutate(., indN = as.integer(individual)) %>%
  mutate(., source = "Petra")
meta2 <- meta %>%
  mutate(., indN = as.integer(individual)) %>%
  mutate(., source = "Museum")
combined_meta <- bind_rows(meta2, test_meta2)

combined.procrustes <- gpagen(combined_data)
plot(combined.procrustes)

# make into dataframe
combined.tidy <- data.frame(combined_meta, combined.procrustes$data) %>%
  tidyr::gather(., "label", "coord", 8:45) %>%
  tidyr::separate(., label, c("ignore", "landmark", "axis")) %>%
  dplyr::select(., -ignore) %>%
  spread(., axis, coord)

combined.tidy %>%
  ggplot(aes(x = X, y = Y, color = source)) +
  geom_point()
# which bees are in both datasets?
in.test <- paste(test_meta$colony, test_meta$individual, sep = ".")
duplicates.tidy <- combined.tidy %>%
  mutate(., Bee_ID = paste(as.integer(colony), indN, sep = ".")) %>%
  filter(., Bee_ID %in% in.test) 

for (i in 1:19){
  p <- duplicates.tidy %>%
    filter(., landmark == as.character(i)) %>%
    ggplot(aes(x = X, y = Y, shape = source, color = Bee_ID)) +
    geom_point() +
    ggtitle(paste("landmark ", i))
  plot(p)
}

# predict using lda some of the museum data
predict(lda.wings, 
        newdata = dplyr::select(wings.procrustes$data, -Csize)[1:5, ])  
# now predict using new 'test' data
predict(lda.wings, 
        newdata = dplyr::select(test.procrustes2$data, 
                                -Csize))$class  
test_meta$subspeciesN


# get triangle distance measurements between bees


