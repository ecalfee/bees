library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggpubr)
source("../colors.R")
#ancestry = "A"
#beagle_prefix = "all_CA_AR_ACM"
#beagle_prefix = "A_CA_AR"
# get bee list (must match order in beagle.gz file that went to PCAngsd!)
IDs <- read.table(paste0("../bee_samples_listed/combined_sept19.list"), stringsAsFactors = F,
                  header = F)
colnames(IDs) <- c("Bee_ID")
# get bee metadata
meta <- read.table("../bee_samples_listed/all.meta", stringsAsFactors = F, 
                   header = T, sep = "\t")
#b <- left_join(IDs, meta, by = "Bee_ID")
#which(b$group == "C")-1
#which(b$group == "M")-1
#which(b$group == "A")-1
#which(!(b$group %in% c("A", "C", "M")))-1
# function to get pca in dataframe for bees included and focal ancestry
get_pca <- function(ancestry, beagle_prefix, ref_incl){ # included reference bee pops can be a vector, e.g ("A", "M", "C")
  # get covariance matrix from PCAngsd output
  cov_data = read.table(paste0("results/combined_sept19/", ancestry, 
                               "/PCA/", beagle_prefix, ".cov"),
                        header = F, stringsAsFactors = F)
  
  # join bee IDs and metadata
  bees <- dplyr::left_join(IDs, meta, by = "Bee_ID") %>%
    filter(group %in% c(unlist(ref_incl), "S_CA", "N_CA", "CA_2018", "AR_2018"))
  
  # PCA from covariance matrix:
  pca <- eigen(cov_data) # take PCA of covariance matrix
  m = 10 # make small dataframe w/ only first m eigenvectors
  pca_small <- data.frame(pca$vectors[ , 1:m])
  colnames(pca_small) = paste0("PC", 1:m)
  
  # join bams and firstPCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
  d <- bind_cols(bees, pca_small)  %>%
    arrange(., population) %>%
    arrange(., strain) %>%
    mutate(., label = ifelse(group %in% c("A", "C", "M"), 
                             group,
                             ifelse(group == "AR_2018", "S. America", "N. America")),
           geo_label = ifelse(label %in% c("S. America", "N. America"),
                              label,
                              ifelse(geographic_location_short == "Pag Island, Croatia",
                                     "Croatia", 
                                     ifelse(geographic_location_short == "South Africa",
                                            "S. Africa", geographic_location_short))),
           geo_label = factor(geo_label, ordered = T, levels = names(col_geo_label)),
           beagle_prefix = beagle_prefix,
           ancestry_4pca = ancestry)
  
  # rounded eigen values
  PC_var_explained = 100*pca$values/sum(pca$values)
  return(list(d = d, PC_var_explained = PC_var_explained))
}

ancestry_list <- c("A", "C", "M", "A", "C", "M")
beagle_prefix_list <- c("A_CA_AR", "C_CA_AR", "M_CA_AR",
                        rep("all_CA_AR_ACM", 3))
ref_incl_list <- list("A", "C", "M", 
                   c("A", "C", "M"), c("A", "C", "M"), c("A", "C", "M"))
pcas <- #lapply(1:length(ancestry_list), function(i)
  lapply(1:length(ancestry_list), function(i)
               get_pca(ancestry = ancestry_list[i], 
                       beagle_prefix = beagle_prefix_list[i], 
                       ref_incl = ref_incl_list[[i]]))

# plot first PC's
# PC1 and 2, all no filter
pca_plots <- lapply(4:length(pcas), function(i){
  pcas[[i]]$d %>%
    ggplot(., aes(PC1, PC2, 
                  color = abs(lat),
                  shape = label)) + 
    geom_point(alpha = 0.5) +
    xlab(paste0("PC1 (", round(pcas[[i]]$PC_var_explained[1], 2), "%)")) +
    ylab(paste0("PC2 (", round(pcas[[i]]$PC_var_explained[2], 2), "%)")) +
    theme_classic() +
    ggtitle(ancestry_list[[i]]) +
    scale_color_viridis_c(name = "Lat") +
    scale_shape_manual(values = shapes_ACM_NA_SA,
                       name = "Source")
})
arrangeGrob(grobs = pca_plots, nrow = 3) %>%
  plot()
legend1 <- ggpubr::get_legend(pca_plots[[1]])
plot(legend1)
pca_acm_by_anc <- arrangeGrob(grobs = list(pca_plots[[1]] + theme(legend.position = "none"), 
                         pca_plots[[2]] + theme(legend.position = "none"),
                         pca_plots[[3]] + theme(legend.position = "none"),
                         legend1),
            nrow = 3, ncol = 2, 
            layout_matrix = rbind(c(1,4),c(2,4),c(3,4)),
            heights = c(1,1,1), widths = c(4,1))
plot(pca_acm_by_anc)
ggsave("../../bee_manuscript/figures/pca_acm_by_anc.png",
       plot = pca_acm_by_anc,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/pca_acm_by_anc.tif",
       plot = pca_acm_by_anc,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, 
       device = "tiff",
       compression = "lzw", type = "cairo")


# instead plot PCA colored by continent
pca_by_anc <- lapply(1:3, function(i){
  pcas[[i]]$d %>%
    ggplot(., aes(PC1, PC2,
                  color = geo_label)) + 
    geom_point(alpha = 0.75) +
    xlab(paste0("PC1 (", round(pcas[[i]]$PC_var_explained[1], 2), "%)")) +
    ylab(paste0("PC2 (", round(pcas[[i]]$PC_var_explained[2], 2), "%)")) +
    theme_classic() +
    ggtitle(ancestry_list[[i]]) +
    scale_color_manual(values = col_geo_label, 
                       name = element_blank())
})

p_pca_by_anc <- arrangeGrob(grobs = pca_by_anc, nrow = 3)
plot(p_pca_by_anc)
ggsave("../../bee_manuscript/figures/p_pca_by_anc.png",
       plot = p_pca_by_anc,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/p_pca_by_anc.tif",
       plot = p_pca_by_anc,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, 
       device = "tiff",
       compression = "lzw", type = "cairo")

# also plot pc1 by latitude
pc1_by_lat <- lapply(1:3, function(i){
  pcas[[i]]$d %>%
    mutate(lat = ifelse(is.na(lat), 39.3, lat)) %>% # put ref bees @ just outside range of lat
    ggplot(., aes(abs(lat), PC1,
                  color = geo_label)) + 
    geom_point(alpha = 0.75) +
    xlab("Absolute latitude (degrees)") +
    ylab(paste0("PC1 (", round(pcas[[i]]$PC_var_explained[1], 2), "%)")) +
    theme_classic() +
    ggtitle(ancestry_list[[i]]) +
    scale_color_manual(values = col_geo_label, 
                       name = element_blank()) +
    #geom_hline(yintercept = 0, linetype = "dashed")
    geom_segment(x = 28, xend = 38.8, y = 0, yend = 0, 
                     color = "black", linetype = "dashed",
                 size = 0.01)
})

pc1_anc_by_lat <- arrangeGrob(grobs = pc1_by_lat, nrow = 3)
plot(pc1_anc_by_lat)
ggsave("../../bee_manuscript/figures/pc1_anc_by_lat.png",
       plot = pc1_anc_by_lat,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/pc1_anc_by_lat.tif",
       plot = pc1_anc_by_lat,
       width = 5.2, height = 5, 
       units = "in", dpi = 600, 
       device = "tiff",
       compression = "lzw", type = "cairo")



# just get covariance matrix to plot it:
get_cov <- function(ancestry, beagle_prefix, names){
  # get covariance matrix from PCAngsd output
  cov_data = read.table(paste0("results/combined_sept19/", ancestry, 
                               "/PCA/", beagle_prefix, ".cov"),
                        header = F, stringsAsFactors = F) %>%
    as.matrix(.)
  # need to name rows and columns of matrix using bee IDs
  colnames(cov_data) = names
  rownames(cov_data) = names
  return(cov_data)
  }

bees = dplyr::left_join(IDs, meta, by = "Bee_ID") %>%
  filter(group %in% c("S_CA", "N_CA", "CA_2018", "AR_2018"))
covs <- lapply(ACM, function(a) get_cov(ancestry = a, 
                                        beagle_prefix = "CA_AR", 
                                        names = bees$Bee_ID))

melt(covs[[1]]) %>%
  mutate(Var1 = factor(Var1,
                       levels = arrange(bees, lat)$Bee_ID, 
                       ordered = T),
         Var2 = factor(Var2,
                       levels = arrange(bees, lat)$Bee_ID, 
                       ordered = T)) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     limits = c(-.12, .32),
                     name = "Ancestry\ncovariance") +
  ggtitle("covariances within A") +
  xlab("") +
  ylab("")

cov_plots <- lapply(1:3, function(i) {
  melt(cov2cor(covs[[i]])) %>%
    mutate(Var1 = factor(Var1,
                         levels = arrange(bees, lat)$Bee_ID, 
                         ordered = T),
           Var2 = factor(Var2,
                         levels = arrange(bees, lat)$Bee_ID, 
                         ordered = T)) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     limits = c(-.12, .32),
                     name = "Ancestry\ncorrelation") +
  ggtitle(ACM[i]) +
  xlab("") +
  ylab("")
})
plot(cov_plots[[1]])
plot(cov_plots[[2]])
plot(cov_plots[[3]])
