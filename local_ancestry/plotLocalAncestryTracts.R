# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2) # melt function
library(gridExtra)
library(viridis)
library(coda) # for hpdi calc
library(emdbook)
library(ggpubr)
library(gtable)
library(grid)
#library(ggExtra) # for marginal histogram ggmarginal
library(sp) # for polygons
library(GISTools)
source("../colors.R") # for color palette

# Script to make plots of simulated and inferred local ancestry in CA and AR, including Fig 5

# load metadata
load("results/meta.RData")

# load local ancestry data
load("results/A.RData")
load("results/C.RData")
load("results/M.RData")
load("results/mean_ancestry_AR_CA.RData")

# mean ancestry across the genome for each pop
admix_proportions <- data.frame(population = pops_by_lat,
                                A = apply(A, 2, mean),
                                C = apply(C, 2, mean),
                                M = apply(M, 2, mean),
                                stringsAsFactors = F)

# ------------ Plot Timing of admixture ---------
# get time of admixture estimates
dir_results <- "results/ancestry_hmm/combined_sept19/posterior"
time_pops <- read.table(paste0(dir_results, "/", "time_pops.txt"), 
                        stringsAsFactors = F, header = F) %>%
  data.table::setnames("population")
admix_times <- do.call(rbind, lapply(c("A", "C", "M"), function(anc) read.table(paste0(dir_results, "/", "time_", anc, ".txt"), 
                                                                                stringsAsFactors = F, header = F) %>%
                                       data.table::setnames(c("ancestry_n", "time", "proportion")) %>%
                                       mutate(ancestry = anc) %>% 
                                       cbind(time_pops, .)))

# plot ancestry times vs. mean ancestry proportion:
t_admix <- admix_times %>%
  left_join(., meta.pop, by = "population") %>%
  filter(., ancestry != "C") %>% # no time, first ancestry
  ggplot(., aes(x = abs(lat), y = time, color = zone, shape = factor(year))) +
  geom_point(alpha = .75) +
  facet_grid(. ~ ancestry) +
  scale_color_manual(values = col_NA_SA_both, name = "Continent") +
  #ggtitle("Inferred time of admixture pulses from HMM") +
  ylab("Time (generations)") +
  xlab("Degrees latitude from the equator") +
  labs(shape = "Sample") +
  theme_classic() +
  scale_shape_manual(values = c(17, 19))
ggsave("plots/time_of_admixture_vs_latitude.png",
       plot = t_admix,
       height = 3, width = 5.2, units = "in")
ggsave("../../bee_manuscript/figures/time_of_admixture_vs_latitude.png",
       plot = t_admix,
       height = 3, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/time_of_admixture_vs_latitude.tiff",
       plot = t_admix,
       height = 3, width = 5.2, units = "in", dpi = 600)

# range of admixture time estimates for 1st pulse (M) and 2nd pulse (A)
admix_times %>%
  group_by(ancestry) %>%
  summarise(mean = mean(time),
            median = median(time),
            min = min(time),
            max = max(time))

# -------- Ancestry covariances (K matrix) ------------- #

# plot K matrix as correlations (note: order by latitude!)
load("results/zAnc.RData")


# # plot the lower triangle of K
# lower_tri <- zAnc_bees$K
# lower_tri[lower.tri(lower_tri, diag = T)] <- NA
# melt(lower_tri) %>% 
#   ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()

k_plot_all <- melt(cov2cor(zAnc_bees$K)) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     name = "Ancestry\ncorrelation") +
  ggtitle("K correlation matrix - bees") +
  xlab("") +
  ylab("")
#k_plot_all

segment_buffer = 0.15
k_plot_fancy <- k_plot_all +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments,
               aes(x = starts + segment_buffer, xend = ends - segment_buffer,
                   y = -0.5, yend = -0.5,
                   color = Region), size = 2, inherit.aes = F) +
  geom_segment(data = NS_segments,
               aes(x = -0.25, xend = -0.25,
                   y = starts + segment_buffer, yend = ends - segment_buffer,
                   color = Region), size = 2, inherit.aes = F) +
  #add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5)) +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5)) +
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), 
                   labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"),
                   labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A)

#k_plot_fancy
ggsave("plots/k_correlation_matrix_all_pops.png", 
       plot = k_plot_fancy,
       height = 3, width = 4, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures/k_correlation_matrix_all_pops.png", 
       plot = k_plot_fancy,
       height = 3, width = 4, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures_main/k_correlation_matrix_all_pops.tif", 
       plot = k_plot_fancy,
       height = 3, width = 4, 
       units = "in", device = "tiff", dpi = 600,
       compression = "lzw", type = "cairo")

# make covariance matrix for the supplement:
k_cov_plot_offdiag <- melt(zAnc_bees$K) %>%
  filter(Var1 != Var2) %>% # omit diagonal
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     option = "viridis",
                     name = "Ancestry\ncovariance") +
  ggtitle("K covariance matrix - bees") +
  xlab("") +
  ylab("") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments,
               aes(x = starts + segment_buffer, xend = ends - segment_buffer,
                   y = -0.5, yend = -0.5,
                   color = Region), size = 2, inherit.aes = F) +
  geom_segment(data = NS_segments,
               aes(x = -0.25, xend = -0.25,
                   y = starts + segment_buffer, yend = ends - segment_buffer,
                   color = Region), size = 2, inherit.aes = F) +
  # add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5), color = "white") +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5), color = "white") +
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A)
#k_cov_plot_offdiag

k_cov_plot_diag <- melt(zAnc_bees$K) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  coord_equal() + # ensures aspect ratio makes a square
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1,
                     option = "inferno",
                     discrete = F,
                     name = "Ancestry\ncovariance") +
  
  ggtitle("K covariance matrix - bees") +
  xlab("") +
  ylab("") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  xlab("Population 1") +
  ylab("Population 2") +
  ggtitle("") +
  # add lines for low A and high A groups
  geom_segment(data = NS_segments,
               aes(x = starts + segment_buffer, xend = ends - segment_buffer,
                   y = -0.5, yend = -0.5,
                   color = Region), size = 2, inherit.aes = F) +
  geom_segment(data = NS_segments,
               aes(x = -0.25, xend = -0.25,
                   y = starts + segment_buffer, yend = ends - segment_buffer,
                   color = Region), size = 2, inherit.aes = F) +
  # add lines marking division between N. and S. American pops
  geom_segment(aes(x = 21.5, xend = 21.5, y = 0.5, yend = 39.5), color = "white") +
  geom_segment(aes(x = 0.5, xend = 39.5, y = 21.5, yend = 21.5), color = "white") +
  scale_x_discrete("Population 1", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_y_discrete("Population 2", breaks = c("CA09","AR14"), labels = c("N. America", "S. America")) +
  scale_color_manual(values = col_low_high_A)
#k_cov_plot_diag
k_cov_plot_fancy <- gridExtra::arrangeGrob(grobs = list(k_cov_plot_diag + guides(color = F),
                                                        k_cov_plot_offdiag + guides(color = F),
                                                        ggpubr::get_legend(k_cov_plot_diag + guides(fill = F) + 
                                                                             theme(legend.position = "bottom"))),
                                           ncol = 2, nrow = 2, 
                                           layout_matrix = rbind(c(1,2),c(3,3)),
                                           heights = c(10,1), widths = c(1,1))
#plot(k_cov_plot_fancy)

ggsave("plots/k_covariance_matrix_all_pops.png", 
       plot = k_cov_plot_fancy,
       height = 3, width = 7.5, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures/k_covariance_matrix_all_pops.png", 
       plot = k_cov_plot_fancy,
       height = 3, width = 7.5, 
       units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/k_covariance_matrix_all_pops.tif", 
       plot = k_cov_plot_fancy,
       height = 3, width = 7.5, 
       units = "in", device = "tiff", dpi = 600,
       compression = "lzw", type = "cairo")



# -------------- Ancestry simulations ----------
load("results/MVNsim_bounded.RData")
load("results/MVNsim_unbounded.RData") # before truncation
load("results/MVNsim_no_cov_bounded.RData")
load("results/PoiBinsim.RData")

p_perc_sim_freq_out_of_bounds = perc_MVNsim_out_of_bounds %>%
  pivot_longer(data = ., cols = c("Lower", "Upper"), names_to = "bound", values_to = "p") %>% 
  left_join(., admix_proportions, by = "population") %>%
  left_join(., meta.pop, by = "population") %>%
  ggplot(., aes(x = A, y = p, color = zone)) +
  geom_point(alpha = .75) +
  facet_wrap(~bound) +
  xlab("Mean A ancestry proportion (ancestry_hmm)") +
  ylab("Percent sims exceeding bound") +
  #ggtitle("Percent simulated population ancestry frequencies < 0 (MVN)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic()
# p_perc_sim_freq_out_of_bounds

# QQ plots function                  
plot_QQ <- function(NA1, NA2, SA1, SA2, axis1, axis2, ps_qq = seq(0, 1, length.out = 10000)){
  d_qq <- bind_rows(data.frame(p = ps_qq,
                               zone = "N. America", stringsAsFactors = F) %>%
                      mutate(q1 = quantile(NA1, p), 
                             q2 = quantile(NA2, p)),
                    data.frame(p = ps_qq,
                               zone = "S. America", stringsAsFactors = F) %>%
                      mutate(q1 = quantile(SA1, p), 
                             q2 = quantile(SA2, p)))
  ggplot(data = d_qq) + 
    geom_point(aes(x = q1, y = q2, color = p)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab(axis1) +
    ylab(axis2) +
    scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
    theme_classic() +
    facet_grid(.~zone)
}
# make QQ plots for paper:
# effect of truncating at [0,1] is small:
MVN_vs_01_truncated_MVN_qq = plot_QQ(NA1 = meanA_MVNsim_CA_bounded, NA2 = meanA_MVNsim_CA_unbounded, 
                                     SA1 = meanA_MVNsim_AR_bounded, SA2 = meanA_MVNsim_AR_unbounded, 
                                     axis1 = "MVN simulation", axis2 = "MVN simulation unbounded", 
                                     ps_qq = seq(0, 1, length.out = 10000))

effect_truncation_on_MVN <- arrangeGrob(p_perc_sim_freq_out_of_bounds + ggtitle("A"),
                                        MVN_vs_01_truncated_MVN_qq + ggtitle("B"),
                                        nrow = 2,
                                        ncol = 1)
#plot(effect_truncation_on_MVN)
ggsave("../../bee_manuscript/figures/effect_truncation_on_MVN.png", 
       plot = effect_truncation_on_MVN,
       height = 6, width = 5.2, dpi = 600,
       units = "in", device = "png")
ggsave("../../bee_manuscript/figures_supp/effect_truncation_on_MVN.tiff", 
       plot = effect_truncation_on_MVN,
       height = 6, width = 5.2, dpi = 600,
       units = "in", device = "tiff")

# plot fit of data to different simulations:
sim_NA = list(meanA_MVNsim_CA_bounded, meanA_MVNsim_CA_no_cov_bounded, PoiBinsim_CA)
sim_SA = list(meanA_MVNsim_AR_bounded, meanA_MVNsim_AR_no_cov_bounded, PoiBinsim_AR)
sim_names = c("MVN", "MVN zero covariances", "Poisson binomial")
ps_qq = seq(0, 1, length.out = 10000)

d_qq <- do.call(rbind, lapply(1:3, function(i)
  bind_rows(data.frame(p = ps_qq,
                       zone = "N. America",
                       sim = sim_names[i],
                       stringsAsFactors = F) %>%
              mutate(q2 = quantile(meanA_CA, p),
                     q1 = quantile(sim_NA[[i]], p)),
            data.frame(p = ps_qq,
                       zone = "S. America",
                       sim = sim_names[i],
                       stringsAsFactors = F) %>%
              mutate(q2 = quantile(meanA_AR, p), 
                     q1 = quantile(sim_SA[[i]], p)))))
ggplot(data = d_qq) + 
  geom_point(aes(x = q1, y = q2, color = p), size = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Simulated A frequencies") +
  ylab("Observed A frequencies") +
  scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
  theme_classic() +
  facet_grid(zone~sim)#, scales = "free_x")
ggsave("../../bee_manuscript/figures/qq_vs_data_sim_comparison.png",
       width = 7.5, height = 6, 
       units = "in", dpi = 600, device = "png")
ggsave("../../bee_manuscript/figures_supp/qq_vs_data_sim_comparison.tiff",
       width = 7.5, height = 6, 
       units = "in", dpi = 600, device = "tiff")
ggsave("plots/qq_vs_data_sim_comparison.png",
       width = 7.5, height = 6, 
       units = "in", dpi = 600, device = "png")


# Plot MVN 2D density compared to data
# Compute kde2d
kd_data <- with(data = data.frame(x = meanA_CA, y = meanA_AR), 
                MASS::kde2d(y = y, x = x, n = 100,
                            lims = c(c(0,.6),c(0,.6))))
kd_mvn <- with(data = data.frame(x = meanA_MVNsim_CA_unbounded, 
                                 y = meanA_MVNsim_AR_unbounded), 
               MASS::kde2d(y = y, x = x, n = 100,
                           lims = c(c(0,.6),c(0,.6))))

# function to plot
plot_2d_density <- function(){
  contour(kd_mvn, col = viridis(10)[1], nlevels = 7, lwd = 2,
          xlim = c(0.12, 0.35),
          ylim = c(0.32, 0.5),
          xlab = "N. America",
          ylab = "S. America")
  contour(kd_data, add = T, col = viridis(10)[9], cex = 2, nlevels = 7)
  legend("bottomright", legend = c("Observed A frequency", "Simulated A frequency (MVN)"),
         lwd = 1, lty = 1, col = viridis(10)[c(9,1)], cex = 0.5)
}
# save plots
png(file = "../../bee_manuscript/figures/comparison_2d_density_data_mvn.png",
    height = 5, width = 5.2, res = 600, units = "in")
plot_2d_density()
dev.off()
tiff(file = "../../bee_manuscript/figures_supp/comparison_2d_density_data_mvn.tiff",
     height = 5, width = 5.2, res = 600, units = "in")
plot_2d_density()
dev.off()
png(file = "plots/comparison_2d_density_data_mvn.png",
    height = 5, width = 5.2, res = 600, units = "in")
plot_2d_density()
dev.off()


# What are the false discovery rates?
FDRs <- read.table("results/FDRs_MVN_01_high_low_A_percent_cutoffs.txt", 
                   stringsAsFactors = F, header = T, sep = "\t")



# Calculate 99% (and 90%, 75%) HPDI for MVN simulations
b <- coda::as.mcmc(cbind(meanA_MVNsim_CA_bounded, meanA_MVNsim_AR_bounded))
a <- emdbook::HPDregionplot(b, vars = c("meanA_MVNsim_CA_bounded", "meanA_MVNsim_AR_bounded"), 
                            prob = c(0.99, .9, .75), 
                            n = 100, # number of grid points
                            lwd = 3, 
                            h = .05,
                            add = FALSE)

# HPDI for MVN without covariances
b_MVNnoCov <- coda::as.mcmc(cbind(meanA_MVNsim_CA_no_cov_bounded, meanA_MVNsim_AR_no_cov_bounded))
a_MVNnoCov <- emdbook::HPDregionplot(b_MVNnoCov, vars = c("meanA_MVNsim_CA_no_cov_bounded", "meanA_MVNsim_AR_no_cov_bounded"),
                            prob = c(0.99, .9, .75),
                            n = 100,
                            lwd = 3,
                            h = .05,
                            add = FALSE)
# HPDI for poisson binomial
b_poiBin <- coda::as.mcmc(cbind(PoiBinsim_CA, PoiBinsim_AR))
a_poiBin <- emdbook::HPDregionplot(b_poiBin, vars = c("PoiBinsim_CA", "PoiBinsim_AR"),
                                     prob = c(0.99, .9, .75),
                                     n = 100,
                                     lwd = 3,
                                     h = .05,
                                     add = FALSE)


p_density2 <- ggplot() +
  geom_hex(data = data.frame(meanA_CA = meanA_CA, meanA_AR = meanA_AR), 
           aes(x = meanA_CA, y = meanA_AR, color = ..count.., fill = ..count..), bins = 200) +
  scale_fill_viridis(option = "magma", end = 1, begin = 0.1, name = "SNP count") +
  scale_color_viridis(option = "magma", end = 1, begin = 0.1, name = "SNP count") +
  theme_light() +
  coord_fixed() +
  theme(panel.grid.minor = element_blank()) +
  xlab("A ancestry in North America") +
  ylab("A ancestry in South America") +
  geom_polygon(data = data.frame(a[[1]]), 
               aes(x = x, y = y),
               fill = NA, color = "orange", lwd = 1)
plot(p_density2)

ggsave("plots/A_CA_vs_AR_density.png",
       plot = p_density2,
       height = 4.3, width = 5.2, 
       units = "in", dpi = 600,
       device = "png")

# marginal histograms to compare simulations with data
sim_compare_CA <- 
  data.frame(poisson_binomial = PoiBinsim_CA,
             MVN_no_covariance = meanA_MVNsim_CA_no_cov_bounded,
             MVN_with_covariance = meanA_MVNsim_CA_bounded) %>%
  tidyr::gather(., "distribution", "A") %>%
  bind_rows(., data.frame(distribution = "observed_data", A = meanA_CA, stringsAsFactors = F)) %>%
  mutate(distribution = factor(distribution, levels = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial", "observed_data"), ordered = T))

p_sim_compare_CA_with_legend <- sim_compare_CA %>%
  filter(., distribution != "observed_data") %>%
  ggplot(., aes(x = A, color = distribution)) +
  geom_histogram(data = filter(sim_compare_CA, distribution == "observed_data"),
                 aes(y = ..density..), bins = 40, fill = "white", color = "black", lwd = 0.3) +
  geom_line(stat = "density", alpha = 0.75, lwd = 1, aes(linetype = distribution)) +
  theme_classic() +
  xlab("Mean African ancestry") +
  ylab("Density") +
  scale_color_manual(values = c("orange", viridis(4)[2:3]), name = "Model", 
                     limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                     labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_linetype_manual(values=c(1, 2, 3), name = "Model",
                        limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     labels = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     limits = c(0.05, 0.6)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "in"),
        panel.border = element_blank(), 
        panel.spacing = unit(0, "in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#p_sim_compare_CA_with_legend
p_sim_compare_CA <- p_sim_compare_CA_with_legend + guides(color = F, linetype = F)

# marginal histogram for A frequencies in Argentina
sim_compare_AR <- 
  data.frame(poisson_binomial = PoiBinsim_AR,
             MVN_no_covariance = meanA_MVNsim_AR_no_cov_bounded,
             MVN_with_covariance = meanA_MVNsim_AR_bounded) %>%
  tidyr::gather(., "distribution", "A") %>%
  bind_rows(., data.frame(distribution = "observed_data", A = meanA_AR, stringsAsFactors = F)) %>%
  mutate(distribution = factor(distribution, levels = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial", "observed_data"), ordered = T))

p_sim_compare_AR <- sim_compare_AR %>%
  filter(., distribution != "observed_data") %>%
  ggplot(., aes(x = A, color = distribution)) +
  geom_histogram(data = filter(sim_compare_AR, distribution == "observed_data"),
                 aes(y = ..density..), bins = 40, fill = "white", color = "black", lwd = 0.3) +
  geom_line(stat = "density", alpha = 0.75, lwd = 1, aes(linetype = distribution)) +
  theme_classic() +
  xlab("Mean African ancestry") +
  ylab("Density") +
  scale_color_manual(values = c("orange", viridis(4)[2:3]), name = NULL, 
                     limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                     labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  scale_linetype_manual(values=c(1, 2, 3), name = "Model",
                        limits = c("MVN_with_covariance", "MVN_no_covariance", "poisson_binomial"),
                        labels = c("poisson_binomial"="Poisson Binomial", "MVN_no_covariance"="MVN variance only", "MVN_with_covariance"="MVN")) +
  guides(color = F, linetype = F) +
  xlim(c(0.15, 0.7)) +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        axis.line.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "in"),
        panel.border = element_blank(), 
        panel.spacing = unit(0, "in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#p_sim_compare_AR

# construct Fig5: 2D histogram + marginal histograms
# extract legends
legend_p_sim_compare_bottom <- ggpubr::get_legend(p_sim_compare_CA_with_legend +
                                            theme(plot.margin = unit(c(1,1,1,1), "lines"),
                                                  legend.box.margin = unit(c(1,1,1,1), "lines"),
                                                  legend.position = "bottom"))
legend_p_sim_compare <- ggpubr::get_legend(p_sim_compare_CA_with_legend +
                                     theme(plot.margin = unit(c(1,1,1,1), "lines"),
                                           legend.box.margin = unit(c(1,1,1,1), "lines")))
legend_p_density3 <- get_legend(p_density2 +
                                  theme(legend.position = "right",
                                        plot.margin = unit(c(1,1,1,1), "lines"),
                                        legend.box.margin = unit(c(1,1,1,1), "lines")))

p_density3 <- p_density2 + 
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0, unit='cm'),
        legend.title = element_text(size = 10)
  ) +
  guides(fill = guide_colorbar(title.vjust = 0.75)) +
  ylim(c(0.15, 0.7)) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     labels = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     limits = c(0.05, 0.6))
#plot(p_density3)
g <- ggplotGrob(p_density3)
panel_id <- g$layout[g$layout$name == "panel", c("t","l")]
#g <- g0
g <- gtable_add_cols(g, unit(1.5, "in"))
g <- gtable_add_grob(g, ggplotGrob(p_sim_compare_AR + 
                                     theme(axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
                                           axis.line = element_line(colour = "grey87", size = rel(0.5)))),
                     t = panel_id$t, l = ncol(g))
g <- gtable_add_rows(g, unit(1.5, "in"), 0) # add row on top
g <- gtable_add_grob(g, ggplotGrob(p_sim_compare_CA + 
                                     theme(axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
                                           axis.line = element_line(colour = "grey87", size = rel(0.5)))),
                     t = 1, l = panel_id$l)
g <- gtable_add_rows(g, unit(0.4, "in"), -1)
g <- gtable_add_grob(g, legend_p_sim_compare_bottom$grobs[[1]],
                     t = -1, l = panel_id$l) # add row to bottom?

ggsave("plots/A_CA_vs_AR_density_marginal_hist.png",
       plot = g,
       height = 5.4, width = 5.2, 
       units = "in", dpi = 600,
       device = "png")
ggsave("../../bee_manuscript/figures/A_CA_vs_AR_density_marginal_hist.png",
       plot = g,
       height = 5.4, width = 5.2, 
       #height = 6.4, width = 6, 
       units = "in", dpi = 600, 
       device = "png")
ggsave("../../bee_manuscript/figures_main/A_CA_vs_AR_density_marginal_hist.tif",
       plot = g,
       height = 5.4, width = 5.2, 
       units = "in", dpi = 600, 
       device = "tiff",
       compression = "lzw", type = "cairo")

# what % of data points fall within approximate 99% HPDI for 
# each simulated model?
poly_MVN <- Polygons(srl = list(Polygon(coords = with(a[[1]], cbind(x, y)), hole = F)),
                     ID = "MVN")
poly_MVNNoCov <- Polygons(srl = list(Polygon(coords = with(a_MVNnoCov[[1]], cbind(x, y)), hole = F)),
                          ID = "NoCov")
poly_PoiBin <- Polygons(srl = list(Polygon(coords = with(a_poiBin[[1]], cbind(x, y)), hole = F)),
                        ID = "Poibin")
polygons <- SpatialPolygons(Srl = list(poly_MVN, poly_MVNNoCov, poly_PoiBin))
plot(polygons)
some_points <- A_AR_CA %>%
  dplyr::select(CA, AR) %>%
  SpatialPoints(coords = .)
# % of points. note: it's approximate in part because polygons are approximately defined with finite # of points
# and also because I use the sample of 100k points from these distributions to approximate the density rather than the true density.
GISTools::poly.counts(pts = some_points, polys = polygons)/nrow(A_AR_CA)
