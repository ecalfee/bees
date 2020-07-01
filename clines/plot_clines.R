# script to plot clines in ancestry, e.g. with latitude and environment
library(sp)
library(gridExtra)
library(raster)
library(bedr)
library(geosphere)
#library(rethinking)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(purrr)
library(nlstools)
library(viridisLite)
#library(betareg) # alternative ML fitting
library(xtable)
library(nls.multstart)
source("../local_ancestry/calc_FDRs.R")
source("../colors.R") # get color palette
source("cline_functions.R") # loads logistic and stepped clines etc.
load("results/d_A.RData") # load data
load("../wing_analysis/results/wing_fits.RData") # wing linear model fits
load("../local_ancestry/results/meta.RData")
load("../local_ancestry/results/mean_ancestry_AR_CA.RData") # to tag ancestry outlier SNPs

# plot some climate variables:
d_A %>%
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "MinTempColdestMonth", "AnnualPrecip")) %>%
  ggplot(., aes(x = abs(lat), y = value, color = continent, shape = year)) +  
                #color = geographic_location_short)) + 
  geom_point(size = 1, alpha = .75) + 
  facet_wrap(~climate_var, scales = "free_y", ncol = 1) +
  xlab("Degrees latitude from the equator") + 
  theme(legend.position="bottom") + 
  labs(color = "") +
  ylab("") +
  scale_shape_manual(values = c(17, 19), name = "Sample") +
  scale_color_manual(values = col_NA_SA_both, name = "Continent") +
  theme_classic()
ggsave("plots/climate_variables_across_latitude.png", 
       height = 5, width = 5, units = "in")
# save in figures for manuscript:
ggsave("../../bee_manuscript/figures/climate_variables_across_latitude.png", 
       height = 5.2, width = 5, units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/climate_variables_across_latitude.tiff", 
       height = 5.2, width = 5, units = "in", device = "tiff", dpi = 600)
# outliers min and max:
d_A %>%
  filter(continent == "N. America") %>%
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "MinTempColdestMonth", "AnnualPrecip")) %>%
  group_by(climate_var) %>%
  summarise(population = d_A$population[d_A$continent == "N. America"][which.min(value)],
            lat = d_A$lat[d_A$continent == "N. America"][which.min(value)],
            min = min(value))
d_A %>%
  filter(continent == "N. America") %>%
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "MinTempColdestMonth", "AnnualPrecip")) %>%
  group_by(climate_var) %>%
  summarise(population = d_A$population[d_A$continent == "N. America"][which.max(value)],
            lat = d_A$lat[d_A$continent == "N. America"][which.max(value)],
            max = max(value))
d_A %>%
  filter(population == "Riverside_2014") %>%
  tidyr::gather(., "climate_var", "value", c("AnnualMeanTemp", "MeanTempColdestQuarter", "MinTempColdestMonth", "AnnualPrecip")) %>%
  group_by(climate_var, lat) %>%
  summarise(value = mean(value))
# outliers:
# we have a cold/wet mountain top at 33.7 degrees lat
# and a dry/hot desert at 33.8 degrees lat
# in the Riverside 2014 dataset

# how many years for spread?
km_per_degree_lat = 111 # 111.699 km at the poles, 110.567 estimate from the equator
years_tot = 2018-1957
#1957 to 2018 is approx 30 to 60 generations depending on whether it's 1 gen every year or every other year
(2018-1970)/2 # ~24 generations since reaching current cline center in Argentina

# Are known feral bees different?
# enjambre = collected within 5m of a known feral nest
glm.feral = glm(alpha ~ abs_lat + enjambre, data = d_A, family = gaussian(link = "logit"))
summary(glm.feral)

#----------some physical distances---------------#
# cline distances:
#Rio Claro, Sao Paulo Brazil: 22.4149° S, 47.5651° W (Google maps)
sao_paulo <- data.frame(long = -47.5651, lat = -22.4149)

# approximate length of sample transects
dist_NA <- matrix(nrow = nrow(filter(meta.pop, zone == "N. America")),
                  ncol = nrow(filter(meta.pop, zone == "N. America")))
dist_SA <- matrix(nrow = nrow(filter(meta.pop, zone == "S. America")),
                  ncol = nrow(filter(meta.pop, zone == "S. America")))
# calculate distances between population pairs within N. American cline
for (i in 1:nrow(filter(meta.pop, zone == "N. America"))){
  for (j in 1:nrow(filter(meta.pop, zone == "N. America"))){
    dist_NA[i, j] <- distm(filter(meta.pop, zone == "N. America")[i, c("long", "lat")], 
                           filter(meta.pop, zone == "N. America")[j, c("long", "lat")], 
                           fun = distGeo)/1000
  }
}
# and within S. American cline
for (i in 1:nrow(filter(meta.pop, zone == "S. America"))){
  for (j in 1:nrow(filter(meta.pop, zone == "S. America"))){
    dist_SA[i, j] <- distm(filter(meta.pop, zone == "S. America")[i, c("long", "lat")], 
                           filter(meta.pop, zone == "S. America")[j, c("long", "lat")], 
                           fun = distGeo)/1000
  }
}
# how much distance does sampling in these clines span?
max(dist_NA) # 772.5385 - note: CA has more W<->E difference between sampling endpoints
max(dist_SA) # 886.3263

# A more appropriate comparison is just the N<->S length:
group_by(d_A, continent) %>%
  summarise(min_lat = min(lat),
            max_lat = max(lat),
            mean_long = mean(long),
            length_km = distm(c(mean_long, min_lat),
                              c(mean_long, max_lat),
                              fun = distGeo)/1000)

# distance from sao paulo to each hybrid zone
d_A %>%
  group_by(., continent) %>%
  summarise(mean = mean(km_from_sao_paulo),
            min = min(km_from_sao_paulo),
            max = max(km_from_sao_paulo))
#continent   mean   min    max
#N. America 9868. 9561. 10268.
#S. America 1680. 1342.  2018.

# distance in km between hottest and coldest site (CA mountains and desert)
distm(d_A[which.max(d_A$AnnualMeanTemp), c("long", "lat")], 
      d_A[which.min(d_A$AnnualMeanTemp), c("long", "lat")], 
      fun = distGeo)/1000


############-----logistic clines w/ nls()------#########

# fit clines, rescale by 0.84 because that's the maximum A ancestry (proportion in Brazil):
fit_lat_zone_mu_and_b_.84 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                                     b = b + b_SA*S_America, 
                                                                     mu = mu + mu_SA*S_America,
                                                             K = 0.84),
                                           start_lower = list(b = -5, mu = min(d_A$abs_lat), b_SA = -5, mu_SA = -10),
                                           start_upper = list(b = 5, mu = max(d_A$abs_lat), b_SA = 5, mu_SA = 10),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = d_A) 
summary(fit_lat_zone_mu_and_b_.84) # b_SA is not significant
broom::tidy(fit_lat_zone_mu_and_b_.84) 
glance(fit_lat_zone_mu_and_b_.84)
broom::tidy(fit_lat_zone_mu_and_b_.84) %>%
  write.table(., "results/fit_lat_zone_mu_and_b_84.txt", 
              quote = F, col.names = T, row.names = F, sep = "\t")

coef_lat_zone_mu_b_.84 = tidy(fit_lat_zone_mu_and_b_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_NA = abs(4/b),
         width_SA = abs(4/(b + b_SA)))
coef_lat_zone_mu_b_.84

# only allow mu to vary, not slope b
# rescale by 0.84:
fit_lat_zone_mu_.84 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                             b = b, 
                                                             mu = mu + mu_SA*S_America,
                                                             K = 0.84),
                                           start_lower = list(b = -5, mu = min(d_A$abs_lat), mu_SA = -10),
                                           start_upper = list(b = 5, mu = max(d_A$abs_lat), mu_SA = 10),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = d_A)
summary(fit_lat_zone_mu_.84)
tidy(fit_lat_zone_mu_.84)
glance(fit_lat_zone_mu_.84)

coef_lat_zone_mu_.84 = tidy(fit_lat_zone_mu_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_NA = abs(4/b))

# don't allow any slope or center differences between clines:
# fit latitude
fit_lat_.84 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                               b = b, 
                                               mu = mu,
                                               K = 0.84),
                             start_lower = list(b = -5, mu = min(d_A$abs_lat)),
                             start_upper = list(b = 5, mu = max(d_A$abs_lat)),
                             supp_errors = 'Y',
                             iter = 250,
                             convergence_count = 100,
                             data = d_A)
summary(fit_lat_.84)
tidy(fit_lat_.84)
glance(fit_lat_.84)
coef_lat_.84 = tidy(fit_lat_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width = abs(4/b))
coef_lat_.84


# ANOVA compare nested models
anova(fit_lat_zone_mu_and_b_.84, fit_lat_zone_mu_.84) # not sig.
#anova(fit_lat_zone_mu_.84, fit_lat_.84) # is sig.
anova_2clines <- anova(fit_lat_zone_mu_and_b_.84, fit_lat_zone_mu_.84, fit_lat_.84)
#anova_2clines
#xtable(anova_2clines)
print(xtable(anova_2clines, 
             caption = "\\color{Gray} \textbf{Latitudinal clines} Significance test for different cline centers (sig.) and cline center and slopes (not sig.) between North and South America",
             label = "anova_lat_clines",
             type = "latex", 
             latex.environments = NULL), 
      file = "../../bee_manuscript/tables/anova_lat_clines.tex")


#----- at what latitude do these clines cross the 50% A ancestry frequency?---------#
find_x <- function(mu, K, A, b){ # solve cline equation to get x
  log(K/A - 1)/-b + mu
}
fifty_perc_lat <- c(CA = find_x(mu = coef_lat_zone_mu_b_.84$mu, 
                                K = 0.84, A = 0.5, 
                               b =coef_lat_zone_mu_b_.84$b),
                    AR = find_x(mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA, 
                                K = 0.84, A = 0.5, 
                                b =coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA))
fifty_perc_lat

# ------------------------ plot wing and ancestry clines -------------------#
# plot the two clines:
# ggplot version
p_A_cline_loess <- filter(d_A, !is.na(alpha)) %>%
  ggplot(data = ., aes(x = abs_lat, y = alpha)) +
  geom_smooth(aes(fill = continent),
              method = "loess",
              lty = 2, lwd = 0,
              fullrange = T) +
  geom_point(aes(color = continent, 
                 shape = population == "Avalon_2014"), size = 2) +
  xlab("Degrees latitude from the equator") +
  ylab("Proportion African ancestry") +
  scale_x_continuous(position = "bottom", breaks = c(30, 32, 34, 36, 38), labels = waiver()) +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  ylim(c(-0.05, 1)) +
  theme_classic() +
  stat_function(fun = function(x) logistic4(x, 
                                            b = coef_lat_zone_mu_b_.84$b, 
                                            mu = coef_lat_zone_mu_b_.84$mu,
                                            K = 0.84),
                color = col_NA_SA_both["N. America"],
                size = 1) +
  stat_function(fun = function(x) logistic4(x, 
                                            b = coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA, 
                                            mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA,
                                            K = 0.84),
                color = col_NA_SA_both["S. America"],
                size = 1) +
  
  geom_vline(data = data.frame(continent = c("N. America", "S. America"),
                               abs_lat = fifty_perc_lat), 
             aes(xintercept = abs_lat),
             color = col_NA_SA_both[c("N. America", "S. America")],
             lty = 3) +
  scale_fill_manual(values = col_NA_SA_both) +
  scale_shape_manual(values = c(1,6)) +
  guides(shape = "none",
         fill = "none", 
         color = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.position="top")
p_A_cline_loess
ggsave("plots/A_cline_prediction_vs_loess.png", height = 3, width = 5.2)


# phenotypic wing clines for comparison:
p_wing_cline_pops <- 
  d_A %>%
  mutate(., wing_length_mm = wing_cm*10) %>%
  filter(!is.na(wing_length_mm)) %>%
  group_by(., population) %>%
  summarise(wing_length_mm = mean(wing_length_mm, na.rm = T),
            abs_lat = mean(abs_lat)) %>%
  left_join(., meta.pop, by = "population") %>%
  rename(continent = zone) %>%
  ggplot(., aes(x = abs_lat, y = wing_length_mm)) +
  ylab("Wing length (mm)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic() +
  geom_point(aes(shape = continent, fill = continent, 
                 color = continent), size = 2) + # shape = 17
  scale_x_continuous(position = "top", breaks = c(30, 32, 34, 36, 38),
                     limits = range(d_A$abs_lat),
                     labels = waiver()) +
  theme(axis.title.x=element_blank()) +
  scale_y_reverse(
    position = "left", 
    breaks = c(9.0, 8.5, 8.0), 
    labels = c("9.00", "8.50", "8.00"), 
    limits = c(9,8)) +
  guides(color = "none", shape = "none", fill = "none") +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_manual(values = col_NA_SA_both)
p_wing_cline_pops

# plot predicted wing cline
p_wing_cline_loess <- d_A %>%
  mutate(., wing_length = wing_cm*10) %>%
  filter(., !is.na(wing_length)) %>%
  ggplot(., aes(x = abs_lat, y = wing_length)) +
  geom_smooth(aes(fill = continent), method = "loess",
              color = "black", lwd = 0, lty = 2) +
  geom_point(aes(color = continent), shape = 1, 
             size = 2) +
  xlab("Degrees latitude from the equator") +
  ylab("Wing length (mm)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic() +
  scale_y_reverse(
    position = "left", 
    breaks = c(9.0, 8.5, 8.0),# set breaks but not limits
    labels = c("9.00", "8.50", "8.00")) +
  scale_x_continuous(position = "top", breaks = c(30, 32, 34, 36, 38),
                     limits = range(d_A$abs_lat),
                     labels = waiver()) +
  stat_function(fun = function(x) -10*(logistic4(x, # need to plot negative because not affected by reverse y scale
                                            b = coef_lat_zone_mu_b_.84$b, 
                                            mu = coef_lat_zone_mu_b_.84$mu,
                                            K = 0.84)*coefficients(m_wing)[2] +
                  coefficients(m_wing)[1]),
                color = col_NA_SA_both["N. America"],
                size = 1,
                lty = 5) +
  stat_function(fun = function(x) -10*(logistic4(x, 
                                            b = coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA, 
                                            mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA,
                                            K = 0.84)*coefficients(m_wing)[2] +
                  coefficients(m_wing)[1]),
                color = col_NA_SA_both["S. America"],
                size = 1,
                lty = 5) +
  scale_fill_manual(values = col_NA_SA_both) +
  guides(fill = "none", color = "none") +
  theme(axis.title.x=element_blank())
p_wing_cline_loess
ggsave("plots/wing_cline_prediction_vs_loess.png", 
       height = 4, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/wing_cline_prediction_vs_loess.png", 
       height = 4, width = 5.2, units = "in", dpi = 600)

# combine plots -- put genetic and phenotypic (wing) clines together:
grid.arrange(p_A_cline_loess, p_wing_cline_loess, nrow = 2, heights = c(5,3))
ggsave("plots/A_and_wing_clines.png", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/A_and_wing_clines.tif", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600,
       device = "tiff", compression = "lzw", type = "cairo")
ggsave("../../bee_manuscript/figures_main/Fig2.tif", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600,
       device = "tiff", compression = "lzw", type = "cairo")
ggsave("../../bee_manuscript/figures/A_and_wing_clines.png", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)


# --------------- latitude vs. distance vs. climate -------------------------------#
# compare model fits for latitude vs. distance from sao paulo vs. climate vars
model_vars = c("abs_lat", "km_from_sao_paulo", "AnnualMeanTemp", "MinTempColdestMonth", "MeanTempColdestQuarter", "AnnualPrecip")
model_names = c("lat", "dist", "temp", "cold", "winter", "precip")
model_info = data.frame(var = model_vars, 
                        name = model_names, 
                        predictor = c("Latitude", 
                                       "Distance to Sao Paulo", 
                                       "Mean temperature", 
                                       "Minimum temperature of coldest month", 
                                       "Mean temperature of coldest quarter",
                                       "Annual precipitation"),
                        stringsAsFactors = F)


# function to fit cline models
fit_cline_m <- function(var, data) nls_multstart(alpha ~ logistic4(x = x, 
                                                  b = b, 
                                                  mu = mu,
                                                  K = 0.84),
                                start_lower = list(b = -5, mu = min(d_A[ , var])),
                                start_upper = list(b = 5, mu = max(d_A[ , var])),
                                supp_errors = 'Y',
                                iter = 250,
                                convergence_count = 100,
                                data = data %>%
                                  mutate(x = data[ , var]))
# run each model:
m_all <- lapply(model_vars, function(x) fit_cline_m(var = x, data = d_A))
t_all <- do.call(rbind, lapply(1:length(m_all), function(i) glance(m_all[[i]]) %>% 
                                 mutate(model = model_names[i]))) %>%
  arrange(AIC) %>%
  left_join(., model_info, by = c("model"="name")) %>%
  mutate(dAIC =  round(AIC - min(AIC), 1),
         rel_lik = exp(- 0.5 * dAIC),
         weight = round(rel_lik / sum(rel_lik), 3)) %>%
  dplyr::select(predictor, df.residual, deviance, dAIC, weight)
t_all
#xtable(t_all)
print(xtable(t_all, 
             caption = "\\color{Gray} \textbf{Cline model comparison} Model rankings between logistic cline fits for African ancestry predicted by climate and distance variables (n = 313 bees).",
             label = "AIC_climate_clines",
             type = "latex", 
             latex.environments = NULL), 
      file = "../../bee_manuscript/tables/AIC_cline_fits_climate.tex")


### --------------------- INDIVIDUAL SNP CLINES --------------------------------- #####
mvn_zero <- list(info = read.table("results/ind_snp_nls_multistart_clines/MVNsim_zero_bounded.info", # data simulation
                             header = T, stringsAsFactors = F),
            params = read.table("results/ind_snp_nls_multistart_clines/MVNsim_zero_bounded.params",
                header = T, stringsAsFactors = F))
mvn <- list(info = read.table("results/ind_snp_nls_multistart_clines/MVNsim_bounded.info", # data simulation
                                   header = T, stringsAsFactors = F),
                 params = read.table("results/ind_snp_nls_multistart_clines/MVNsim_bounded.params",
                                     header = T, stringsAsFactors = F))
sets <- paste0("A", 1:5)
clines <- list(info = do.call(rbind, lapply(sets, function(s)
                                            read.table(paste0("results/ind_snp_nls_multistart_clines/", s, ".info"), # data subset
                                 header = T, stringsAsFactors = F))),
               params = do.call(rbind, lapply(sets, function(s)
                 read.table(paste0("results/ind_snp_nls_multistart_clines/", s, ".params"), # data subset
                            header = T, stringsAsFactors = F))))
# check all models converged -- great!
table(mvn$info$isConv)
table(mvn_zero$info$isConv)
table(clines$info$isConv)

# mean and spread of parameter values across ind. snps
clines$params %>%
  mutate(term = ifelse(term == "b", "w", term),
         estimate = ifelse(term == "w", 4/estimate, estimate)) %>%
  group_by(term) %>%
  summarise(
    mean = mean(estimate),
    min = min(estimate),
    low = quantile(estimate, 0.025),
    median = quantile(estimate, 0.5),
    high = quantile(estimate, 0.975),
    max = max(estimate))

# estimated in km:
clines$params %>%
  mutate(term = ifelse(term == "b", "w_km", term),
         estimate = ifelse(term == "w_km", 4/estimate*km_per_degree_lat, estimate)) %>%
  group_by(term) %>%
  summarise(
    mean = mean(estimate),
    min = min(estimate),
    low = quantile(estimate, 0.025),
    median = quantile(estimate, 0.5),
    high = quantile(estimate, 0.975),
    max = max(estimate))
# 95% range
diff(quantile(clines$params[clines$params$term == "mu", "estimate"], c(0.025, 0.975))) # 95% of clines have centers within this distance

# violin plot for spread of individual SNP clines
d_clines_combined <- rbind(mutate(mvn$params, data = "Simulated SNPs (MVN)"), 
      mutate(clines$params, data = "Observed SNPs")) %>%
  mutate(term = ifelse(term == "b", "Width", "Center"),
         estimate = ifelse(term == "Width", 4/estimate, estimate)) %>%
  left_join(., mutate(A_AR_CA, snp_index = 1:nrow(A_AR_CA), data = "Observed SNPs"), by = c("snp_index", "data")) %>%
  mutate(data = factor(data, levels = c("Simulated SNPs (MVN)", "Observed SNPs"), ordered = T)) %>%
  mutate(outlier_status = ifelse(data == "Simulated SNPs (MVN)", "", 
                          ifelse(!is.na(FDR_AR_high), "(High A)", ifelse(!is.na(FDR_AR_low), "(Low A)", "(Non-outliers)"))))

violin_ind_clines <- d_clines_combined %>%
  bind_rows(., filter(d_clines_combined, data == "Observed SNPs") %>% mutate(outlier_status = "(All)")) %>%
  mutate(group = paste(data, outlier_status)) %>%
  mutate(group = factor(group, levels = c("Simulated SNPs (MVN) ", "Observed SNPs (All)", "Observed SNPs (Non-outliers)", "Observed SNPs (High A)", "Observed SNPs (Low A)"), ordered = T)) %>%
  ggplot(.) +
  geom_violin(aes(y = estimate, x = group, fill = group, color = group)) +
  facet_wrap(~term, scales = "free_y") +
  ylab("Degrees latitude") +
  xlab("") +
  theme_classic() +
  scale_fill_manual(values = c(viridis(4)[1:2], magma(3, end = .9)), name = "") +
  scale_color_manual(values = c(viridis(4)[1:2], magma(3, end = .9)), name = "") +
  guides(color = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

violin_ind_clines
ggsave("plots/violin_ind_snp_clines.png",
       plot = violin_ind_clines,
       height = 5.2, width = 5.2, units = "in")
# save in figures for manuscript:
ggsave("../../bee_manuscript/figures/violin_ind_snp_clines.png", 
       plot = violin_ind_clines,
       height = 5.2, width = 5.2, units = "in", device = "png", dpi = 600)
ggsave("../../bee_manuscript/figures_supp/violin_ind_snp_clines.tiff", 
       plot = violin_ind_clines,
       height = 5.2, width = 5.2, units = "in", device = "tiff", dpi = 600)

# --------------- Any individual SNP outliers for cline steepness? ------------------------#
# what is the FDR? Uses script ../local_ancestry/calc_FDRs.R
# make seq of values to calculate FDR at
b_range = seq(min(clines$params$estimate[clines$params$term == "b"]), 
              max(clines$params$estimate[clines$params$term == "b"]), 
              length.out = 1000)
mu_range = seq(min(clines$params$estimate[clines$params$term == "mu"]), 
             max(clines$params$estimate[clines$params$term == "mu"]), 
             length.out = 1000)

# false discovery rate calculations
fdr_high_b <- sapply(b_range, function(b) fdr_1pop_high(a = b, pop = clines$params$estimate[clines$params$term == "b"], 
                                                        sims = mvn_zero$params$estimate[mvn_zero$params$term == "b"]))
fdr_low_b <- sapply(b_range, function(b) fdr_1pop_low(a = b, pop = clines$params$estimate[clines$params$term == "b"], 
                                                        sims = mvn_zero$params$estimate[mvn_zero$params$term == "b"]))
fdr_high_mu <- sapply(mu_range, function(m) fdr_1pop_high(a = m, pop = clines$params$estimate[clines$params$term == "mu"], 
                                                        sims = mvn_zero$params$estimate[mvn_zero$params$term == "mu"]))
fdr_low_mu <- sapply(mu_range, function(m) fdr_1pop_low(a = m, pop = clines$params$estimate[clines$params$term == "mu"], 
                                                      sims = mvn_zero$params$estimate[mvn_zero$params$term == "mu"]))
# combine to set FDR threshold
FDR_clines = data.frame(FDR_values = FDR_values,
                  high_b = sapply(FDR_values, function(p) ifelse(min(fdr_high_b) <= p, # are there any sig snps?
                                                                     max(b_range[fdr_high_b > p], na.rm = T),
                                                                 NA)),
                  low_b = sapply(FDR_values, function(p) ifelse(min(fdr_low_b) <= p, # are there any sig snps?
                                                                min(b_range[fdr_low_b > p], na.rm = T),
                                                                NA)),
                  high_mu = sapply(FDR_values, function(p) ifelse(min(fdr_high_mu) <= p, # are there any sig snps?
                                                                  max(mu_range[fdr_high_mu > p], na.rm = T),
                                                                  NA)),
                  low_mu = sapply(FDR_values, function(p) ifelse(min(fdr_low_mu) <= p, # are there any sig snps?
                                                                 min(mu_range[fdr_low_mu > p], na.rm = T),
                                                                 NA)),
                  stringsAsFactors = F)
FDR_clines # note: NA false-discovery-rate value for steep clines (high_b) means there are no sig. outliers meeting that FDR

write.table(FDR_clines, "results/FDRs_ind_snp_nls_multistart_clines.txt", quote = F, col.names = T, row.names = F, sep = "\t")


##### -------------------- ENRICHMENT OF STEEP CLINES IN REGIONS OF THE GENOME WITH LOW RECOMB? ---------------######


# calculation of enrichment for steep clines in regions of the genome with low recombination rates
load("../local_ancestry/results/sites_r.RData") # load recombination rate data by sites
# this is run 2 ways with qualitatively the same results (top 1% vs. top 5% steepest clines), so top 5% is reported in the manuscript
steepest_clines <- clines$params$estimate[clines$params$term == "b"] > 
                           quantile(clines$params$estimate[clines$params$term == "b"], .99)
steepest_clines05 <- clines$params$estimate[clines$params$term == "b"] > # top 5% 
                           quantile(clines$params$estimate[clines$params$term == "b"], .95)
table(sites_r$r_bin5[steepest_clines])
table(sites_r$r_bin5)
# % of snps in each recombination rate quintile that are part of the steepest clines
table(sites_r$r_bin5[steepest_clines])/table(sites_r$r_bin5)
table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5)
# (will bootstrap conf. intervals)
# how are the outliers distributed across chromosomes? 
# Top 5% on all chr. Top 1% on all chr except 16.
sites_r %>%
  mutate(steepest_clines = steepest_clines,
         steepest_clines05 = steepest_clines05) %>%
  group_by(chr) %>%
  summarise(n = n(),
            steepest_clines05 = sum(steepest_clines05),
            steepest_clines01 = sum(steepest_clines))

# plot across genome:
sites_r %>%
  mutate(., b = clines$params$estimate[clines$params$term == "b"],
         percentile = sapply(1:nrow(sites_r), function(i) 
           ifelse(steepest_clines[i], "1%", 
                  ifelse(steepest_clines05[i], "5%", "n.s.")))) %>%
  filter(., percentile != "n.s.") %>% # don't plot non-sig points
  ggplot(., aes(x = pos, y = b, color = percentile)) +
  geom_point(size = 0.1) +
  facet_wrap(~chr)

# for bootstrap we divide the genome into 0.2cM windows (see script divide_map_0.2cM_windows.R):
rmap2 <- read.table("../data/recomb_map/Wallberg_HAv3.1/map_10kb_in_0.2cM_windows.bed", stringsAsFactors = F , 
                    header = T , sep = "\t")
# associate each site with an ancestry call with a 0.2cM window
sites_r_windows0 <- bedr(
  engine = "bedtools", 
  input = list(a = dplyr::select(sites_r, chr, pos) %>% # take just center position & find window
                 mutate(start = pos - 1) %>%
                 mutate(end = pos) %>%
                 dplyr::select(., chr, start, end),
               b = rmap2), 
  method = "map", 
  params = "-g ../data/honeybee_genome/chr.lengths.by.name -c 4 -o mean",
  check.chr = F
) %>%
  data.table::setnames(c("chr", "ignore", "pos", "window_0.2cM")) %>%
  dplyr::select(., chr, pos, window_0.2cM) %>%
  mutate(pos = as.integer(pos), window_0.2cM = as.integer(window_0.2cM)) %>%
  left_join(sites_r, ., 
            by = c("chr", "pos"))
sites_r_windows = sites_r_windows0 %>%
  group_by(window_0.2cM) %>%
  summarise(n = n()) %>%
  left_join(sites_r_windows0, ., by = "window_0.2cM") %>%
  mutate(b = clines$params$estimate[clines$params$term == "b"]) # add slope

# citation for block bootstrap with genomic data: Bickel, P. J., Boley, N., Brown, J. B., Huang, H., & Zhang, N. R. (2010). Subsampling Methods for Genomic Inference.The Annals ofApplied Statistics,4(4), 1660-1697.http://dx.doi.org/10.1214/10-AOAS363
# this follows ENCODE Strategy..generally though they deal with homogeneous blocks..and it's overly conservative if blocks are not homogeneous...better to segment the blocks in this case.
# other possible publication to look at : https://www.jstor.org/stable/2290993?seq=1#metadata_info_tab_contents

top_1perc = quantile(sites_r_windows$b, 0.99)
top_5perc = quantile(sites_r_windows$b, 0.95)

save(file = "results/sites_r_windows_steep_slopes.RData",
     list = c("sites_r_windows", "top_1perc", "top_5perc"))


# bootstrapping moved to it's own script -- bootstrap_steep_clines.R                                  
# read in and summarise bootstraps
boot <- do.call(rbind, lapply(1:100, function(b) 
  read.table(paste0("results/BOOT_r_steep_clines/boot_", b, ".txt"), 
             header = T, stringsAsFactors = F))) %>%
  mutate(., top5_orig_lowr_perc = top5_orig_1/n_1,
         top5_sample_lowr_perc = top5_sample_1/n_1,
         top5_orig_highr_perc = top5_orig_5/n_5,
         top5_sample_highr_perc = top5_sample_5/n_5,
         top1_sample_lowr_perc = top1_sample_1/n_1,
         top1_sample_highr_perc = top1_sample_5/n_5,
         top1_orig_lowr_perc = top1_orig_1/n_1,
         top1_orig_highr_perc = top1_orig_5/n_5)

# 'basic' bootstrap otherwise known as 'pivot conf. intervals' = 2*estimate - 95% bootstrap endpoints:
# note: does not subtract mean bootstrap, but full sample mean, because we expect bias to be the same from pop -> sample as from sample -> bootstrap, i.e.
# 'The population is to the sample as the sample is to the bootstrap samples.'
# estimate lowest recomb. bin bootstrap for % of SNPs that are in the top 5% steepest SNPs genomewide:
lowr_estimate = (table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5))[[1]]

# basic bootstrap conf. intervals
2*lowr_estimate - quantile(boot$top5_sample_lowr_perc, c(.975, .025))
lowr_estimate - quantile(boot$top5_sample_lowr_perc - lowr_estimate, c(0.975, 0.025)) # same

# plot bootstrap lowr
hist(boot$top5_sample_lowr_perc, main = "low r bootstrap - % steep clines")
abline(v = lowr_estimate, col = "orange")
abline(v = 2*lowr_estimate - quantile(boot$top5_sample_lowr_perc, c(.975, .025)), col = "blue")
abline(v = 0.05, col = "red")

# estimate highest recomb. bin bootstrap for % of SNPs that are in the top 5% steepest SNPs genomewide:
highr_estimate = (table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5))[[5]]
# basic bootstrap intervals for highest recomb. bin:
2*highr_estimate - quantile(boot$top5_sample_highr_perc, c(.975, .025))
# plot bootstrap highr
hist(boot$top5_sample_highr_perc, main = "high r bootstrap - % steep clines", xlim = range(c(0.05, boot$top5_sample_highr_perc)))
abline(v = highr_estimate, col = "orange")
abline(v = 2*highr_estimate - quantile(boot$top5_sample_highr_perc, c(.975, .025)), col = "blue")
abline(v = 0.05, col = "red")

with(boot, hist(top5_cutoff)) # 5% cutoff is very similar across bootstraps, 
# so makes little diff. if I use original fixed cutoff or bootstrap cutoff.


# bootstrap difference in mean slope between high and low r bins:
mean_b_r_bin5 = sites_r %>%
  mutate(., b = clines$params$estimate[clines$params$term == "b"]) %>%
  group_by(r_bin5) %>%
  summarise(mean_b = mean(b),
            mean_w = mean(4/b))
diff_mean_slope_estimate = as.numeric(mean_b_r_bin5[5, "mean_b"] - mean_b_r_bin5[1, "mean_b"])
diff_mean_slope_estimate
# basic bootstrap conf. intervals
2*diff_mean_slope_estimate - quantile(boot$mean_b_5 - boot$mean_b_1, c(.975, .025))
4/(2*diff_mean_slope_estimate - quantile(boot$mean_b_5 - boot$mean_b_1, c(.975, .025))) # wrong
(2*4/diff_mean_slope_estimate - quantile(4/boot$mean_b_5 - 4/boot$mean_b_1, c(.975, .025))) # also wrong
diff_mean_w_estimate = mean_b_r_bin5[5, "mean_w"] - mean_b_r_bin5[1, "mean_w"]
diff_mean_w_estimate*111
# I'd need to redo the bootstrap if I want to get CI for means of w instead of means of b for each boot

# plot bootstrap lowr
hist(boot$mean_b_5 - boot$mean_b_1, main = "diff mean cline slope high to low r", xlim = range(c(0, boot$mean_b_5 - boot$mean_b_1)))
abline(v = diff_mean_slope_estimate, col = "orange")
abline(v = 2*diff_mean_slope_estimate - quantile(boot$mean_b_5 - boot$mean_b_1, c(.975, .025)), col = "blue")
abline(v = 0, col = "red")
# mean slope has a small but sig. difference across recombination rate quintiles
mutate(sites_r, b = clines$params$estimate[clines$params$term == "b"]) %>%
  group_by(r_bin5) %>%
  summarise(mean_b = mean(b),
            mean_w = mean(4*km_per_degree_lat/b))
mutate(sites_r, b = clines$params$estimate[clines$params$term == "b"]) %>%
  mutate(w = 4*km_per_degree_lat/b) %>%
  lm(data = ., w ~ r_bin5) %>%
  summary()

# plot steep outliers across the genome, colored by recombination rate
cbind(sites_r, clines$params[clines$params$term == "b", ]) %>%
  filter(estimate >= quantile(estimate, .99)) %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = estimate, color = r_bin5_factor)) +
  facet_wrap(~chr) +
  labs(color = "recombination rate") +
  ggtitle("Top 1% steepest SNP clines by r")