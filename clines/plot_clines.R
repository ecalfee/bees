# script to plot clines in ancestry, e.g. with latitude and environment
library(sp)
library(gridExtra)
library(raster)
library(geosphere)
library(rethinking)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(purrr)
library(nlstools)
library(betareg) # alternative ML fitting
library(xtable)
source("../colors.R") # get color palette
#install.packages("nls.multstart")
library(nls.multstart)
source("../local_ancestry/calc_FDRs.R")
old.par <- par() # save default
source("cline_functions.R") # loads logistic and stepped clines etc.
load("results/d_A.RData") # load data
load("../wing_analysis/results/wing_fits.RData") # wing linear model fits
load("../local_ancestry/results/meta.RData")

iHot <- which.max(d_A$AnnualMeanTemp)
iCold <- which.min(d_A$AnnualMeanTemp)

distm(d_A[iHot, c("long", "lat")], d_A[iCold, c("long", "lat")], 
      fun = distGeo)/1000 # distance in km between hottest and coldest site (CA mountains and desert)

# plots
d_A %>%
  ggplot(., aes(y = MeanTempColdestQuarter, x = AnnualMeanTemp, color = continent, shape = year)) +
  geom_point()
d_A %>%
  ggplot(., aes(y = MeanTempColdestQuarter, x = abs(lat), color = continent, shape = year)) +
  geom_point()
d_A %>%
  ggplot(., aes(y = AnnualMeanTemp, x = abs(lat), color = continent, shape = year)) +
  geom_point()
d_A %>%
  ggplot(., aes(x = abs(lat), y = AnnualPrecip, color = continent, shape = year)) +
  geom_point() # basically S. America is wet and N. America is dry

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
# we have a cold/wet mountain top at 33.7 degrees lat
# and a dry/hot desert at 33.8 degrees lat



# A ancestry vs. latitude:
d_A %>%
  ggplot(., aes(x = abs(lat), y = alpha, color = continent)) +
  geom_point() +
  scale_color_manual(values = col_NA_SA_both)

############-----logistic clines w/ nls()------#########


# what width would be expect for neutral diffusion?
km_per_degree_lat = 111 # 111.699 km at the poles, 110.567 estimate from the equator
curve((coef_lat$b*4*km_per_degree_lat)^2/(2*pi*x), from = 20, to = 90,
      xlab = "generations since admixture",
      ylab = "neutral s.d. parent-offspring dispersal")
years_tot = 2018-1957
dispersal_kernel <- (coef_lat$b*4*km_per_degree_lat)^2/(2*pi*years_tot)
dispersal_kernel # variance in distance between parent and offspring
660/(1.68*sqrt(years_tot))
# a dispersal kernel variance of ~130km/year seems like a lot,  
# but may still be low for the rate of spread we saw in the Africanized honey bee invasion:
# in a gaussian dispersal kernel, this would only gain bees ~5 km a year
x = matrix(0, years_tot+1, 1000)
for (i in 2:nrow(x)){
  x[i, ] <- sapply(x[i - 1, ], function(m) rnorm(n = 1, mean = m, sd = sqrt(dispersal_kernel)))
}
max_pos = apply(x, 1, max)
year = 1:nrow(x)
plot(year, max_pos)
abline(lm(max_pos ~ year))
max(max_pos)/nrow(x) # mean distance gained per year

#1957 to 2018 is approx 30 to 60 generations depending on whether you believe 1 gen every year or every other year
(2018-1970)/2 # ~24 generations if you think of time since reaching current cline center


# ok what about fitting 1 model with a SA variable? joint fit. No rescaling
fit_lat_zone_mu_and_b <- nls_multstart(alpha ~ logistic3(x = abs_lat, 
                                                                 b = b + b_SA*S_America, 
                                                                 mu = mu + mu_SA*S_America),
                                       start_lower = list(b = -5, mu = min(d_A$abs_lat), b_SA = -5, mu_SA = -10),
                                       start_upper = list(b = 5, mu = max(d_A$abs_lat), b_SA = 5, mu_SA = 10),
                                       supp_errors = 'Y',
                                       iter = 250,
                                       convergence_count = 100,
                                       data = d_A) #%>%
                                         #filter(year != 2014))
summary(fit_lat_zone_mu_and_b) # b_SA is not significant, drop:
coef_lat_zone_mu_b = tidy(fit_lat_zone_mu_and_b) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_NA = abs(4/b),
         width_SA = abs(4/(b + b_SA)))
coef_lat_zone_mu_b

# rescale by 0.84:
fit_lat_zone_mu_and_b_.84 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                                     b = b + b_SA*S_America, 
                                                                     mu = mu + mu_SA*S_America,
                                                             K = 0.84),
                                           start_lower = list(b = -5, mu = min(d_A$abs_lat), b_SA = -5, mu_SA = -10),
                                           start_upper = list(b = 5, mu = max(d_A$abs_lat), b_SA = 5, mu_SA = 10),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = d_A) #%>%
#filter(year != 2014))
summary(fit_lat_zone_mu_and_b_.84) # b_SA is not significant, drop:
tidy(fit_lat_zone_mu_and_b_.84)
glance(fit_lat_zone_mu_and_b)
glance(fit_lat_zone_mu_and_b_.84)

coef_lat_zone_mu_b_.84 = tidy(fit_lat_zone_mu_and_b_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_NA = abs(4/b),
         width_SA = abs(4/(b + b_SA)))
coef_lat_zone_mu_b_.84
coef_lat_zone_mu_b_.84$mu
sum(coef_lat_zone_mu_b_.84[ , c("mu", "mu_SA")])
coef_lat_zone_mu_b

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
                                           data = d_A) #%>%
#filter(year != 2014))
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
anova(fit_lat_zone_mu_.84, fit_lat_.84) # is sig.
anova_2clines <- anova(fit_lat_zone_mu_and_b_.84, fit_lat_zone_mu_.84, fit_lat_.84)
anova_2clines
xtable(anova_2clines)
print(xtable(anova_2clines, 
             caption = "\\color{Gray} \textbf{Latitudinal clines} Significance test for different cline centers (sig.) and cline center and slopes (not sig.) between North and South America",
             label = "anova_lat_clines",
             type = "latex", 
             latex.environments = NULL), 
      file = "../../bee_manuscript/tables/anova_lat_clines.tex")



# bootstraps for uncertainty in cline parameter estimates:
# (individuals re-sampled within populations)
# bootstraps:
boots_logistic <- read.table("results/bootstrap_logistic_cline_seed100.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
table(boots_logistic$isConv) # all 1000 converged -- NOTE re-doing for 10000 bootstraps
# pivot=basic bootstrap CI's= estimate - quantile(boots - estimate, c(0.975, 0.025))
# NA center
CI_NA_mu = coef_lat_zone_mu_b_.84$mu - quantile(boots_logistic$mu - coef_lat_zone_mu_b_.84$mu, c(0.975, 0.025))
# SA center
CI_SA_mu = (coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA) - quantile(boots_logistic$mu + boots_logistic$mu_SA - (coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA), c(0.975, 0.025))
# diff center
CI_diff_mu = coef_lat_zone_mu_b_.84$mu_SA - quantile(boots_logistic$mu_SA - coef_lat_zone_mu_b_.84$mu_SA, c(0.975, 0.025))
# NA slope
CI_NA_b = coef_lat_zone_mu_b_.84$b - quantile(boots_logistic$b - coef_lat_zone_mu_b_.84$b, c(0.975, 0.025))
# SA slope
CI_SA_b = (coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA) - quantile(boots_logistic$b + boots_logistic$b_SA - (coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA), c(0.975, 0.025))
# diff slope
CI_diff_b = coef_lat_zone_mu_b_.84$b_SA - quantile(boots_logistic$b_SA - coef_lat_zone_mu_b_.84$b_SA, c(0.975, 0.025))

coef_lat_zone_mu_b_.84
CI_SA_mu
CI_NA_mu
CI_diff_mu
CI_SA_b
CI_NA_b
CI_diff_b

4/CI_NA_b
4/CI_SA_b
4/CI_diff_b # not right



CI_NA_w = 4/coef_lat_zone_mu_b_.84$b - quantile(4/boots_logistic$b - 4/coef_lat_zone_mu_b_.84$b, c(0.975, 0.025))
CI_NA_w
# SA slope
CI_SA_w = (4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA)) - quantile(4/(boots_logistic$b + boots_logistic$b_SA) - 4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA), c(0.975, 0.025))
CI_SA_w
CI_diff_w2 = (4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA) - 4/(coef_lat_zone_mu_b_.84$b)) - 
                quantile(4/(boots_logistic$b + boots_logistic$b_SA) - 4/boots_logistic$b - 
                           (4/(coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA) - 4/coef_lat_zone_mu_b_.84$b), c(0.975, 0.025))
CI_diff_w2
curve(-4/(-0.5 + x), from = -3, to = 0)
CI_SA_w-CI_NA_w


# change parameterization to estimate w directly with bootstrap:
# no rescale:
fit_lat_zone_mu_and_w_1 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                             b = -4/(w + w_SA*S_America), 
                                                             mu = mu + mu_SA*S_America,
                                                             K = 1),
                                           start_lower = list(w = 0, mu = min(d_A$abs_lat), w_SA = -10, mu_SA = -5),
                                           start_upper = list(w = 15, mu = max(d_A$abs_lat), w_SA = 10, mu_SA = 5),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = d_A)
summary(fit_lat_zone_mu_and_w_1)
tidy(fit_lat_zone_mu_and_w_1)
glance(fit_lat_zone_mu_and_w_1)

# need to flip the sign of b for correct direction of cline:
curve(logistic4(x, b = .5, mu = 32, K = 0.84), from = 25, to = 40)
curve(logistic4(x, b = -.5, mu = 32, K = 0.84), from = 25, to = 40, add = T, lty = 2)
curve(logistic4(-x, b = .5, mu = -32, K = 0.84), from = 25, to = 40, add = T, lty = 3)

# rescale by 0.84:
fit_lat_zone_mu_and_w_.84 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                             b = -4/(w + w_SA*S_America), # if I try to fit this cline for + w, it finds a suboptimal solution 
                                                             mu = mu + mu_SA*S_America,
                                                             K = 0.84),
                                           start_lower = list(w = -15, mu = min(d_A$abs_lat), w_SA = -10, mu_SA = -5),
                                           start_upper = list(w = 15, mu = max(d_A$abs_lat), w_SA = 10, mu_SA = 5),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = d_A)
fit_lat_zone_mu_and_w_.84_SA_only <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                             b = 4/w, 
                                                             mu = mu,
                                                             K = 0.84),
                                           start_lower = list(w = -15, mu = min(d_A$abs_lat)),
                                           start_upper = list(w = 0, mu = max(d_A$abs_lat)),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = filter(d_A, continent == "S. America"))
summary(fit_lat_zone_mu_and_w_.84_SA_only)
summary(fit_lat_zone_mu_and_w_.84)
tidy(fit_lat_zone_mu_and_w_.84)
glance(fit_lat_zone_mu_and_w_.84)

coef_lat_zone_mu_w_.84 = tidy(fit_lat_zone_mu_and_w_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)
coef_lat_zone_mu_w_.84
4/coef_lat_zone_mu_b_.84$b
coef_lat_zone_mu_w_.84$mu_SA
coef_lat_zone_mu_w_.84$mu
sum(coef_lat_zone_mu_w_.84[ , c("mu", "mu_SA")])
coef_lat_zone_mu_w_.84$w
sum(coef_lat_zone_mu_w_.84[ , c("w", "w_SA")])
4/sum(coef_lat_zone_mu_b_.84[ , c("b", "b_SA")]) # same. 

curve(logistic4(x, b = 4/(coef_lat_zone_mu_w_.84$w), mu = coef_lat_zone_mu_w_.84$mu, K = 0.84),
      from = 25, to = 40, col = col_NA_SA_both[["S. America"]])
curve(logistic4(x, b = coef_lat_zone_mu_b_.84$b, mu = coef_lat_zone_mu_b_.84$mu, K = 0.84),
      from = 25, to = 40, add = T, col = "blue", lty = 2)
curve(logistic4(x, b = 4/(coef_lat_zone_mu_w_.84$w + coef_lat_zone_mu_w_.84$w_SA), mu = coef_lat_zone_mu_w_.84$mu + coef_lat_zone_mu_w_.84$mu_SA, K = 0.84),
      from = 25, to = 40, add = T, col = col_NA_SA_both[["N. America"]])

# at what latitude do these clines cross the 50% A ancestry frequency?
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
fifty_perc_lat_mu <- c(CA = find_x(mu = coef_lat_zone_mu_.84$mu, 
                                K = 0.84, A = 0.5, 
                                b = coef_lat_zone_mu_.84$b),
                    AR = find_x(mu = coef_lat_zone_mu_.84$mu + coef_lat_zone_mu_.84$mu_SA, 
                                K = 0.84, A = 0.5, 
                                b =coef_lat_zone_mu_.84$b))
fifty_perc_lat_mu # still <.6 degrees lat diff.
twenty_perc_lat <- c(CA = find_x(mu = coef_lat_zone_mu_b_.84$mu, 
                                 K = 0.84, A = 0.2, 
                                 b =coef_lat_zone_mu_b_.84$b),
                     AR = find_x(mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA, 
                                 K = 0.84, A = 0.2, 
                                 b =coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA))
twenty_perc_lat
fortytwo_perc_lat <- c(CA = find_x(mu = coef_lat_zone_mu_b_.84$mu, 
                                K = 0.84, A = 0.42, 
                                b =coef_lat_zone_mu_b_.84$b),
                    AR = find_x(mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA, 
                                K = 0.84, A = 0.42, 
                                b =coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA))
fortytwo_perc_lat



png("results/sigma_from_est_cline_width.png")
# approximate cline width = slope at steepest part near center:
width_SA = diff(sapply(c(0.41, 0.43), function(p) find_x(mu = coef_lat_zone_mu_b_.84$mu, 
       K = 0.84, A = p, 
       b = coef_lat_zone_mu_b_.84$b)))/-.02 * 111 # 111km/degree
# estimate sigma from a neutral diffusion model of cline width:
curve(width_SA/sqrt(2*pi*x), from = 20, to = 120, n = 100,
      ylim = c(0, 80),
      xlab = "generations t", ylab = "sigma (km/gen)", 
      main = "cline width = sqrt(2*pi*sigma^2*t)")
curve(500/sqrt(2*pi*x), from = 20, to = 120, n = 100, col = "blue", add = T)
legend("topright", legend = c(paste0(round(width_SA), "km"), "500km"), lty = 1, col = c("black", "blue"))
dev.off()


# fit wing length phenotypic cline:
fit_wing_lat_zone_mu_and_b <- nls_multstart(wing_rescaled ~ logistic3(x = abs_lat, 
                                                                 b = b + b_SA*S_America, 
                                                                 mu = mu + mu_SA*S_America),
                                            start_lower = list(b = -5, mu = min(d_A$abs_lat), b_SA = -5, mu_SA = -10),
                                            start_upper = list(b = 5, mu = max(d_A$abs_lat), b_SA = 5, mu_SA = 10),
                                       supp_errors = 'N', # 'Y'
                                       iter = 250,
                                       convergence_count = 100,
                                       data = d_A %>%
                                         mutate(wing_rescaled = 
                                                  (wing_cm - min(d_A$wing_cm, na.rm = T))/diff(range(d_A$wing_cm, na.rm = T))))
summary(fit_wing_lat_zone_mu_and_b)
coef_wing_lat_zone_mu_b = tidy(fit_wing_lat_zone_mu_and_b) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate)  %>%
  mutate(width_NA = abs(4/b),
         width_SA = abs(4/(b + b_SA)))
coef_wing_lat_zone_mu_b # super wide
#d_A %>%
#  mutate(wing_rescaled = 
#           (wing_cm - sum(coefficients(m_wing)*c(1, 0.5)))/abs(coefficients(m_wing)[2])) %>%
d_A %>%
  #mutate(., wing_rescaled = (wing_cm - min(d_A$wing_cm, na.rm = T))/diff(range(d_A$wing_cm, na.rm = T))) %>%
  ggplot(., aes(x = abs(lat), y = wing_cm, color = continent)) +
  geom_point() +
  stat_function(fun = function(x) diff(range(d_A$wing_cm, na.rm = T))*logistic3(x, # need to plot negative because not affected by reverse y scale
                                                 b = coef_wing_lat_zone_mu_b$b, 
                                                 mu = coef_wing_lat_zone_mu_b$mu) + 
                  min(d_A$wing_cm, na.rm = T),
                color = col_NA_SA_both["N. America"],
                size = 1,
                lty = 5) +
  scale_color_manual(values = col_NA_SA_both) +
  geom_smooth() +
  stat_function(fun = function(x) diff(range(d_A$wing_cm, na.rm = T))*logistic3(x, # need to plot negative because not affected by reverse y scale
                                                                                b = coef_wing_lat_zone_mu_b$b + coef_wing_lat_zone_mu_b$b_SA, 
                                                                                mu = coef_wing_lat_zone_mu_b$mu + coef_wing_lat_zone_mu_b$mu_SA) + 
                  min(d_A$wing_cm, na.rm = T),
                color = col_NA_SA_both["S. America"],
                size = 1,
                lty = 5) +
  ggtitle("attempt to fit separate clines to (rescaled) wing lengths")



# plot the two clines:
#plot
logistic3(10, # need to plot negative because not affected by reverse y scale
          b = 4/(fit_wing_lat_zone_mu_and_w[["w"]]), 
          mu = coef(fit_wing_lat_zone_mu_and_w)[["mu"]])

# plot wing cline and genomic cline together:
png("../wing_analysis/plots/fit_wing_A_clines.png", height = 6, width = 8, units = "in", res = 300)
par(mfrow = c(2,1))
par(mar=c(5.1, 4.1, 4.1, 8.1))
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     xlab = "Degrees latitude from equator",
     ylab = "African ancestry proportion")
points(alpha ~ abs_lat, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x, b = coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA, 
                         mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA,
                K = 0.84), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 1, add = T)
#abline(h=0.84, lty = 2) # clutters plot
curve(logistic4(x, b = coef_lat_zone_mu_b_.84$b, 
                mu = coef_lat_zone_mu_b_.84$mu,
                K = 0.84), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 1, add = T)
# lines for 50% theshold
abline(v = fifty_perc_lat, lty = 2,
       col = col_NA_SA_both[c("N. America", "S. America")])
# lines for c, cline center
#abline(v = c(coef_lat_zone_mu_b_.84$mu, coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA), 
#       lty = 2,
#       col = col_NA_SA_both[c("N. America", "S. America")])
# legend
#legend(39, 0.6, c("N. America", "S. America"), xpd = TRUE,
#       pch = c(1, 1), col = col_NA_SA_both[c("N. America", "S. America")])
legend("topright", c("N. America", "S. America"),
       pch = c(1, 1), col = col_NA_SA_both[c("N. America", "S. America")])

# 2nd plot
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     xlab = "Degrees latitude from equator",
     ylab = "A-like Wing Length")
points(1-wing_rescaled ~ abs_lat, data = d_A %>%
         mutate(wing_rescaled = 
                  (wing_cm - min(d_A$wing_cm, na.rm = T))/diff(range(d_A$wing_cm, na.rm = T))), # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]),
       pch = 4)
curve(1-logistic3(x, b = coef_wing_lat_zone_mu_b$b, mu = coef_wing_lat_zone_mu_b$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 2, add = T)
curve(1-logistic3(x, b = coef_wing_lat_zone_mu_b$b + coef_wing_lat_zone_mu_b$b_SA, 
                  mu = coef_wing_lat_zone_mu_b$mu + coef_wing_lat_zone_mu_b$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 2, add = T)


dev.off()
par(mfrow = c(1,1))
par(old.par)

# ggplot version
p_A_cline_loess <-
  ggplot(data = d_A, aes(x = abs_lat, y = alpha)) +
  geom_smooth(data = filter(d_A, !is.na(alpha)), aes(fill = continent), 
              method = "loess",
              lty = 2, lwd = 0,
              fullrange = T) +
  geom_point(aes(color = continent), shape = 1, size = 2) +
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
  guides(fill = "none", 
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
  geom_point(aes(shape = continent, fill = continent, color = continent), size = 2) + # shape = 17
  scale_x_continuous(position = "top", breaks = c(30, 32, 34, 36, 38),
                     limits = range(d_A$abs_lat),
                     labels = waiver()) +
  theme(axis.title.x=element_blank()) +
  #scale_y_continuous(
  #  position = "left", 
  #  breaks = c(8.0, 8.5, 9.0), 
  #  labels = c("8.00", "8.50", "9.00"), 
  #  limits = c(8,9)) +
  scale_y_reverse(
    position = "left", 
    breaks = c(9.0, 8.5, 8.0), 
    labels = c("9.00", "8.50", "8.00"), 
    limits = c(9,8)) +
  
  guides(color = "none", shape = "none", fill = "none") +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_manual(values = col_NA_SA_both) #+
#  xlim(range(d_A$abs_lat))
p_wing_cline_pops
sqrt((500^2)/2/pi/50)

# plot predicted wing cline
p_wing_cline_loess <- 
  d_A %>%
  mutate(., wing_length = wing_cm*10) %>%
  ggplot(., aes(x = abs_lat, y = wing_length)) +

  geom_smooth(aes(fill = continent), method = "loess", 
              color = "black", lwd = 0, lty = 2) +
  geom_point(aes(color = continent), shape = 1, size = 2) +
  xlab("Degrees latitude from the equator") +
  ylab("Wing length (mm)") +
  scale_color_manual(values = col_NA_SA_both, name = NULL) +
  theme_classic() +
  scale_y_reverse(
    position = "left", 
    breaks = c(9.0, 8.5, 8.0), 
    labels = c("9.00", "8.50", "8.00"), 
    limits = c(9, 8)) +
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
#grid.arrange(p_A_cline_loess, p_wing_cline_pops, nrow = 2, heights = c(5,1.5))
grid.arrange(p_A_cline_loess, p_wing_cline_pops, nrow = 2, heights = c(5,3))
ggsave("plots/A_and_wing_pop_clines.png", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_pops, 
                           nrow = 2, heights = c(2,1)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/A_and_wing_pop_clines.tiff", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_pops, 
                           nrow = 2, heights = c(2,1)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/A_and_wing_pop_clines.png", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_pops, 
                           nrow = 2, heights = c(2,1)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)
# alt. version:
grid.arrange(p_A_cline_loess, p_wing_cline_loess, nrow = 2, heights = c(5,3))
ggsave("plots/A_and_wing_clines.png", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures_main/A_and_wing_clines.tiff", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)
ggsave("../../bee_manuscript/figures/A_and_wing_clines.png", 
       plot = grid.arrange(p_A_cline_loess, p_wing_cline_loess, 
                           nrow = 2, heights = c(3,2.5)),
       height = 7.2, width = 5.2, units = "in", dpi = 600)





# overlapping plots:
png("../wing_analysis/plots/fit_wing_A_clines_1_plot_together.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
points(alpha ~ abs_lat, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic3(x, b = coef_lat_zone_mu_b$b + coef_lat_zone_mu_b$b_SA, 
                mu = coef_lat_zone_mu_b$mu + coef_lat_zone_mu_b$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 1, add = T)
curve(logistic3(x, b = coef_lat_zone_mu_b$b, mu = coef_lat_zone_mu_b$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 1, add = T)
curve(0.84*logistic3(x, b = coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA, 
                     mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 3, add = T)
curve(0.84*logistic3(x, b = coef_lat_zone_mu_b_.84$b, mu = coef_lat_zone_mu_b_.84$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 3, add = T)
points(1-wing_rescaled ~ abs_lat, data = d_A %>%
         mutate(wing_rescaled = 
                  (wing_cm - min(d_A$wing_cm, na.rm = T))/diff(range(d_A$wing_cm, na.rm = T))), # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]),
       pch = 4)
curve(1-logistic3(x, b = coef_wing_lat_zone_mu_b$b, mu = coef_wing_lat_zone_mu_b$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 2, add = T)
curve(1-logistic3(x, b = coef_wing_lat_zone_mu_b$b + coef_wing_lat_zone_mu_b$b_SA, 
                  mu = coef_wing_lat_zone_mu_b$mu + coef_wing_lat_zone_mu_b$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 2, add = T)

legend("topright", c("N. America", "S. America"), pch = c(1, 1), col = col_NA_SA_both[c("N. America", "S. America")])

dev.off()

# no wing data:
# overlapping plots:
png("../wing_analysis/plots/fit_wing_A_clines_1_plot_together_no_wing_points.png", height = 6, width = 8, units = "in", res = 300)
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     xlab = "Degrees latitude from equator",
     ylab = "A ancestry proportion")
points(alpha ~ abs_lat, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic3(x, b = coef_lat_zone_mu_b$b + coef_lat_zone_mu_b$b_SA, 
                mu = coef_lat_zone_mu_b$mu + coef_lat_zone_mu_b$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 1, add = T)
curve(logistic3(x, b = coef_lat_zone_mu_b$b, mu = coef_lat_zone_mu_b$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 1, add = T)
curve(0.84*logistic3(x, b = coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA, 
                     mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 3, add = T)
curve(0.84*logistic3(x, b = coef_lat_zone_mu_b_.84$b, mu = coef_lat_zone_mu_b_.84$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 3, add = T)

curve(1-logistic3(x, b = coef_wing_lat_zone_mu_b$b, mu = coef_wing_lat_zone_mu_b$mu), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["N. America"], 
      lwd = 2, lty = 2, add = T)
curve(1-logistic3(x, b = coef_wing_lat_zone_mu_b$b + coef_wing_lat_zone_mu_b$b_SA, 
                  mu = coef_wing_lat_zone_mu_b$mu + coef_wing_lat_zone_mu_b$mu_SA), 
      range(d_A$abs_lat), n = 1000,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 2, add = T)

legend("topright", c("N. America", "S. America"), pch = c(1, 1), col = col_NA_SA_both[c("N. America", "S. America")])
dev.off()

p_A <- ggplot(d_A, aes(x = abs(lat), y = alpha, color = continent)) + 
  geom_point() +
  geom_smooth() +
  ylab("A ancestry proportion")
p_wing <- ggplot(d_A, aes(x = abs(lat), y = -wing_cm*10, color = continent)) + 
  geom_point() +
  geom_smooth() +
  ylab("(negative) Wing length (mm)")
grid.arrange(p_A, p_wing, nrow = 2)
ggsave("../wing_analysis/plots/wing_A_cline.png", 
       plot = grid.arrange(p_A, p_wing, nrow = 2),
       height = 6, width = 8)


# climate vs. latitude model comparison:
# don't allow SA and NA to have different intercepts
# I fit latitude above
# plot fit_lat_.84 predictions:
plot(alpha ~ abs_lat, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by latitude", 
     xlab = "Degrees latitude from the equator",
     ylab = "A ancestry proportion")
points(alpha ~ abs_lat, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x = x, 
                b = coef(fit_lat_.84)[["b"]], 
                mu = coef(fit_lat_.84)[["mu"]],
                K = 0.84),
      range(d_A$abs_lat), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)

# distance from brazil
fit_dist_.84 <- nls_multstart(alpha ~ logistic4(x = km_from_sao_paulo, 
                                               b = b, 
                                               mu = mu,
                                               K = 0.84),
                             start_lower = list(b = -5, mu = min(d_A$km_from_sao_paulo)),
                             start_upper = list(b = 5, mu = max(d_A$km_from_sao_paulo)),
                             supp_errors = 'Y',
                             iter = 250,
                             convergence_count = 100,
                             data = d_A)
summary(fit_dist_.84)
tidy(fit_dist_.84)
glance(fit_dist_.84)
coef_dist_.84 = tidy(fit_dist_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width = abs(4/b))
coef_dist_.84
plot(alpha ~ km_from_sao_paulo, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by dispersal distance", 
     xlab = "Distance from Sao Paulo (km)",
     ylab = "A ancestry proportion")
points(alpha ~ km_from_sao_paulo, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x = x, 
                        b = coef(fit_dist_.84)[["b"]], 
                        mu = coef(fit_dist_.84)[["mu"]],
                        K = 0.84),
      range(d_A$km_from_sao_paulo), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)
# clearly this is a silly exercise

# mean temp
fit_temp_.84 <- nls_multstart(alpha ~ logistic4(x = AnnualMeanTemp, 
                                                b = b, 
                                                mu = mu,
                                                K = 0.84),
                              start_lower = list(b = -5, mu = min(d_A$AnnualMeanTemp)),
                              start_upper = list(b = 5, mu = max(d_A$AnnualMeanTemp)),
                              supp_errors = 'Y',
                              iter = 250,
                              convergence_count = 100,
                              data = d_A)
summary(fit_temp_.84)
tidy(fit_temp_.84)
glance(fit_temp_.84)
coef_temp_.84 = tidy(fit_temp_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width = abs(4/b))
coef_temp_.84
plot(alpha ~ AnnualMeanTemp, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by annual mean temperature", 
     xlab = "Annual mean temperature",
     ylab = "A ancestry proportion")
points(alpha ~ AnnualMeanTemp, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x = x, 
                b = coef(fit_temp_.84)[["b"]], 
                mu = coef(fit_temp_.84)[["mu"]],
                K = 0.84),
      range(d_A$AnnualMeanTemp), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)
# riverside 2014:
#filter(d_A, population == "Riverside_2014") %>%
filter(d_A, year == 2014) %>%
  points(alpha ~ AnnualMeanTemp, data = ., # plot points again on top
       col = "black", cex = 1, pch = 4)



# coldest temp
fit_cold_.84 <- nls_multstart(alpha ~ logistic4(x = MinTempColdestMonth, 
                                                b = b, 
                                                mu = mu,
                                                K = 0.84),
                              start_lower = list(b = -5, mu = min(d_A$MinTempColdestMonth)),
                              start_upper = list(b = 5, mu = max(d_A$MinTempColdestMonth)),
                              supp_errors = 'Y',
                              iter = 250,
                              convergence_count = 100,
                              data = d_A)
summary(fit_cold_.84)
tidy(fit_cold_.84)
glance(fit_cold_.84)
coef_cold_.84 = tidy(fit_cold_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width = abs(4/b))
coef_cold_.84
plot(alpha ~ MinTempColdestMonth, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by minimum temperature", 
     xlab = "Min. temperature coldest month",
     ylab = "A ancestry proportion")
points(alpha ~ MinTempColdestMonth, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x = x, 
                b = coef(fit_cold_.84)[["b"]], 
                mu = coef(fit_cold_.84)[["mu"]],
                K = 0.84),
      range(d_A$MinTempColdestMonth), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)

# winter mean temp
fit_winter_.84 <- nls_multstart(alpha ~ logistic4(x = MeanTempColdestQuarter, 
                                                b = b, 
                                                mu = mu,
                                                K = 0.84),
                              start_lower = list(b = -5, mu = min(d_A$MeanTempColdestQuarter)),
                              start_upper = list(b = 5, mu = max(d_A$MeanTempColdestQuarter)),
                              supp_errors = 'Y',
                              iter = 250,
                              convergence_count = 100,
                              data = d_A)
summary(fit_winter_.84)
tidy(fit_winter_.84)
glance(fit_winter_.84)
coef_winter_.84 = tidy(fit_winter_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width = abs(4/b))
coef_winter_.84
plot(alpha ~ MeanTempColdestQuarter, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by winter temperature", 
     xlab = "Mean temperature coldest quarter",
     ylab = "A ancestry proportion")
points(alpha ~ MeanTempColdestQuarter, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x = x, 
                b = coef(fit_winter_.84)[["b"]], 
                mu = coef(fit_winter_.84)[["mu"]],
                K = 0.84),
      range(d_A$MeanTempColdestQuarter), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)

# annual precipitation:
fit_precip_.84 <- nls_multstart(alpha ~ logistic4(x = AnnualPrecip, 
                                                  b = b, 
                                                  mu = mu,
                                                  K = 0.84),
                                start_lower = list(b = -5, mu = min(d_A$AnnualPrecip)),
                                start_upper = list(b = 5, mu = max(d_A$AnnualPrecip)),
                                supp_errors = 'Y',
                                iter = 250,
                                convergence_count = 100,
                                data = d_A)
summary(fit_precip_.84)
tidy(fit_precip_.84)
glance(fit_precip_.84)
coef_precip_.84 = tidy(fit_precip_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width = abs(4/b))
coef_precip_.84
plot(alpha ~ AnnualPrecip, data = d_A, 
     col = NULL,
     ylim = c(0, 1), 
     main = "African ancestry predicted by annual precipitation", 
     xlab = "Mean annual precipitation (cm)",
     ylab = "A ancestry proportion")
points(alpha ~ AnnualPrecip, data = d_A, # plot points again on top
       col = ifelse(d_A$continent == "S. America", 
                    col_NA_SA_both["S. America"], 
                    col_NA_SA_both["N. America"]))
# plot MAP line from model coefficients:
curve(logistic4(x = x, 
                b = coef(fit_precip_.84)[["b"]], 
                mu = coef(fit_precip_.84)[["mu"]],
                K = 0.84),
      range(d_A$AnnualPrecip), n = 1000,
      col = "black", lwd = 2, add = T, lty = 2)

# compare model fits
models = list(fit_lat_.84, fit_dist_.84, fit_temp_.84, fit_cold_.84, fit_winter_.84, fit_precip_.84)
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
do.call(rbind, lapply(1:length(models), function(i) glance(models[[i]]) %>% mutate(model = model_names[i]))) %>%
  arrange(AIC)

# what is the difference if I exclude outliers?
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
# rerun models all together:
m_all <- lapply(model_vars, function(x) fit_cline_m(var = x, data = d_A))
t_all <- do.call(rbind, lapply(1:length(models), function(i) glance(m_all[[i]]) %>% 
                                 mutate(model = model_names[i]))) %>%
  arrange(AIC) %>%
  left_join(., model_info, by = c("model"="name")) %>%
  mutate(dAIC =  round(AIC - min(AIC), 1),
         rel_lik = exp(- 0.5 * dAIC),
         weight = round(rel_lik / sum(rel_lik), 3)) %>%
  dplyr::select(predictor, df.residual, deviance, dAIC, weight)
t_all
xtable(t_all)
print(xtable(t_all, 
             caption = "\\color{Gray} \textbf{Cline model comparison} Model rankings between logistic cline fits for African ancestry predicted by climate and distance variables (n = 313 bees).",
             label = "AIC_climate_clines",
             type = "latex", 
             latex.environments = NULL), 
      file = "../../bee_manuscript/tables/AIC_cline_fits_climate.tex")


# explore the effects of outliers on these model rankings:
# 2018 data only -- temp clearly best
m_2018only <- lapply(model_vars, function(x) fit_cline_m(var = x, data = filter(d_A, year == 2018)))
do.call(rbind, lapply(1:length(models), function(i) glance(m_2018only[[i]]) %>% mutate(model = model_names[i]))) %>%
  arrange(AIC) %>%
  left_join(., model_info, by = c("model"="name")) %>%
  mutate(dAIC =  round(AIC - min(AIC), 1),
         rel_lik = exp(- 0.5 * dAIC),
         weight = round(rel_lik / sum(rel_lik), 3)) %>%
  dplyr::select(predictor, df.residual, deviance, dAIC, weight)
# just no Riverside 2014 -- temp mostly best
m_noRiv2014 <- lapply(model_vars, function(x) fit_cline_m(var = x, data = filter(d_A, population != "Riverside_2014")))
t_noRiv2014 <- do.call(rbind, lapply(1:length(models), function(i) glance(m_noRiv2014[[i]]) %>% mutate(model = model_names[i]))) %>%
  arrange(AIC) %>%
  left_join(., model_info, by = c("model"="name")) %>%
  mutate(dAIC =  round(AIC - min(AIC), 1),
         rel_lik = exp(- 0.5 * dAIC),
         weight = round(rel_lik / sum(rel_lik), 3)) %>%
  dplyr::select(predictor, df.residual, deviance, dAIC, weight)
t_noRiv2014

# just no skyvalley (desert) -- lat is best
m_noSkyValley <- lapply(model_vars, function(x) fit_cline_m(var = x, data = filter(d_A, !(population == "Riverside_2014" & lat == 33.84))))
do.call(rbind, lapply(1:length(models), function(i) glance(m_noSkyValley[[i]]) %>% mutate(model = model_names[i]))) %>%
  arrange(AIC)
# just no idyllwild (mountain top) -- lat is best
m_noIdyllwild <- lapply(model_vars, function(x) fit_cline_m(var = x, data = filter(d_A, !(population == "Riverside_2014" & lat != 33.84))))
do.call(rbind, lapply(1:length(models), function(i) glance(m_noIdyllwild[[i]]) %>% mutate(model = model_names[i]))) %>%
  arrange(AIC)


AIC(fit_lat_zone_mu_and_b_.84_2018)
AIC(nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                    b = b + b_SA*S_America, 
                                    mu = mu + mu_SA*S_America,
                                    K = 0.84),
                  start_lower = list(b = -5, mu = min(d_A$abs_lat), b_SA = -5, mu_SA = -5),
                  start_upper = list(b = 5, mu = max(d_A$abs_lat), b_SA = 5, mu_SA = 5),
                  supp_errors = 'Y',
                  iter = 250,
                  convergence_count = 100,
                  data = filter(d_A, population != "Riverside_2014")))

fit_lat_2014_2018_.84 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                             b = b + b_2014*from2014, 
                                                             mu = mu + mu_2014*from2014,
                                                             K = 0.84),
                                           start_lower = list(b = -5, mu = min(d_A$abs_lat), b_2014 = -5, mu_2014 = -10),
                                           start_upper = list(b = 5, mu = max(d_A$abs_lat), b_2014 = 5, mu_2014 = 10),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 500,
                                           data = d_A)
                                           #data = filter(d_A, continent == "N. America"))
summary(fit_lat_2014_2018_.84)
tidy(fit_lat_2014_2018_.84)
glance(fit_lat_2014_2018_.84)

coef_lat_2014_2018_.84 = tidy(fit_lat_2014_2018_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_2018 = abs(4/b),
         width_2014 = abs(4/(b + b_2014)))
coef_lat_2014_2018_.84


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
4/quantile(clines$params[clines$params$term == "b", "estimate"], c(0, 0.025, 0.1, 0.5, 0.9, 0.975, 1))
# in km:
111*4/quantile(clines$params[clines$params$term == "b", "estimate"], c(0, 0.025, 0.1, 0.5, 0.9, 0.975, 1))
quantile(clines$params[clines$params$term == "mu", "estimate"], c(0, 0.025, 0.1, 0.5, 0.9, 0.975, 1))
summary(clines$params[clines$params$term == "mu", "estimate"])
diff(quantile(clines$params[clines$params$term == "mu", "estimate"], c(0.025, 0.975))) # 95% of clines have centers within this distance
diff(quantile(clines$params[clines$params$term == "mu", "estimate"], c(0.025, 0.975))) # 95% of clines have centers within this distance
diff(quantile(clines$params[clines$params$term == "mu", "estimate"], c(0.25, 0.75))) # interquartile range


# cline width, under diffusion expected dispersal kernel: w = sqrt(2*pi*t*sigma^2)
mean_width = 4/mean(clines$params[clines$params$term == "b", "estimate"])*111 # km
median_t_admix = median(admix_times[admix_times$ancestry == "A", "time"])
sigma_2 <- mean_width^2/(2*pi*median_t_admix)
sqrt(sigma_2)
summary(abs(rnorm(n = 1000, mean = 0, sd = sqrt(sigma_2))))
sigma_2_years <- mean_width^2/(2*pi*(2018-1957))
sqrt(sigma_2_years)

# add in avg. block size to discussion, and in results range of ind snp cline widths and widths (95% range?)
mean_recombination_rate = 23.94/100 # Jones 2019 recombination rate estimate
# add a plot showing these estimates, some sense of the dist, and also the steepest cline (& other outliers? sure)
10^6/median_t_admix*mean_recombination_rate

# plot density plot to visualize spread of simulated clines
rbind(mutate(mvn_zero$params, data = "MVN_sim_zero_bounded"), 
      mutate(mvn$params, data = "MVN_sim_bounded"), 
      mutate(clines$params, data = "observed")) %>%
  mutate(term = ifelse(term == "b", "w", term),
         estimate = ifelse(term == "w", 4/estimate, estimate)) %>%
  ggplot(., aes(x = estimate, fill = data, color = data)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~term, scales = "free_x") +
  geom_vline(data = data.frame(term = c("w", "mu", "w", "mu"),
                               model = c("0.84 A asym", "0.84 A asym", "0.0 A asym", "0.0 A asym"),
                               estimate = -1*c(4/(coef_lat_zone_mu_b_.84$b_SA + coef_lat_zone_mu_b_.84$b),
                                               coef_lat_zone_mu_b_.84$mu_SA + coef_lat_zone_mu_b_.84$mu,
                                               4/(coef_lat_zone_mu_b_.84$b_SA + coef_lat_zone_mu_b_.84$b),
                                               coef_lat_zone_mu_b_.84$mu_SA + coef_lat_zone_mu_b_.84$mu)),
            aes(xintercept = estimate, linetype = model))
ggsave("plots/compare_ind_snp_to_genomewide_cline_model_estimates.png", height = 3, width = 7.5, units = "in", dpi = 600, device = "png")

# tag outliers:
load("../local_ancestry/results/mean_ancestry_AR_CA.RData")
# violin plot alternative
d_clines_combined <- rbind(mutate(mvn$params, data = "Simulated SNPs (MVN)"), 
      mutate(clines$params, data = "Observed SNPs")) %>%
  mutate(term = ifelse(term == "b", "Width", "Center"),
         estimate = ifelse(term == "Width", 4/estimate, estimate)) %>%
  left_join(., mutate(A_AR_CA, snp_index = 1:nrow(A_AR_CA), data = "Observed SNPs"), by = c("snp_index", "data")) %>%
  mutate(data = factor(data, levels = c("Simulated SNPs (MVN)", "Observed SNPs"), ordered = T)) %>%
  mutate(outlier_status = ifelse(data == "Simulated SNPs (MVN)", "", 
                          ifelse(!is.na(FDR_AR_high), "(High A)", ifelse(!is.na(FDR_AR_low), "(Low A)", "(Non-outliers)"))))

d_clines_combined %>%
  filter(data == "Observed SNPs") %>%
  pivot_wider(data = ., id_cols = c(snp_index, outlier_status), names_from = term, values_from = estimate) %>%
  ggplot(., aes(y = Width, x = Center, color = outlier_status)) +
  geom_point() +
  ylim(0, 15)

violin_ind_clines <- d_clines_combined %>%
  bind_rows(., filter(d_clines_combined, data == "Observed SNPs") %>% mutate(outlier_status = "(All)")) %>%
  mutate(group = paste(data, outlier_status)) %>%
  #group_by(group, data) %>%
  #summarise(n = n())
  #mutate(outlier_status = factor(outlier_status, levels = c("Combined", "High A", "Low A", "Non-outliers"), ordered = T)) %>%
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
        axis.ticks.x=element_blank()) #+
  #geom_point(data = d_clines_combined %>% 
  #             group_by(data, term) %>%
  #             summarise(low = quantile(estimate, 0.025),
   #                      high = quantile(estimate, 0.975),
  #                       mean = mean(estimate)) %>%
  #             pivot_longer(cols = c("low", "high", "mean"), names_to = "summary", values_to = "estimate"),
  #           aes(y = estimate, x = data),
  #           shape = 18, size = 1, color = "white")#  +
  #guides(color = "none") +
  #scale_color_manual(values = c("white", "white", "black"))
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



qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "mu"], clines$params$estimate[clines$params$term == "mu"])
#qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "mu"], mvn$params$estimate[mvn$params$term == "mu"])
abline(a = 0, b = 1, col = "blue")
qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "b"], clines$params$estimate[clines$params$term == "b"])
abline(a = 0, b = 1, col = "blue")
qqplot(mvn_zero$params$estimate[mvn_zero$params$term == "b"], mvn$params$estimate[mvn$params$term == "b"])
abline(a = 0, b = 1, col = "blue")
qqplot(mvn$params$estimate[mvn_zero$params$term == "b"], clines$params$estimate[mvn$params$term == "b"])
abline(a = 0, b = 1, col = "blue")

par(mfrow=c(1,3))
# mu vs. b at the same snps
plot(clines$params$estimate[clines$params$term == "mu"], clines$params$estimate[clines$params$term == "b"], ylim = c(0,1))
# in the neutral model mu vs. b
plot(mvn_zero$params$estimate[clines$params$term == "mu"], mvn_zero$params$estimate[clines$params$term == "b"], ylim = c(0,1))
plot(mvn$params$estimate[clines$params$term == "mu"], mvn$params$estimate[clines$params$term == "b"], ylim = c(0,1))
par(mfrow=c(1,1))
# mu and b across the genome

# what is the FDR?
#FDR_values = c(.1, .05, .01)
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
# visualize tails:
#plot(mu_range[850:1000], fdr_high_mu[850:1000])
#points(mvn_zero$params$estimate[mvn_zero$params$term == "mu"], rep(.01, 100000), col = "blue")

write.table(FDR_clines, "results/FDRs_ind_snp_nls_multistart_clines.txt", quote = F, col.names = T, row.names = F, sep = "\t")



##########3 ---------------------- MAKE SOME PLOTS OF CLINES ---------------- ##############


# plot a set of clines from simulations:
mean_cline_mvn <- group_by(mvn$params, term) %>%
  summarise(mean = mean(estimate))
plot(meta.AR.order.by.lat$lat, 
     rep(0, length(meta.AR.order.by.lat$lat)),
     col = NULL,
     ylim = c(0, 1),
     xlim = range(d_A$lat[d_A$group == "AR_2018"]),
     ylab = "A ancestry",
     xlab = "latitude")
curve(logistic3(mu = mean_cline_mvn$mean[mean_cline_mvn$term == "mu"], 
                b = mean_cline_mvn$mean[mean_cline_mvn$term == "b"], 
                x), 
      from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
      to = range(d_A$lat[d_A$group == "AR_2018"])[2],
      n = 100,
      add = T,
      lwd = 2, 
      lty = 2)
# sample 10 random clines
set.seed(100)
sample_10 <- sample(1:(nrow(mvn$params)/2), 100, replace = F)
sapply(sample_10, function(i)
       curve(logistic3(mu = mvn$params$estimate[mvn$params$term == "mu" & mvn$params$snp_index == i], 
                       b = mvn$params$estimate[mvn$params$term == "b" & mvn$params$snp_index == i], 
                       x), 
             from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
             to = range(d_A$lat[d_A$group == "AR_2018"])[2],
             n = 100,
             col = "grey",
             add = T,
             lwd = 1, 
             lty = 1))

# plot 1 steep cline from my data vs. 1 steep cline in the simulated data and median snp cline
outlier_clines <- c(which.max(clines$params$estimate[clines$params$term == "b"]),
                    which.min(clines$params$estimate[clines$params$term == "b"]),
                    which.max(clines$params$estimate[clines$params$term == "mu"]),
                    which.min(clines$params$estimate[clines$params$term == "mu"]))
sapply(1:length(outlier_clines), function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = rainbow(4)[i],
        add = T,
        lwd = 1, 
        lty = 1))
outlier_mvn_clines <- c(which.max(mvn$params$estimate[mvn$params$term == "b"]),
                    which.min(mvn$params$estimate[mvn$params$term == "b"]),
                    which.max(mvn$params$estimate[mvn$params$term == "mu"]),
                    which.min(mvn$params$estimate[mvn$params$term == "mu"]))
sapply(1:length(outlier_mvn_clines), function(i)
  curve(logistic3(mu = mvn$params$estimate[mvn$params$term == "mu" & mvn$params$snp_index == outlier_mvn_clines[i]], 
                  b = mvn$params$estimate[mvn$params$term == "b" & mvn$params$snp_index == outlier_mvn_clines[i]], 
                  x), 
        from = range(d_A$lat[d_A$group == "AR_2018"])[1], 
        to = range(d_A$lat[d_A$group == "AR_2018"])[2],
        n = 100,
        col = rainbow(4)[i],
        add = T,
        lwd = 1, 
        lty = 2))
plot_random_clines <- function(d, n, seed = 500, color = "darkgrey"){
  # set frame
  plot(abs(meta.AR.order.by.lat$lat), 
       rep(0, length(meta.AR.order.by.lat$lat)),
       col = NULL,
       ylim = c(0, 1),
       xlim = range(d_A$abs_lat),#[d_A$group == "AR_2018"]),
       ylab = "African ancestry frequency",
       xlab = "Degrees latitude from equator")
  # set seed
  set.seed(seed)
  # sample data
  sample_data <- sample(1:(nrow(d$params)/2), n, replace = F)
  # plot
  sapply(sample_data, function(i)
    curve(logistic3(mu = d$params$estimate[d$params$term == "mu" & d$params$snp_index == i], 
                    b = d$params$estimate[d$params$term == "b" & d$params$snp_index == i], 
                    -x), 
          from = range(d_A$abs_lat)[1],
          to = range(d_A$abs_lat)[2],
          #from = range(d_A$abs_lat[d_A$group == "AR_2018"])[1], 
          #to = range(d_A$abs_lat[d_A$group == "AR_2018"])[2],
          n = 100,
          col = "darkgrey",
          add = T,
          lwd = 1, 
          lty = 1))
}
# just the background random clines
png("plots/random_clines_100.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
dev.off()
viridis::magma(n = 3)

png("plots/random_clines_100_plus_steep.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1, function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],

        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1))
dev.off()

png("plots/random_clines_100_plus_introgress.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(2:3, function(i)
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1))
dev.off()

# all 3 outliers
png("plots/random_clines_100_plus_3outliers.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = -clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = -clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
        to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
}
  )
dev.off()

# ADD Points to all of these plots!
png("plots/random_clines_100_plus_steep_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
  })
dev.off()

png("plots/random_clines_100_plus_introgress_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(2:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
       A[outlier_clines[i], meta.AR.order.by.lat$population],
       pch = 20,
       col = viridis::viridis(n = 3)[i])
})
dev.off()

# all 3 outliers
png("plots/random_clines_100_plus_3outliers_points.png", units = "in", res = 300, height = 6, width = 6)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  -x), 
        from = range(d_A$abs_lat)[1], 
        to = range(d_A$abs_lat)[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
}
)
dev.off()


# what does the N. American cline look like for the 'steep' SNP? Also looks steep.
# should exclude Avalon -- such a consistent outlier point can really affect cline shape
png("plots/random_clines_100_plus_3outliers_points_and_NA_at_steep_cline_with_and_without_Avalon.png", units = "in", res = 300, height = 6, width = 6)

plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
sapply(1:3, function(i){
  curve(logistic3(mu = -clines$params$estimate[clines$params$term == "mu" & clines$params$snp_index == outlier_clines[i]], 
                  b = -clines$params$estimate[clines$params$term == "b" & clines$params$snp_index == outlier_clines[i]], 
                  x), 
        from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
        to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::viridis(n = 3)[i],
        add = T,
        lwd = 3, 
        lty = 1)
  points(abs(meta.AR.order.by.lat$lat), 
         A[outlier_clines[i], meta.AR.order.by.lat$population],
         pch = 20,
         col = viridis::viridis(n = 3)[i])
}
)

points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[1], meta.pop$zone == "N. America"],
       pch = 20,
       col = "deeppink")
meta.CA.order.by.lat <- meta.pop %>%
  arrange(lat) %>%
  filter(zone == "N. America")# & population != "Avalon_2014") # exclude Avalon (major outlier) from cline analysis

fit_pink_with_Avalon = nls_multstart(A ~ logistic3(x = lat, 
                                       b = b, mu = mu),
                                       start_lower = list(b = -1, mu = min(d_A$abs_lat)),
                                       start_upper = list(b = 1, mu = max(d_A$abs_lat)),
                                       supp_errors = 'Y',
                                       iter = 250,
                                       convergence_count = 100,
                                       data = data.frame(A = unname(t(A[outlier_clines[1], 
                                                                        meta.CA.order.by.lat$population])),
                                                         lat = meta.CA.order.by.lat$lat, 
                                                         stringsAsFactors = F))
fit_pink = nls_multstart(A ~ logistic3(x = lat, 
                                                   b = b, mu = mu),
                                     start_lower = list(b = -1, mu = min(d_A$abs_lat)),
                                     start_upper = list(b = 1, mu = max(d_A$abs_lat)),
                                     supp_errors = 'Y',
                                     iter = 250,
                                     convergence_count = 100,
                                     data = data.frame(A = unname(t(A[outlier_clines[1], 
                                                                      meta.CA.order.by.lat$population[meta.CA.order.by.lat$population != "Avalon_2014"]])), 
                                                       lat = meta.CA.order.by.lat$lat[meta.CA.order.by.lat$population != "Avalon_2014"], 
                                                       stringsAsFactors = F))

curve(logistic3(mu = unlist(tidy(fit_pink)[tidy(fit_pink)$term == "mu", "estimate"]), 
                b = unlist(tidy(fit_pink)[tidy(fit_pink)$term == "b", "estimate"]), 
                x), 
      from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
      to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
      n = 100,
      col = "deeppink",
      add = T,
      lwd = 3, 
      lty = 1)
curve(logistic3(mu = unlist(tidy(fit_pink_with_Avalon)[tidy(fit_pink_with_Avalon)$term == "mu", "estimate"]), 
                b = unlist(tidy(fit_pink_with_Avalon)[tidy(fit_pink_with_Avalon)$term == "b", "estimate"]), 
                x), 
      from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
      to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
      n = 100,
      col = "deeppink",
      add = T,
      lwd = 3, 
      lty = 2)
dev.off()
# possibly all the clines are a lot steeper without Avalon though .. need to look at the distribution for NA
# overall N. America is going to be a lot more sensitive to drift in 1 pop because it has fewer pops to fit the cline



# I can alternatively plot absolute latitude:
curve(logistic3(mu = -mean_cline_mvn$mean[mean_cline_mvn$term == "mu"], 
                b = -mean_cline_mvn$mean[mean_cline_mvn$term == "b"], 
                x), 
      from = abs(range(d_A$lat[d_A$group == "AR_2018"])[2]), 
      to = abs(range(d_A$lat[d_A$group == "AR_2018"])[1]),
      n = 100,
      col = "blue",
      add = F,
      lwd = 2, 
      lty = 2)

# the top clines with data from N. America too:
meta.pop.A <- meta.pop %>%
  mutate(meanA = apply(A[ , meta.pop$population], 2, mean)) %>%
  mutate(abs_lat = abs(lat))
png("plots/high_A_SA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$lat[meta.pop.A$zone == "S. America"], 
       meta.pop.A$meanA[meta.pop.A$zone == "S. America"],
       ylim = c(0, 1),
       pch = 20,
     cex = 1,
       col = col_blind[1],
        ylab = "A ancestry",
        xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "S. America"], 
       A[outlier_clines[2], meta.pop$zone == "S. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["S. America"])
dev.off()

png("plots/high_A_NA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$abs_lat[meta.pop.A$zone == "N. America"], 
       meta.pop.A$meanA[meta.pop.A$zone == "N. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[2], meta.pop$zone == "N. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["N. America"])
dev.off()

# low SA outlier:
png("plots/low_A_SA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$lat[meta.pop.A$zone == "S. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "S. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "S. America"], 
       A[outlier_clines[3], meta.pop$zone == "S. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["S. America"])
dev.off()

png("plots/low_A_NA_clines.png", units = "in", res = 300, height = 3, width = 3)
plot(meta.pop.A$abs_lat[meta.pop.A$zone == "N. America"], 
     meta.pop.A$meanA[meta.pop.A$zone == "N. America"],
     ylim = c(0, 1),
     pch = 20,
     cex = 1,
     col = col_blind[1],
     ylab = "A ancestry",
     xlab = "latitude")
points(meta.pop$lat[meta.pop$zone == "N. America"], 
       A[outlier_clines[3], meta.pop$zone == "N. America"],
       pch = 20,
       cex = 1,
       col = col_NA_SA_both["N. America"])
dev.off()






##### -------------------- ENRICHMENT OF STEEP CLINES IN REGIONS OF THE GENOME WITH LOW RECOMB? ---------------######


# re-do calculation of enrichment for steep clines in regions of the genome with low recombination rates
load("../local_ancestry/results/sites_r.RData") # load recombination rate data by sites
steepest_clines <- clines$params$estimate[clines$params$term == "b"] > 
                           quantile(clines$params$estimate[clines$params$term == "b"], .99)
steepest_clines05 <- clines$params$estimate[clines$params$term == "b"] > # top 5% 
                           quantile(clines$params$estimate[clines$params$term == "b"], .95)
table(sites_r$r_bin5[steepest_clines])
table(sites_r$r_bin5)
# calc. enrichment of steepest clines in regions of the genome with low recombination rates:
table(sites_r$r_bin5[steepest_clines])/table(sites_r$r_bin5)
table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5)
# (will bootstrap conf. intervals)
# does this pattern hold across all chromosomes? Small # outliers.. this doesn't really make sense to do with small #s
sites_r %>%
  mutate(steepest_clines = steepest_clines,
         steepest_clines05 = steepest_clines05) %>%
  group_by(chr, r_bin5) %>%
  summarise(n = n(),
            steepest_clines05_perc = sum(steepest_clines05)/n,
            steepest_clines01_perc = sum(steepest_clines)/n) %>%
  arrange(r_bin5) %>%
  View(.)
# how are the outliers distributed across chromosomes? Top 5% on all chr. Top 1% on all chr except 16.
sites_r %>%
  mutate(steepest_clines = steepest_clines,
         steepest_clines05 = steepest_clines05) %>%
  group_by(chr) %>%
  summarise(n = n(),
            steepest_clines05 = sum(steepest_clines05),
            steepest_clines01 = sum(steepest_clines))
table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5)

# show across genome:
sites_r %>%
  mutate(., b = clines$params$estimate[clines$params$term == "b"],
         percentile = sapply(1:nrow(sites_r), function(i) 
           ifelse(steepest_clines[i], "1%", 
                  ifelse(steepest_clines05[i], "5%", "n.s")))) %>%
  #filter(chr %in% c("Group1", "Group11")) %>%
  ggplot(., aes(x = pos, y = b, color = percentile)) +
  geom_point() +
  facet_wrap(~chr)

# for bootstrap divide genome into 0.2cM windows (see script divide_map_0.2cM_windows.R):
rmap2 <- read.table("../data/recomb_map/Wallberg_HAv3.1/map_10kb_in_0.2cM_windows.bed", stringsAsFactors = F , 
                    header = T , sep = "\t")
# associate each site with an ancestry call with a 0.2cM window
sites_r_windows <- bedr(
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
head(sites_r_windows)
sites_r_windows = sites_r_windows %>%
  group_by(window_0.2cM) %>%
  summarise(n = n()) %>%
  left_join(sites_r_windows, ., by = "window_0.2cM") %>%
  mutate(b = clines$params$estimate[clines$params$term == "b"]) # add slope

# citation for block bootstrap with genomic data: Bickel, P. J., Boley, N., Brown, J. B., Huang, H., & Zhang, N. R. (2010). Subsampling Methods for Genomic Inference.The Annals ofApplied Statistics,4(4), 1660-1697.http://dx.doi.org/10.1214/10-AOAS363
# this follows ENCODE Strategy..generally though they deal with homogeneous blocks..and it's overly conservative if blocks are not homogeneous...better to segment the blocks in this case.
# other possible publication to look at : https://www.jstor.org/stable/2290993?seq=1#metadata_info_tab_contents

top_1perc = quantile(sites_r_windows$b, 0.99)
top_5perc = quantile(sites_r_windows$b, 0.95)

save(file = "results/sites_r_windows_steep_slopes.RData",
     list = c("sites_r_windows", "top_1perc", "top_5perc"))


# bootstrapping moved to it's own script -- bootstrap_steep_clines.R                                  
#set.seed(200)
#b = bootstrap_steep_clines(windows_with_data = windows_with_data, n_windows = n_windows, name = "test_boot")
#bootstrap_steep_clines(windows_with_data = split(max_clines, max_clines$window_0.2cM), n_windows = nrow(max_clines), name = "test_boot")

# can alternatively bootstrap using just the steepest cline for any 0.2cM window:
# or one random cline for a 0.2cM window.
#max_clines <- sites_r_windows %>%
#  mutate(b = clines$params$estimate[clines$params$term == "b"]) %>% # add slope
#  group_by(., sites_r_windows$window_0.2cM) %>%
#  arrange(., sample(1:nrow(.), replace = F)) %>% # randomly permute order
#  arrange(., desc(b)) %>% # move highest b (steepest slope) to top
#  filter(!duplicated(window_0.2cM)) # keep the first site with the highest slope
#max_clines %>%
#  group_by(., r_bin5) %>%
#  summarise(perc1 = sum(b > quantile(max_clines$b, .99))/n(),
#            perc5 = sum(b > quantile(max_clines$b, .95))/n(),
#            n = n())

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
# # estimate lowest recomb. bin:
lowr_estimate = (table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5))[[1]]
lowr_estimate
# basic bootstrap conf. intervals
2*lowr_estimate - quantile(boot$top5_sample_lowr_perc, c(.975, .025))
lowr_estimate - quantile(boot$top5_sample_lowr_perc - lowr_estimate, c(0.975, 0.025)) # same

# plot bootstrap lowr
hist(boot$top5_sample_lowr_perc, main = "low r bootstrap - % steep clines")
abline(v = lowr_estimate, col = "orange")
abline(v = 2*lowr_estimate - quantile(boot$top5_sample_lowr_perc, c(.975, .025)), col = "blue")
abline(v = 0.05, col = "red")

# estimate highest recomb. bin:
highr_estimate = (table(sites_r$r_bin5[steepest_clines05])/table(sites_r$r_bin5))[[5]]
highr_estimate
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


plot(sites_r$cM_Mb, clines$params$estimate[clines$params$term == "b"], main = "slope b by r (cM/Mb)")
non_outlier_cline_center = abs(clines$params$estimate[clines$params$term == "mu"] - mean(clines$params$estimate[clines$params$term == "mu"])) < 1
plot(sites_r$cM_Mb[non_outlier_cline_center], 
     clines$params$estimate[clines$params$term == "b"][non_outlier_cline_center], 
     main = "non-outlier centers: slope b by r (cM/Mb)")
# % of overall points with ancestry calls in lowest recombination bin vs. % steepest clines in lowest recombination bin
table(sites_r$r_bin5)/nrow(sites_r)
table(sites_r$r_bin5[steepest_clines])/nrow(sites_r[steepest_clines,])
table(sites_r$r_bin5[steepest_clines05])/nrow(sites_r[steepest_clines05,])
table(sites_r$r_bin5[steepest_clines05 & non_outlier_cline_center])/table(sites_r$r_bin5[non_outlier_cline_center])


mutate(sites_r, b = clines$params$estimate[clines$params$term == "b"]) %>%
  group_by(r_bin5) %>%
  summarise(mean_b = mean(b),
            mean_w = mean(4/b))
mutate(sites_r, b = clines$params$estimate[clines$params$term == "b"]) %>%
  mutate(w = 4/b) %>%
  lm(data = ., w ~ r_bin5) %>%
  summary()
mutate(sites_r, b = clines$params$estimate[clines$params$term == "b"]) %>%
  lm(data = ., b ~ r_bin5) %>%
  anova()
mutate(sites_r, b = clines$params$estimate[clines$params$term == "b"]) %>% # still true if I exclude clines with somewhat skewed centers
  filter(non_outlier_cline_center) %>%
  group_by(r_bin5) %>%
  summarise(mean = mean(b))

cbind(sites_r, clines$params[clines$params$term == "b", ]) %>%
  filter(estimate >= quantile(estimate, .99)) %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = estimate, color = r_bin5_factor)) +
  facet_wrap(~chr)

# (!) Here high b is steep, but for slopes the other way (on + lat), high b will be large & negative
clines_plus_sites <- sites_r %>%
  mutate(b = clines$params[clines$params$term == "b", "estimate"]) %>%
  mutate(mu = clines$params[clines$params$term == "mu", "estimate"]) %>%
  mutate(top_b = b >= quantile(b, .99)) %>%
  mutate(top_mu = mu >= quantile(mu, .99)) %>%
  mutate(bottom_mu = mu <= quantile(mu, .01))
mean_b_by_r <- clines_plus_sites %>%
  group_by(r_bin5_factor) %>%
  summarise(mean_b = mean(b),
            mean_mu = mean(mu))
clines_plus_sites %>%
  group_by(r_bin5_factor) %>%
  summarise(mean_b = mean(mu))
clines_plus_sites %>%
  group_by(r_bin5_factor) %>%
  summarise(perc_top_b = mean(top_b),
            perc_top_mu = mean(top_mu),
            perc_bottom_mu = mean(bottom_mu))
b_by_r <- with(clines_plus_sites, lm(b ~ cM_Mb))
summary(b_by_r)
mu_by_r <- with(clines_plus_sites, lm(mu ~ cM_Mb))
summary(mu_by_r)
# subsample to plot faster
dplyr::sample_frac(clines_plus_sites, size = .1) %>%
  ggplot(., aes(x = cM_Mb, y = b)) +
  geom_point() +
  geom_smooth()
# ok so the slopes go down with higher recomb, but is this a big enough effect ot be interesting? not sure. let's plot:
#mean_b_by_r
png("plots/lower_r_regions_have_very_slightly_steeper_mean_slopes.png", height = 6, width = 6, units = "in", res = 300)
plot_random_clines(d = clines, n = 100, seed = 500, color = "darkgrey")
for (i in 1:5){
  curve(logistic3(mu = -mean_b_by_r$mean_mu[i], 
                  #mu = -mean(clines_plus_sites$mu), 
                  b = -mean_b_by_r$mean_b[i], 
                  x), 
        from = range(abs(d_A$lat)[d_A$group == "AR_2018"])[1], 
        to = range(abs(d_A$lat)[d_A$group == "AR_2018"])[2],
        n = 100,
        col = viridis::magma(5)[i],
        add = T,
        lwd = 2, 
        lty = 1)
}
legend(x = "topright", title = "cM/Mb", legend = levels(sites_r$r_bin5_factor), 
       cex = 0.5,
       col = viridis::magma(5), lty = 1, lwd = 3)
dev.off()
# plot with ggplot if needed later

mean_b_by_r$mean_mu

# Group11 possibly has some large inversion or something, not the largest peak, but the widest.
# also not colocalized with the region of high M
cbind(sites_r, clines$params[clines$params$term == "b", ]) %>%
  mutate(AR = apply(A[ , meta.AR.order.by.lat$population], 1, mean)) %>%
  filter(estimate >= quantile(estimate, .99)) %>%
  ggplot(.) +
  geom_point(aes(x = pos, y = AR, color = r_bin5_factor)) +
  facet_wrap(~chr)

load("../local_ancestry/results/A.RData")
A_mean <- apply(A, 2, mean)
A_CA <- apply(A[ , meta.pop$population[meta.pop$zone == "N. America"]], 1, mean)
A_CA_2018 <- apply(A[ , meta.pop$population[meta.pop$year == 2018 & meta.pop$zone == "N. America"]], 1, mean)
A_CA_most_northern <- apply(A[ , meta.pop$population[meta.pop$lat >= 35 & meta.pop$zone == "N. America"]], 1, mean)
summary(A_CA)
summary(A_CA_2018)
summary(A_CA_most_northern)
# low in top 1/4 AR and CA?
A_CA_top <- apply(A[ , meta.pop$population[meta.pop$lat >= twentyfive_perc_lat["CA"] & meta.pop$zone == "N. America"]], 1, mean)
A_AR_top <- apply(A[ , meta.pop$population[abs(meta.pop$lat) >= twentyfive_perc_lat["AR"] & meta.pop$zone == "S. America"]], 1, mean)
plot(A_CA_top, A_AR_top, pch = 20, col = alpha("black", 0.1), xlim = c(0,1), ylim = c(0,1))
abline(1, 0, col = "blue")
abline(v = 0, col = "red")
abline(h = 0, col = "orange")
A_AR_tophalf <- apply(A[ , meta.pop$population[abs(meta.pop$lat) >= fifty_perc_lat["AR"] & meta.pop$zone == "S. America"]], 1, mean)
plot(A_CA, A_AR_tophalf, pch = 20, col = alpha("black", 0.1), xlim = c(0,1), ylim = c(0,1), main = "all CA vs. lower half Arg")
abline(0, 1, col = "blue")
abline(v = 0, col = "red")
abline(h = 0, col = "orange")
summary(A_CA_top)
summary(A_AR_top)
summary(A_AR_tophalf)


# more complicated cline with 2 exponential tails:
cline_barton_szymura <- nls_multstart(alpha ~ ifelse(abs_lat > 35, 0.1, 
                                                     logistic4(x = abs_lat,
                                                         b = b,
                                                         mu = mu,
                                                         K = 0.84)), 
    data = filter(d_A, continent == "S. America"), 
    start_lower = list(b = -1, mu = min(d_A$abs_lat)),
    start_upper = list(b = 0, mu = max(d_A$abs_lat)),
    iter = 100)


curve(szymura_barton_center(x, y = 0, w = 10), from = -10, to = 10)
curve(1-szymura_barton_edge(x, y = 0, theta = 1, w = 10), from = 0, to = 10, col = "blue", add = T)
curve(szymura_barton_edge(x, y = 0, theta = 1, w = -10), from = -10, to = 0, col = "red", add = T)
curve(szymura_barton_center(x, y = -30, w = 7), from = -40, to = -25)
curve(1 - szymura_barton_edge(x, y = -30, theta = 1, w = 7), from = -30, to = -25, col = "blue", add = T)
curve(szymura_barton_edge(x, y = -30, theta = 1, w = -7), from = -40, to = -30, col = "red", add = T)
# closer to scenario here:
curve(szymura_barton_center(x, y = -30, w = 7), from = -40, to = -25)
curve(1 - szymura_barton_edge(x, y = -30, theta = 1, w = 7), from = -30, to = -25, col = "blue", add = T)
curve(szymura_barton_edge(x, y = -30, theta = 1, w = -7), from = -40, to = -30, col = "red", add = T)


curve(stepped_cline_center(x, y = 0, w = 10) - .5, from = -10, to = 10)
curve(szymura_barton_center(x, y = 0, w = 10) - .5, from = -10, to = 10, add = T, col = "purple", lty = 2)
abline(a = 0, b = 1/10, col = "blue") # just to confirm expected slope

# draw clines with tails
curve(stepped_cline_center(x, y = -30, w = 7), from = -40, to = -25, col = "purple")
curve(szymura_barton_center(x, y = -30, w = 7), from = -40, to = -25, add = T, col = "white", lty = 2)
#curve(szymura_barton_edge(x, y = -30, theta = 1, w = -7), from = -40, to = -30, col = "red", add = T)
curve(stepped_cline_edge(x, y = -30, theta = 1/2, w = 7, d = -2), #left tail
      from = -40, to = -32, col = "red", lty = 2, add = T)
curve(1 - stepped_cline_edge(x, y = -30, theta = 1/2, w = -7, d = 2), # right tail
      from = -28, to = -25, col = "blue", lty = 2, add = T)

# try to fit 3-part cline model:
# full cline:
# symmetric:
fit_stepped_cline_symmetric <- function(data, rescale = 0.84){
  cline_stepped <- nls_multstart(alpha_rescaled ~ ifelse(lat > y + d,
                                                     1 - stepped_cline_edge(x = lat, 
                                                                            y = y, 
                                                                            theta = theta, 
                                                                            w = -w, 
                                                                            d = d),
                                                     ifelse(lat < y - d,
                                                            stepped_cline_edge(x = lat, 
                                                                               y = y, 
                                                                               theta = theta, 
                                                                               w = w, 
                                                                               d = -d),
                                                            stepped_cline_center(x = lat, 
                                                                                 y = y, 
                                                                                 w = w))), 
                                 data =  data %>%
                                   mutate(alpha_rescaled = alpha/rescale), 
                                 lower = c(y = -40, d = 0, theta = 0, w = 0),
                                 upper = c(y = -25, d = 5, theta = 1, w = 15),
                                 start_lower = list(y = -40, d = 0, theta = 0, w = 0),
                                 start_upper = list(y = -25, d = 5, theta = 1, w = 15), # width here is in degrees latitude
                                 iter = 100)
  return(cline_stepped)
}

# constrain theta = 1
fit_stepped_cline_symmetric_theta1 <- function(d_thresh, data, rescale = 0.84){
  cline_stepped <- nls_multstart(alpha_rescaled ~ ifelse(lat > y + d,
                                                     1 - stepped_cline_edge(x = lat, 
                                                                            y = y, 
                                                                            theta = 1, 
                                                                            w = -w, 
                                                                            d = d),
                                                     ifelse(lat < y - d,
                                                            stepped_cline_edge(x = lat, 
                                                                               y = y, 
                                                                               theta = 1, 
                                                                               w = w, 
                                                                               d = -d),
                                                            stepped_cline_center(x = lat, 
                                                                                 y = y, 
                                                                                 w = w))), 
                                 data =  data %>%
                                   mutate(alpha_rescaled = alpha/rescale), 
                                 lower = c(y = -40, d = 0, w = 0),
                                 upper = c(y = -25, d = 5, w = 15),
                                 start_lower = list(y = -40, d = 0, w = 0),
                                 start_upper = list(y = -25, d = 5, w = 15), # width here is in degrees latitude
                                 iter = 100)
  return(cline_stepped)
}

# non-symmetric
fit_stepped_cline <- function(dL, dR, data, rescale = 0.84){
  cline_stepped <- nls_multstart(alpha_rescaled ~ ifelse(lat > y + dR,
                                                     1 - stepped_cline_edge(x = lat, 
                                                                            y = y, 
                                                                            theta = thetaR, 
                                                                            w = -w, 
                                                                            d = dR),
                                                     ifelse(lat < y - dL,
                                                            stepped_cline_edge(x = lat, 
                                                                               y = y, 
                                                                               theta = thetaL, 
                                                                               w = w, 
                                                                               d = -dL),
                                                            stepped_cline_center(x = lat, 
                                                                                 y = y, 
                                                                                 w = w))), 
                                 data =  data %>%
                                   mutate(dL = dL) %>%
                                   mutate(dR = dR) %>%
                                   mutate(alpha_rescaled = alpha/rescale), 
                                 lower = c(y = -40, thetaR = 0, w = 0, thetaL = 0),
                                 upper = c(y = -25, thetaR = 1, w = 15, thetaL = 1),
                                 start_lower = list(y = -40, thetaR = 0, w = 0, thetaL = 0),
                                 start_upper = list(y = -25, thetaR = 1, w = 15, thetaL = 1), # width here is in degrees latitude
                                 iter = 100)
  return(cline_stepped)
}

# fit d also:
fit_stepped_cline_d <- function(data, rescale = 0.84, supp_errors = 'N'){
  cline_stepped <- nls_multstart(alpha_rescaled ~ ifelse(lat > y + dR,
                                                         1 - stepped_cline_edge(x = lat, 
                                                                                y = y, 
                                                                                theta = thetaR, 
                                                                                w = -w, 
                                                                                d = dR),
                                                         ifelse(lat < y - dL,
                                                                stepped_cline_edge(x = lat, 
                                                                                   y = y, 
                                                                                   theta = thetaL, 
                                                                                   w = w, 
                                                                                   d = -dL),
                                                                stepped_cline_center(x = lat, 
                                                                                     y = y, 
                                                                                     w = w))), 
                                 data =  data %>%
                                   mutate(alpha_rescaled = alpha/rescale), 
                                 lower = c(y = -40, dR = 0, thetaR = 0, w = 0, dL = 0, thetaL = 0),
                                 upper = c(y = -25, dR = 5, thetaR = 1, w = 15, dL = 5, thetaL = 1),
                                 start_lower = list(y = -40, dR = 0, thetaR = 0, w = 0, dL = 0, thetaL = 0),
                                 start_upper = list(y = -25, dR = 5, thetaR = 1, w = 15, dL = 5, thetaL = 1), # width here is in degrees latitude
                                 iter = 100,
                                 supp_errors = supp_errors)
  return(cline_stepped)
}



# simple center cline only
fit_stepped_cline_center <- function(data, rescale = 0.84){
  stepped_cline <- nls_multstart(alpha_rescaled ~ stepped_cline_center(x = lat, 
                                                                        y = y, 
                                                                        w = w), 
                                 data =  data %>%
                                   mutate(alpha_rescaled = alpha/rescale), 
                                 lower = c(y = -40, w = 0),
                                 upper = c(y = -25, w = 40),
                                 start_lower = list(y = -40, w = 0),
                                 start_upper = list(y = -25, w = 15), # width here is in degrees latitude
                                 iter = 100)
  return(stepped_cline)
}

d_A_SA <- filter(d_A, continent == "S. America")
cline_stepped_center <- fit_stepped_cline_center(data = d_A_SA)
cline_stepped_center_1 <- fit_stepped_cline_center(data = d_A_SA, rescale = 1)
cline_stepped_symmetric <- fit_stepped_cline_symmetric(data = d_A_SA)
cline_stepped_symmetric_1 <- fit_stepped_cline_symmetric(data = d_A_SA, rescale = 1)
cline_stepped_symmetric_theta1 <- fit_stepped_cline_symmetric_theta1(data = d_A_SA)
cline_stepped <- fit_stepped_cline(dL = 2, dR = 2, data = d_A_SA)
cline_stepped_1 <- fit_stepped_cline(dL = 2, dR = 2, data = d_A_SA, rescale = 1)
cline_stepped_d <- fit_stepped_cline_d(data = d_A_SA, rescale = 0.84)
cline_stepped_d_1 <- fit_stepped_cline_d(data = d_A_SA, rescale = 1)
# only fit middle of the data:
cline_stepped_center_truncated <- nls_multstart(alpha_0.84 ~ stepped_cline_center(x = lat, 
                                                                                  y = y, 
                                                                                  w = w), 
                                                data =  filter(d_A, continent == "S. America") %>%
                                                  mutate(alpha_0.84 = alpha/0.84) %>%
                                                  filter(lat > -35 & lat < -30), 
                                                lower = c(y = -40, w = 0),
                                                upper = c(y = -25, w = 40),
                                                start_lower = list(y = -40, w = 0),
                                                start_upper = list(y = -25, w = 15), # width here is in degrees latitude
                                                iter = 100)


coefficients(cline_stepped_center)
coefficients(cline_stepped_center_1)
coefficients(cline_stepped_symmetric_theta1)
coefficients(cline_stepped_symmetric)
coefficients(cline_stepped)
coefficients(cline_stepped_1)
coefficients(cline_stepped_center_truncated)
coefficients(cline_stepped_d)
coefficients(cline_stepped_d_1)

# plot residuals of cline_stepped_center
with(d_A_SA, plot(lat, alpha - 0.84*predict(cline_stepped_symmetric),
     ylab = "residuals", main = "stepped cline symmetric"))
with(d_A_SA, plot(lat, alpha - 0.84*predict(cline_stepped_d),
                  ylab = "residuals", main = "stepped cline unconstrained"))
with(d_A_SA, plot(lat, alpha - 0.84*predict(cline_stepped_center),
                  ylab = "residuals", main = "center cline"))
plot(nlsResiduals(cline_stepped_d))
test.nlsResiduals(nlsResiduals(cline_stepped_d))
plot(nlsResiduals(cline_stepped_symmetric))
test.nlsResiduals(nlsResiduals(cline_stepped_symmetric))
plot(nlsResiduals(cline_stepped_center))
test.nlsResiduals(nlsResiduals(cline_stepped_center))
cont <- nlstools::nlsContourRSS(cline_stepped_d)
# test model fits:
# cline with two different tails fits the best. BUT the right 'tail' starts at the center of the cline, which is nonstandard and complicates interpretation.
anova(cline_stepped_d, cline_stepped_center, test = "F") # sig.
anova(cline_stepped_symmetric, cline_stepped_center, test = "F") # not sig.
anova(cline_stepped_symmetric, cline_stepped_d, test = "F") # sig.

AIC(cline_stepped_d, cline_stepped_d_1, cline_stepped_symmetric, cline_stepped_symmetric_1, cline_stepped_center, cline_stepped_center_1)
# plot cline results:
png(filename = "plots/SA_stepped_cline_fits.png", height = 5, width = 5.2, units = "in", res = 600)
with(filter(d_A, continent == "S. America"), 
     plot(lat, alpha,
          xlim = c(-37, -27), ylim = c(0,1),
          main = "S. American genomewide cline",
          ylab = "A ancestry proportion"))
# 1 cline
curve(0.84*stepped_cline_center(x, 
                                y = coefficients(cline_stepped_center)["y"], 
                                w = coefficients(cline_stepped_center)["w"]), 
      from = -40, to = -25, 
      col = "purple", lwd = 2, add = T)
# stepped cline
# symmetric
curve(0.84*ifelse(x > coefficients(cline_stepped_symmetric)["y"] + coefficients(cline_stepped_symmetric)["d"],
                  1 - stepped_cline_edge(x, 
                                         y = coefficients(cline_stepped_symmetric)["y"], 
                                         theta = coefficients(cline_stepped_symmetric)["theta"], 
                                         w = -coefficients(cline_stepped_symmetric)["w"], 
                                         d = coefficients(cline_stepped_symmetric)["d"]),
                  ifelse(x < coefficients(cline_stepped_symmetric)["y"] - coefficients(cline_stepped_symmetric)["d"],
                         stepped_cline_edge(x, 
                                            y = coefficients(cline_stepped_symmetric)["y"], 
                                            theta = coefficients(cline_stepped_symmetric)["theta"], 
                                            w = coefficients(cline_stepped_symmetric)["w"], 
                                            d = -coefficients(cline_stepped_symmetric)["d"]),
                         stepped_cline_center(x = x, 
                                              y = coefficients(cline_stepped_symmetric)["y"], 
                                              w = coefficients(cline_stepped_symmetric)["w"]))),
      from = -40, to = -25, add = T, col = "green", lty = 1, lwd = 2)
curve(0.84*stepped_cline_center(x = x, 
                                y = coefficients(cline_stepped_symmetric)["y"], 
                                w = coefficients(cline_stepped_symmetric)["w"]),
      from = -40, to = -25, add = T, col = "green", lty = 2, lwd = 2)
# free tails annealed stepped cline:
curve(0.84*ifelse(x > coefficients(cline_stepped_d)["y"] + coefficients(cline_stepped_d)["dR"],
                  1 - stepped_cline_edge(x, 
                                         y = coefficients(cline_stepped_d)["y"], 
                                         theta = coefficients(cline_stepped_d)["thetaR"], 
                                         w = -coefficients(cline_stepped_d)["w"], 
                                         d = coefficients(cline_stepped_d)["dR"]),
                  ifelse(x < coefficients(cline_stepped_d)["y"] - coefficients(cline_stepped_d)["dL"],
                         stepped_cline_edge(x, 
                                            y = coefficients(cline_stepped_d)["y"], 
                                            theta = coefficients(cline_stepped_d)["thetaL"], 
                                            w = coefficients(cline_stepped_d)["w"], 
                                            d = -coefficients(cline_stepped_d)["dL"]),
                         stepped_cline_center(x = x, 
                                              y = coefficients(cline_stepped_d)["y"], 
                                              w = coefficients(cline_stepped_d)["w"]))),
      from = -40, to = -25, add = T, col = "blue", lty = 1, lwd = 2)

if (FALSE){
  curve(ifelse(x > coefficients(cline_stepped_d_1)["y"] + coefficients(cline_stepped_d_1)["dR"],
               1 - stepped_cline_edge(x, 
                                      y = coefficients(cline_stepped_d_1)["y"], 
                                      theta = coefficients(cline_stepped_d_1)["thetaR"], 
                                      w = -coefficients(cline_stepped_d_1)["w"], 
                                      d = coefficients(cline_stepped_d_1)["dR"]),
               ifelse(x < coefficients(cline_stepped_d_1)["y"] - coefficients(cline_stepped_d_1)["dL"],
                      stepped_cline_edge(x, 
                                         y = coefficients(cline_stepped_d_1)["y"], 
                                         theta = coefficients(cline_stepped_d_1)["thetaL"], 
                                         w = coefficients(cline_stepped_d_1)["w"], 
                                         d = -coefficients(cline_stepped_d_1)["dL"]),
                      stepped_cline_center(x = x, 
                                           y = coefficients(cline_stepped_d_1)["y"], 
                                           w = coefficients(cline_stepped_d_1)["w"]))),
        from = -40, to = -25, add = T, col = "orange", lty = 1, lwd = 2)
}


# continue middle cline for comparison for whole range
curve(0.84*stepped_cline_center(x = x, 
                                y = coefficients(cline_stepped_d)["y"], 
                                w = coefficients(cline_stepped_d)["w"]),
      from = -40, to = -25, add = T, col = "blue", lty = 2, lwd = 2)

#curve(logistic4(x = -x, # same curve as stepped_cline_center
#                mu = -coefficients(cline_stepped_d)["y"], 
#                b = -4/coefficients(cline_stepped_d)["w"],
#                K = 0.84),
#      from = -40, to = -25, add = T, col = "black", lty = 2, lwd = 2)
legend("topleft", c("simple logistic cline", "barton 'stepped' cline w/ symmetric tails", 
                    "barton 'stepped' cline w/ free tails", "center of cline extended"),
       col = c("purple", "green", "blue", "black"), lty = c(1,1,1,2), cex = .5)
dev.off()



coefficients(cline_stepped_d) # theta can be interpreted as s_e/s*, or selection directly on a locus/effective selection (e.g. due to ld with selected loci)
# here the right tail has theta = 0.77 but the left tail has theta = 0
nlstools::confint2(cline_stepped_d, method = "profile")
nlstools::confint2(cline_stepped_d)
nlstools::confint2(cline_stepped_center, method = "profile")
nlstools::confint2(cline_stepped_center)
nlstools::confint2(cline_stepped_symmetric, method = "profile")
plot(profile(cline_stepped_center))
plot(profile(cline_stepped_d)) # hmm...singular gradient.

# simple center fit only cline:
cline_barton_szymura_center <- nls_multstart(alpha/0.84 ~ szymura_barton_center(x = lat, 
                                                                                y = y, 
                                                                                w = w), 
                                             data = filter(d_A, continent == "S. America"), 
                                             lower = c(y = -40, w = 0),
                                             upper = c(y = -25, w = 15),
                                             start_lower = list(y = -40, w = 0),
                                             start_upper = list(y = -25, w = 15), # width here is in degrees latitude
                                             iter = 100)
# full cline:
cline_barton_szymura <- nls_multstart(alpha/0.84 ~ ifelse(lat > -31 + 1, 
                                                           1 - szymura_barton_edge(x = lat, 
                                                                               y = y,
                                                                               theta = theta, 
                                                                               w = w),
                                                           ifelse(lat < -31 - 1,
                                                                  szymura_barton_edge(x = lat,
                                                                                      y = y,
                                                                                      theta = theta, 
                                                                                      w = -w),
                                                                  szymura_barton_center(x = lat, 
                                                                                        y = y, 
                                                                                        w = w))), 
                                      data = filter(d_A, continent == "S. America"), 
                                      lower = c(y = -40, w = 0, theta = 0),
                                      upper = c(y = -25, w = 15, theta = 1),
                                      start_lower = list(y = -40, w = 0, theta = 0),
                                      start_upper = list(y = -25, w = 15, theta = 1), # width here is in degrees latitude
                                      iter = 100)

with(filter(d_A, continent == "S. America"), plot(alpha/0.84 ~ lat))
curve(szymura_barton_center(x, 
                                y = coefficients(cline_barton_szymura)["y"], 
                                w = coefficients(cline_barton_szymura)["w"]), 
      from = -40, to = -25, col = "blue", add = T)
curve(1 - szymura_barton_edge(x, 
                              y = coefficients(cline_barton_szymura)["y"],
                              theta = coefficients(cline_barton_szymura["theta"]), 
                              w = coefficients(cline_barton_szymura)["w"]), 
      from = -30, to = -25, col = "orange", add = T)
curve(szymura_barton_edge(x, 
                          y = coefficients(cline_barton_szymura)["y"],
                          theta = coefficients(cline_barton_szymura["theta"]), 
                          w = -coefficients(cline_barton_szymura)["w"]), 
      from = -40, to = -34, col = "red", add = T)


# add asymptote free parameter
# why K?
# because bees in brazil have ~ 84% A ancestry, the true
# asymptote is 84% A not 100% A for the cline, so
# I can optionally use this to rescale the inferred logistic curves


curve(logistic4(x, mu = 32, b = -0.6, K = 0.84), from = 25, to = 40)
curve(0.84*(1 - logistic4(x, mu = 32, b = 0.6, K = 0.84)/0.84), from = 25, to = 40, ylim = 0:1)
curve(logistic4(x, mu = 32, b = -0.6, K = 0.84), from = 25, to = 40, add = T, col = "blue", lty = 2)
curve(logistic4(-x, mu = -32, b = 0.6, K = 0.84), from = 25, to = 40, add = T, col = "orange", lty = 2)

curve(logistic4(x, mu = -32, b = 0.6, K = 0.84), from = -40, to = -25, col = "blue", lty = 2)

# these are the same:
curve(logistic4(x, mu = 32, b = -0.6, K = 1), from = 25, to = 40, col = "purple")
curve(1 - logistic4(x, mu = 32, b = 0.6, K = 1), from = 25, to = 40, col = "orange", lty = 2, add = T)
# parameterization equivalence gets a little messier with the 0.84 scaling:
curve(logistic4(x, mu = 32, b = -0.6, K = 0.84), from = 25, to = 40, col = "blue", add = T)
curve(0.84*(1 - logistic4(x, mu = 32, b = 0.6, K = 0.84)/0.84), from = 25, to = 40, col = "green", add = T, lty = 2)


# plot bootstrapped clines:
boots_stepped <- read.table("results/bootstrap_stepped_cline_SA_seed100.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

nrow(boots_stepped) # most, but not all, bootstraps could be fit to the model (missing = no fit)
table(boots_stepped$isConv) # plus one more didn't converge
boots_stepped <- filter(boots_stepped, isConv) # filter for just converged
png("results/bootstrap_results_SA_stepped_cline.png", height = 5.2, width = 5.2, units = "in", dpi = 600)
# note this is just the raw bootstraps, not a CI
with(filter(d_A, continent == "S. America"), 
     plot(lat, alpha,
          xlim = c(-37, -27), ylim = c(0,1),
          main = "S. American genomewide cline",
          ylab = "A ancestry proportion"))
apply(boots_stepped, 1, function(boot){
  curve(stepped_cline_3parts(x = x, 
                             y = boot["y"], 
                             w = boot["w"], 
                             thetaR = boot["thetaR"], 
                             thetaL = boot["thetaL"], 
                             dR = boot["dR"], 
                             dL = boot["dL"], 
                             rescale = 0.84),
        from = -40, to = -25, add = T, col = alpha("black", 0.01), lty = 1, lwd = 2)
})
curve(stepped_cline_3parts(x = x, 
                           y = coefficients(cline_stepped_d)["y"], 
                           w = coefficients(cline_stepped_d)["w"], 
                           thetaR = coefficients(cline_stepped_d)["thetaR"], 
                           thetaL = coefficients(cline_stepped_d)["thetaL"], 
                           dR = coefficients(cline_stepped_d)["dR"], 
                           dL = coefficients(cline_stepped_d)["dL"], 
                           rescale = 0.84),
      from = -40, to = -25, add = T, col = "blue", lty = 1, lwd = 2)
curve(logistic4(-x, b = coef_lat_zone_mu_b_.84$b + coef_lat_zone_mu_b_.84$b_SA, 
                mu = coef_lat_zone_mu_b_.84$mu + coef_lat_zone_mu_b_.84$mu_SA,
                K = 0.84), 
      from = -40, to = -25,
      col = col_NA_SA_both["S. America"], 
      lwd = 2, lty = 2, add = T)
curve(0.84*stepped_cline_center(x = x, # same
                                y = coefficients(cline_stepped_center)["y"], 
                                w = coefficients(cline_stepped_center)["w"]),
      from = -40, to = -25, add = T, col = col_NA_SA_both["S. America"], lty = 3, lwd = 2)
dev.off()
# calculate pivot confidence interval:
coefficients(cline_stepped_d)
#ci_boots_stepped = estimate - quantile(boots - estimate, c(0.975, 0.025))



# exclude all 2014 data -- we get a very similar answer about parallel clines
fit_lat_zone_mu_and_b_.84_2018 <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                             b = b + b_SA*S_America, 
                                                             mu = mu + mu_SA*S_America,
                                                             K = 0.84),
                                           start_lower = list(b = -5, mu = min(d_A$abs_lat), b_SA = -5, mu_SA = -10),
                                           start_upper = list(b = 5, mu = max(d_A$abs_lat), b_SA = 5, mu_SA = 10),
                                           supp_errors = 'Y',
                                           iter = 250,
                                           convergence_count = 100,
                                           data = filter(d_A, year == 2018))
summary(fit_lat_zone_mu_and_b_.84_2018)
coef_lat_zone_mu_b_.84_2018 = tidy(fit_lat_zone_mu_and_b_.84_2018) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_NA = abs(4/b),
         width_SA = abs(4/(b + b_SA)))
coef_lat_zone_mu_b_.84
coef_lat_zone_mu_b
coef_lat_zone_mu_b_.84_2018


