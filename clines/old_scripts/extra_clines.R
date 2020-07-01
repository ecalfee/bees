# load data from plot_clines.R. Removed this old code to de-clutter.
library(betareg) # alternative ML fitting
#visualize individual model fits for each latitude/climate variable:

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
# put all these models in one object
models = list(fit_lat_.84, fit_dist_.84, fit_temp_.84, fit_cold_.84, fit_winter_.84, fit_precip_.84)
do.call(rbind, lapply(1:length(models), function(i) glance(models[[i]]) %>% mutate(model = model_names[i]))) %>%
  arrange(AIC)

#-----------------------2014 vs. 2018 cline comparion--------------------#
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
summary(fit_lat_2014_2018_.84)
tidy(fit_lat_2014_2018_.84)
glance(fit_lat_2014_2018_.84)

coef_lat_2014_2018_.84 = tidy(fit_lat_2014_2018_.84) %>%
  dplyr::select(., term, estimate) %>%
  spread(., term, estimate) %>%
  mutate(width_2018 = abs(4/b),
         width_2014 = abs(4/(b + b_2014)))
coef_lat_2014_2018_.84

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
coef_lat_zone_mu_b_.84_2018


# feral hives
fit_lat_zone_mu_and_b_.84_feral <- nls_multstart(alpha ~ logistic4(x = abs_lat, 
                                                                   b = b + b_feral*enjambre, 
                                                                   mu = mu + mu_feral*enjambre,
                                                                   K = 0.84),
                                                 start_lower = list(b = -5, mu = min(d_A$abs_lat), b_feral = -5, mu_feral = -10),
                                                 start_upper = list(b = 5, mu = max(d_A$abs_lat), b_feral = 5, mu_feral = 10),
                                                 supp_errors = 'Y',
                                                 iter = 250,
                                                 convergence_count = 100,
                                                 data = d_A)
summary(fit_lat_zone_mu_and_b_.84_feral)

# or just simple logistic
fit_logistic_lat_feral <- nls_multstart(alpha ~ logistic(x = a + b_lat*abs_lat + b_feral*enjambre),
                                        start_lower = list(a = -30, b_lat = -5, b_feral = -5),
                                        start_upper = list(a = 30, b_lat = 5, b_feral = 5),
                                        supp_errors = 'Y',
                                        iter = 250,
                                        convergence_count = 100,
                                        data = d_A)
summary(fit_logistic_lat_feral)

# ok so the slopes go down with higher recomb, but is this a big enough effect ot be interesting? not sure. let's plot:
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

