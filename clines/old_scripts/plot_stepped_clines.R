# code to fit various clines 
# load data first from plot_clines.R
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


