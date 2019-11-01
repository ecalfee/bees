library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)

# load in example data set
data("Chlorella_TRC")

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}

# subset dataset
d_1 <- subset(Chlorella_TRC, curve_id == 1)


# run nls_multstart with shotgun approach
fit <- nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                     data = d_1,
                     iter = 250,
                     start_lower = c(lnc=-10, E=0.1, Eh=0.5, Th=285),
                     start_upper = c(lnc=10, E=2, Eh=5, Th=330),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
fit_bee <- nls_multstart(A ~ logistic3(x = lat, b = b, mu = mu),
                         start_lower = list(b = 0, mu = -40),
                         start_upper = list(b = 1, mu = -25),
                         supp_errors = 'Y',
                         iter = 250,
                         convergence_count = 100,
                         data = d0)


fit
fit_bee
glance(fit_bee)
tidy(fit_bee)
# get confidence intervals using nlstools
CI_bee <- confint2(fit_bee) %>%
  data.frame() %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..)

# bind params and confidence intervals
params_bee <- bind_cols(tidy(fit_bee), CI_bee)
select(params_bee, -c(statistic, p.value))
preds_bee = augment(fit_bee) # predict
# plot fit
ggplot() +
  geom_point(data = d_1, aes(K, ln.rate)) +
  geom_line(data = preds, aes(K, .fitted))
ggplot() +
  geom_point(data = d0, aes(lat, A)) +
  geom_line(data = preds_bee, aes(lat, .fitted))
with(preds_bee, cor(A, .fitted))

# fit multiple curves:
# fit over each set of groupings
fits <- Chlorella_TRC %>%
  group_by(., flux, growth.temp, process, curve_id) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                                data = .x,
                                                iter = 1000,
                                                start_lower = c(lnc=-1000, E=0.1, Eh=0.5, Th=285),
                                                start_upper = c(lnc=1000, E=2, Eh=10, Th=330),
                                                supp_errors = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))

fits_bee <- data.frame(index = 1:10) %>%
  mutate(fit = lapply(index, 
                      function(i) nls_multstart(A ~ logistic3(x = lat, 
                                                     b = b, mu = mu),
                      start_lower = list(b = 0, mu = -40),
                      start_upper = list(b = 1, mu = -25),
                      supp_errors = 'Y',
                      iter = 250,
                      convergence_count = 100,
                      data = data.frame(A = unname(t(A[i, meta.AR.order.by.lat$population])), 
                                        lat = meta.AR.order.by.lat$lat, 
                                        stringsAsFactors = F))))

fits_bee <- lapply(1:10, 
                   function(i) nls_multstart(A ~ logistic3(x = lat, 
                                                           b = b, mu = mu),
                                             start_lower = list(b = 0, mu = -40),
                                             start_upper = list(b = 1, mu = -25),
                                             supp_errors = 'Y',
                                             iter = 250,
                                             convergence_count = 100,
                                             data = data.frame(A = unname(t(A[i, meta.AR.order.by.lat$population])), 
                                                               lat = meta.AR.order.by.lat$lat, 
                                                               stringsAsFactors = F)))


# look at output object
select(fits, curve_id, data, fit)
select(fits_bee, index, fit)

# look at a single fit
summary(fits$fit[[1]])
summary(fits_bee$fit[[1]])

# get summary
info <- fits %>%
  mutate(summary = map(fit, glance)) %>%
  unnest(summary)
info_bee <- fits_bee %>%
  mutate(summary = map(fit, glance)) %>%
  unnest(summary)

# get params
params <- fits %>%
  mutate(., p = map(fit, tidy)) %>%
  unnest(p)
params_bee <- fits_bee %>%
  mutate(., p = map(fit, tidy)) %>%
  unnest(p)

# get confidence intervals
CI <- fits %>%
  mutate(., cis = map(fit, confint2),
         cis = map(cis, data.frame)) %>%
  unnest(cis) %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(., curve_id) %>%
  mutate(., term = c('lnc', 'E', 'Eh', 'Th')) %>%
  ungroup() %>%
  select(., -data, -fit)

CI_bee <- fits_bee %>%
  mutate(., cis = map(fit, confint2),
         cis = map(cis, data.frame)) %>%
  unnest(cis) %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(., index) %>%
  mutate(., term = c('mu', 'b')) %>%
  ungroup() %>%
  select(., -fit)

# merge parameters and CI estimates
params2 <- merge(params, CI, by = intersect(names(params), names(CI)))
params2_bee <- merge(params_bee, CI_bee, by = intersect(names(params_bee), names(CI_bee)))


# get predictions
preds <- fits %>%
  mutate(., p = map(fit, augment)) %>%
  unnest(p)
preds_bee <- fits_bee %>%
  mutate(., p = map(fit, augment)) %>%
  unnest(p)

# did models converge?
select(info, curve_id, logLik, AIC, BIC, deviance, df.residual)
select(info_bee, index, logLik, AIC, BIC, deviance, df.residual)

# plot:
# new data frame of predictions
new_preds <- Chlorella_TRC %>%
  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Chlorella_TRC, curve_id) %>%
  summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create new predictions
preds2 <- fits %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'curve_id') %>%
  group_by(., curve_id) %>%
  filter(., K > unique(min_K) & K < unique(max_K)) %>%
  rename(., ln.rate = .fitted) %>%
  ungroup()

preds2_bee <- fits_bee %>%
  mutate(., p = map(fit, augment, newdata = data.frame(lat = seq(from = min(meta.AR.order.by.lat$lat),
                                                to = max(meta.AR.order.by.lat$lat),
                                                length.out = 100)))) %>%
  unnest(p) %>%
  group_by(., index) %>%
  rename(., A = .fitted) %>%
  ungroup()

# plot
ggplot() +
  #geom_point(aes(A, lat, col = index), size = 2, head(A)) +
  geom_line(aes(A, lat, col = factor(index)), alpha = 0.5, preds2_bee) +
    theme_bw(base_size = 12) +
  ylab('A') +
  xlab('lat')

ggplot() +
  geom_point(aes(K - 273.15, ln.rate, col = flux), size = 2, Chlorella_TRC) +
  geom_line(aes(K - 273.15, ln.rate, col = flux, group = curve_id), alpha = 0.5, preds2) +
  facet_wrap(~ growth.temp + process, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4', 'black')) +
  theme_bw(base_size = 12) +
  ylab('log Metabolic rate') +
  xlab('Assay temperature (ºC)') +
  theme(legend.position = c(0.9, 0.15))

