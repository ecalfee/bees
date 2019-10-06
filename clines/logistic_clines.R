# clines
library(dplyr)

# logistic
logistic <- function(x){
  1/(1 + exp(-x))
}
logistic3 <- function(x, mu, b){
  1/(1 + exp(-b*(x - mu)))
}
logistic2 <- function(x, Asym, xmid, scal){
  Asym/(1 + exp(-(x - xmid)/scal))
}
# simulate some data:
set.seed(101)
x = runif(100, -10, 10)
y = logistic(x + rnorm(100, 0, .1))
plot(x, y)
curve(logistic(x), -10, 10, add = T, col = 'blue')

# fit a logistic:
fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), 
                start=list(a=1,b=.5,c=25))
nls(y ~ logistic(b*(x-mu)),
    start = list(b = 1, mu = 0))
nls(y ~ logistic3(x, mu, alpha, scale),
    start = list(alpha = 1, mu = 0.1, scale = 1),
    trace = T)
nls(y ~ logistic3(x, mu, alpha, scale),
    start = ,
    trace = T)
nls(y ~ logistic2(x, Asym, xmid, scal),
    start = getInitial(y ~ SSlogis(x, Asym, 
                                   xmid, scal), 
                       d = list(x = x, y = y)),
    trace = T)

# now try to get initial self-starting values
getInitial(y ~ SSlogis(x, Asym, xmid, scal), 
           d = list(x = x, y = y))

# load my hybrid zone data
# load population metadata
load("../local_ancestry/results/pops_by_lat.RData") # contains objects pops_by_lat meta.pop and meta.AR.order.by.lat 

# load mean ancestry by population
# .RData file contains objects 'sites' which has columns chr and pos and 'DATA_NAME' which has ancestry where sites are rows and pops are columns
load(paste0("../local_ancestry/results/A.RData")) 

# subset data to pops of interest, only SA
A_snp_all <- A[ , meta.AR.order.by.lat$population]

snp1 <- A_snp_all[1, ] %>%
  tidyr::gather(., "population", "A") %>%
  left_join(meta.AR.order.by.lat, ., by = "population")

nls(A ~ logistic2(lat, Asym, xmid, scal),
    start = getInitial(A ~ SSlogis(lat, Asym, 
                                   xmid, scal), 
                       d = snp1),
    data = snp1,
    trace = T)
start1 <- getInitial(A ~ SSlogis(lat, Asym, 
                                 xmid, scal), 
                     d = snp1)
start2 <- c(b = unname(1/start1["scal"]), 
            mu = unname(start1["xmid"]))
nls(A ~ logistic2(lat, Asym, xmid, scal),
    start = start1,
    data = snp1,
    trace = T)
nls(A ~ logistic3(x = lat, b = b, mu = mu),
    start = list(b = 1/unname(start1["scal"]),
                 mu = unname(start1["xmid"])),
    data = snp1,
    trace = T)
nls(A ~ logistic3(x = lat, b = 1/scal, mu = xmid),
    start = start1[c("xmid", "scal")],
    data = snp1,
    trace = T)
nls(A ~ logistic2(x = lat, Asym = 1, 
                  xmid = xmid, scal = scal),
    start = start1[c("xmid", "scal")],
    data = snp1,
    trace = T)
nls(A ~ logistic3(x = lat, b, mu),
    start = start2,
    data = snp1,
    trace = T)
logistic3 <- function(x, mu, b){
  1/(1 + exp(-b*(x - mu)))
}
logistic2 <- function(x, Asym, xmid, scal){
  Asym/(1 + exp(-(x - xmid)/scal))
}

nls(y ~ logistic2(x, Asym, xmid, scal),
    start = getInitial(y ~ SSlogis(x, Asym, 
                                   xmid, scal), 
                       d = list(x = x, y = y)),
    trace = T)

fit_cline0 <- function(snp){
  start0 <- getInitial(A ~ SSlogis(lat, Asym, 
                                   xmid, scal), 
                       d = snp)
  fit0 <- nls(A ~ logistic3(x = lat, b = b, mu = mu),
      start = list(b = 1/unname(start1["scal"]),
                   mu = unname(start1["xmid"])),
      data = snp,
      trace = F)
  sum0 <- summary(fit0)
  fit <- c(converged = sum0$convInfo$isConv,
           mu = sum0$coefficients["mu", "Estimate"],
           b = sum0$coefficients["b", "Estimate"],
           residual_error = sum0$sigma)
  return(fit)
}
fit_cline <- function(snp){ # wrapper to return an error message and empty vector if fit fails
  tryCatch(fit_cline0(snp),
           error=function(e) {
             print(e)
             return(c(converged = NA,
                      mu = NA,
                      b = NA,
                      residual_error = NA))
           }) 
}
A_snp_subset <- dplyr::sample_n(A_snp_all, 100)
head(A_snp_subset)
A_snp_subset$AR01[1] <- 0
A_snp_subset$AR27[1] <- 1
some_clines <- lapply(1:10, function(i){
  A_snp_subset[i, ] %>%
    tidyr::gather(., "population", "A") %>%
    left_join(meta.AR.order.by.lat, ., by = "population") %>%
    fit_cline(.)
})
system.time(fit_cline(snp = snp1))

  

