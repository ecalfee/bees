# define cline shapes:

# define logistic curve
logistic3 <- function(x, mu, b){
  # mu is the cline center and b = 4/w where w is the cline width (inverse maximum slope)
  1/(1 + exp(-b*(x - mu)))
}

logistic4 <- function(x, mu, b, K){ # K is the maximum value
  K/(1 + exp(-b*(x - mu))) # b = 4/w. If K = 1, this is the same cline as 'stepped cline center' below
}
# logistic4() adds asymptote as a free parameter
# why K?
# because bees in brazil have ~ 84% A ancestry, the true
# asymptote is 84% A not 100% A for the cline, so
# I can optionally use this to rescale the inferred logistic curves
# note: interpretation of mu and b with the scaling parameter K is a little less direct
