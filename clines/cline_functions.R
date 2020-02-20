# define cline shapes:

# define logistic curve
logistic3 <- function(x, mu, b){
  # mu is the cline center and b = 4/w where w is the cline width (inverse maximum slope) 
  1/(1 + exp(-b*(x - mu)))
}
logistic4 <- function(x, mu, b, K){ # K is the maximum value
  K/(1 + exp(-b*(x - mu))) # b = 4/w. If K = 1, this is the same cline as 'stepped cline center' below
} # add asymptote free parameter
# why K?
# because bees in brazil have ~ 84% A ancestry, the true
# asymptote is 84% A not 100% A for the cline, so
# I can optionally use this to rescale the inferred logistic curves
# note: interpretation of mu and b with the scaling parameter K is a little less direct


# from Barton & Hewitt 1985 review, but cites Barton & Bengtsson in prep: barrier strength = width*(W_edge/W_center)^(2/r) where r is harmonic mean recomb. rate between selected & neutral loci
# the width of a cline is proportional to sqrt(var dispersal/s), where s is proportional to selection in a cline maintained by selection, or inverse time if a neutral cline spreading over time
# at the edges, theta = 1 if it's single locus sel. and you're observing the selected locus,
# otherwise the edges are 'shallower' than the middle and so theta < 1, and 
# you can interpret theta as the ratio between direct sel. on the marker itself and effective or indirect sel. due to LD at the center of the cline
# note Szymura and Barton, 1986 has incorrect eq. corrected by Szymura and Barton 1991 p. 11. Cite both.

szymura_barton_center <- function(x, y, w){
  # note w is cline width = inverse maximum slope
  # and y is the cline center
  p = (1 + tanh(2*(x-y)/w))/2
  return(p)
}

szymura_barton_edge <- function(x, y, theta, w, d){ # note: I made it x-y instead of just x to center. correct formula used in text of 1986 (but NOT fig 5) and in text 1991 paper (w/ correction stated for fig 5)
  p = exp(-4*(x - y)*sqrt(theta)/w)
  return(p)
}


# alternative formulation of the stepped cline from Gay 2008, eq. 1 (https://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2008.00491.x)
# note that y is the center, but W/4 is the maximum slope (so a rescaling of inverse of w above)
stepped_cline_center <- function(x, y, w){
  W = 4/w # note: W/4 is the maximum slope, so w is the 'cline width' = inverse maximum slpe
  p = 1/(1 + exp(-W*(x - y)))
  return(p)
}

stepped_cline_edge <- function(x, y, w, theta, d){ 
  W = 4/w # note: W/4 is the maximum slope, so w is the 'cline width'
  p = 1/(1 + exp(-W*d))*exp(W*theta*(x - y - d)/(1 + exp(W*d)))
  return(p)
}

stepped_cline_3parts <- function(x, y, w, thetaR, thetaL, dR, dL, rescale){
  p = ifelse(x > y + dR,
              1 - stepped_cline_edge(x, 
                                     y = y, 
                                     theta = thetaR, 
                                     w = -w, 
                                     d = dR),
              ifelse(x < y - dL,
                     stepped_cline_edge(x, 
                                        y = y, 
                                        theta = thetaL, 
                                        w = w, 
                                        d = -dL),
                     stepped_cline_center(x = x, 
                                          y = y, 
                                          w = w)))
  a = rescale*p
  return(a)
}
