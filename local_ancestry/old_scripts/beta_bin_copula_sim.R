# simulating binomial ancestry draws that are correlated across individuals
# the beta part adds variance (drift) within individuals, the copula
# adds shared variance between, and the binomial does the sampling
# for 2 alleles
library(simstudy) # for R 3.5.2 but probs ok for my version R 3.5.1
# fake data, mean alpha and covariance structure K
K = rbind(c(.2, .1, .05), c(.1, .5, .2), c(.05, .2, 1))
alpha = c(.2, .5, .7)
pop_n = 1:3
#def <- simstudy::defData(varname = "pop_n", dist = "nonrandom", formula = pop_n)
#def <- simstudy::defData(varname = "alpha", dist = "nonrandom", formula = alpha)
def <- simstudy::defData(varname = "pop1", dist = "nonrandom", formula = alpha[1])
def <- simstudy::defData(def, varname = "pop2", dist = "nonrandom", formula = alpha[2])
def <- simstudy::defData(def, varname = "pop3", dist = "nonrandom", formula = alpha[3])
def <- simstudy::defData(def, varname = "A1", dist = "normal", var = .01, formula = "pop1")
def <- simstudy::defData(def, varname = "A2", dist = "normal", var = .01, formula = "pop2")
def <- simstudy::defData(def, varname = "A3", dist = "normal", var = .01, formula = "pop3")
def <- simstudy::defData(def, varname = "B1", dist = "binomial", var = .01, formula = "pop1")
def <- simstudy::defData(def, varname = "B2", dist = "binomial", var = .01, formula = "pop2")
def <- simstudy::defData(def, varname = "B3", dist = "binomial", var = .01, formula = "pop3")
sim <- genData(100, def)
apply(sim, 2, mean)
head(sim)
# simulate correlated data under MVN
# e.g. https://www.rdatagen.net/page/correlated/
sim2 <- genCorData(100, mu = alpha, sigma = sqrt(diag(K)), #sd 
                   corMatrix = cov2cor(K))
# show it worked:
sim2[, round(cor(cbind(V1, V2, V3)), 1)]
sim2[, round(sqrt(diag(var(cbind(V1, V2, V3)))), 1)]

# simulate correlated data under binomial:
sim3 <- genCorGen(100, 
          nvars = 3, 
          params1 = alpha, 
          params2 = 2, # 2 binary draws
          dist = "binomial", 
          corMatrix = cov2cor(K))

