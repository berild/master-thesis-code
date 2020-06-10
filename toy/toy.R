library(parallel)
library(mvtnorm)
library(MASS)
library(coda)
library(INLA)

sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 1 + 1*x1 -1*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}

# sampling dataset
set.seed(1)
df = sample.linreg()

inla_mod = inla(y~x, data=df) 
save(inla_mod, file = "./sims/toy-inla.Rdata")


fit.inla <- function(data, beta){
  data$oset = data$x%*%t(beta)
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = res$marginals.fixed[[1]], 
                           tau = res$marginals.hyperpar[[1]]),
              stats = list(intercept = as.numeric(res$summary.fixed[1]),
                           tau  = as.numeric(res$summary.hyperpar[1]))))
}


prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(x, mean = 0, sd= sigma, log = log))
}

source("./toy/general_functions.R")

init = list(mu = c(0,0),cov = diag(5,2,2))

rq.beta <- function(x, sigma = diag(5,2,2)) {
  rmvnorm(1,mean=x,sigma = sigma)
}

dq.beta <- function(y, x, sigma = diag(5,2,2), log =TRUE) {
  dmvnorm(y,mean = x, sigma = sigma,log = log)
}

source("./toy/amis_w_inla.R")
amis_w_inla_mod <- amis.w.inla(data = df, init = init, prior.beta,
                               dq.beta, rq.beta, fit.inla,
                               N_t = seq(25,50)*10, N_0 = 250)
save(amis_w_inla_mod, file = "./sims/toy-amis-w-inla.Rdata")

source("./toy/is_w_inla.R")
is_w_inla_mod <- is.w.inla(data = df, init = init, prior.beta,
                           dq.beta, rq.beta, fit.inla, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./sims/toy-is-w-inla.Rdata")

rq.beta <- function(x, sigma = .75) {
  rnorm(length(x), mean = x, sd = sigma)
}

dq.beta <- function(y, x, sigma = .75, log =TRUE) {
  sum(dnorm(x, mean = y, sd = sigma, log = log))
}

init = list(mu = c(0,0),cov = 0.5*diag(2))
source("./toy/mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = init, prior.beta,
                               dq.beta, rq.beta, fit.inla, n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_w_inla_mod, file = "./sims/toy-mcmc-w-inla.Rdata")
