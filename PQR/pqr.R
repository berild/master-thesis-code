library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(SemiPar)

# loading general functions
source("./PQR/general_functions.R")
# loading amis with inla functions
source("./PQR/amis_w_inla.R")
# loading is with inla functions
source("./PQR/is_w_inla.R")
# loading mcmc with inla functions
source("./PQR/mcmc_w_inla.R")

# generating gamma pqr dataset
ggamma__model <- function(params, n = 200){
  x = runif(n)
  mu = params[[1]] + params[[2]]*x
  k = exp(params[[3]] + params[[4]]*x)
  theta = exp(mu)/k
  scale = theta
  y = rgamma(n, shape = k, scale = scale)
  return(list(
    data = data.frame(x = x, y = y),
    params = params
  ))
}

# function fitting conditional models with INLA with gamma likelihood
fit.inla.ggamma <- function(data,eta){
  res = inla(y ~ 1 + x, 
             data = data,
             scale = exp(eta[1] + eta[2]*data$x), 
             family = "gamma",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(a = res$marginals.fixed[[1]],
                           b = res$marginals.fixed[[2]])))
}

# function fitting conditional models with INLA with gaussian likelihood
fit.inla.gaussian <- function(data,eta){
  res = inla(y ~ 1 + x,
             data = data, 
             family = "gaussian",
             scale = exp(eta[1] + eta[2]*data$x),
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(a = res$marginals.fixed[[1]],
                           b = res$marginals.fixed[[2]])))
}

# function fitting conditional RW2 models with INLA 
fit.inla.rw2 <- function(data,eta){
  formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
  res = inla(formula, 
             data = data,
             control.predictor = list(compute = T),
             scale = exp(eta[1] + eta[2]*data$range), 
             family = "gaussian",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           prec.range = res$marginals.hyperpar[[1]])))
}

# generating gaussian pqr dataset
gaussian_model <- function(params, n = 200){
  x = runif(n)
  mu = params[[1]] + params[[2]]*x
  tau = exp(params[[3]] + params[[4]]*x)
  y = rnorm(n,mean = mu, sd = 1/sqrt(tau))
  return(list(
    data = data.frame(x = x, y = y),
    params = params
  ))
}

# vague prior for the parameters
prior.param <- function(x, log = TRUE) {
  sum(dnorm(x, 0, 100, log = log))
}

# setting parameters for a vague proposal distribution AMIS and IS
init = list(mu = c(0,0),cov = 10*diag(2))

# proposal distribution amis and is
dq.param <- function(y, x, sigma = init$cov, log =TRUE) {
  dmvt(y,sigma = sigma, df=3, delta = x, type = "shifted",log=log)
  #dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.param <- function(x, sigma = init$cov) {
  as.vector(rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted"))
  #as.vector(rmvnorm(1, mean = x, sigma = sigma))
}


# M1 Gaussian
params = list(a = 1,b = -0.1,c = -1.5,d = -2)
# Simulating data
set.seed(2)
mod = gaussian_model(params,n=500)
#fit model using AMIS w/ INLA
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.gaussian,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gaussian",kde = T)
amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
amis_w_inla_mod$mod = mod
save(amis_w_inla_mod, file = "./sims/pqr-m1-gaussian-amis-w-inla.Rdata")
# fit model using IS w/ INLA
is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
                          dq.param, rq.param, fit.inla.gaussian,
                          N_0 = 2000, N = 10000, pqr = "gaussian", kde = T)
is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
is_w_inla_mod$mod = mod
save(is_w_inla_mod, file = "./sims/pqr-m1-gaussian-is-w-inla.Rdata")


# M2 Gaussian
params = list(a=-2,b=-5,c=-2,d=3.5)
# Simulating data
set.seed(2)
mod = gaussian_model(params,n=500)
#fit model using AMIS w/ INLA
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.gaussian,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gaussian",kde = T)
amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
amis_w_inla_mod$mod = mod
save(amis_w_inla_mod, file = "./sims/pqr-m2-gaussian-amis-w-inla.Rdata")
# fit model using IS w/ INLA
is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
                          dq.param, rq.param, fit.inla.gaussian,
                          N_0 = 2000, N = 10000, pqr = "gaussian", kde = T)
is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
is_w_inla_mod$mod = mod
save(is_w_inla_mod, file = "./sims/pqr-m2-gaussian-is-w-inla.Rdata")

# M1 Gamma
params = list(a=3,b=-1,c=-1,d=6)
# Simulating data
set.seed(2)
mod = ggamma__model(params,n=500)
# fit model using AMIS w/ INLA
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.ggamma,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gamma",kde = T)
amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
amis_w_inla_mod$mod = mod
save(amis_w_inla_mod, file = "./sims/pqr-m1-gamma-amis-w-inla.Rdata")
# fit model using IS w/ INLA
is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
                          dq.param, rq.param, fit.inla.ggamma,
                          N_0 = 2000, N = 10000, pqr = "gamma", kde = T)
is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
is_w_inla_mod$mod = mod
save(is_w_inla_mod, file = "./sims/pqr-m1-gamma-is-w-inla.Rdata")


# M2 Gamma
params = list(a=-3,b=2,c=4,d=-1)
# Simulating data
set.seed(2)
mod = ggamma__model(params,n=500)
# fit model using AMIS w/ INLA
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.ggamma,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gamma",kde = T)
amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
amis_w_inla_mod$mod = mod
save(amis_w_inla_mod, file = "./sims/pqr-m2-gamma-amis-w-inla.Rdata")
# fit model using IS w/ INLA
is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
                          dq.param, rq.param, fit.inla.ggamma,
                          N_0 = 2000, N = 10000, pqr = "gamma", kde = T)
is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
is_w_inla_mod$mod = mod
save(is_w_inla_mod, file = "./sims/pqr-m2-gamma-is-w-inla.Rdata")

# Lidar data
data(lidar)
#scaling data for easier convergence with current proposal
r_max = max(lidar$range)
r_min = min(lidar$range)
lidar$range = lidar$range/r_max

# changing intital proposal distribution amis and is
init = list(mu = c(10,-10),cov = 3*diag(2))

# running amis with inla on lidar data
amis_w_inla_mod = amis.w.inla(data = lidar, init = init, prior.param,
                              dq.param, rq.param, fit.inla.rw2,
                              N_t = seq(25,50,1)*10, N_0 = 250,kde = T)
amis_w_inla_mod$mod = lidar
amis_w_inla_mod$scale = c(r_max,r_min)
save(amis_w_inla_mod, file = "./sims/pqr-gaussian-lidar-amis-w-inla.Rdata")

# running is with inla on lidar data
is_w_inla_mod = is.w.inla(data = lidar, init = init, prior.param,
                          dq.param, rq.param, fit.inla.rw2,
                          N_0 = 800, N = 10000)
is_w_inla_mod$mod = lidar
is_w_inla_mod$scale = c(r_max,r_min)
save(is_w_inla_mod, file = "./sims/pqr-gaussian-lidar-is-w-inla.Rdata")

# initalizing proposal distribution for mcmc with inla 
init_mcmc = list(mu = c(0,0),cov = diag(2))
dq.param.mcmc <- function(y, x, sigma = init$cov, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.param.mcmc  <- function(x, sigma = init$cov) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
}
# running mcmc with inla on lidar data
mcmc_w_inla_mod <- mcmc.w.inla(data = lidar, init = init_mcmc,
                               prior.param, dq.param.mcmc, rq.param.mcmc, fit.inla.rw2,
                               n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_w_inla_mod, file = "./sims/pqr-gaussian-lidar-mcmc-w-inla.Rdata")
eta_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1,nrow(mcmc_w_inla_mod$eta))/nrow(mcmc_w_inla_mod$eta), n = 100, lims = c(-1,1,-1,1))
mcmc_w_inla_mod$eta_kern_joint = data.frame(expand.grid(x=eta_kern_mcmc$x, y=eta_kern_mcmc$y), z=as.vector(eta_kern_mcmc$z))
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                        kernel = "gaussian")[c(1,2)])
})
mcmc_w_inla_mod$mod = lidar
mcmc_w_inla_mod$scale = c(r_max,r_min)
save(mcmc_w_inla_mod, file = "./sims/pqr-gaussian-lidar-mcmc-w-inla.Rdata")

# # ImmunogG data
data("ImmunogG")
# changing intial proposal distribution amis and is
init = list(mu = c(2,0),cov = diag(c(1,0.5)))

#running amis with inla on ImmunogG data
amis_w_inla_mod = amis.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG),
                              init = init, prior.param, dq.param, rq.param, fit.inla.ggamma,
                              N_t = seq(25,50,1)*10, N_0 = 250, pqr = "gamma",kde = T)
amis_w_inla_mod$mod = list(x=ImmunogG$Age,y = ImmunogG$IgG)
save(amis_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-amis-w-inla.Rdata")

# running is with inla on ImmunogG data
is_w_inla_mod = is.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG),
                          init = init, prior.param, dq.param, rq.param, fit.inla.ggamma, N_0 = 800,N = 10000)
is_w_inla_mod$mod = list(x=ImmunogG$Age,y = ImmunogG$IgG)
save(is_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-is-w-inla.Rdata")

# running mcmc with inla on ImmunogG data
init_mcmc = list(mu = c(0,0),cov = diag(c(0.3,0.15)))
mcmc_w_inla_mod <- mcmc.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG), init = init_mcmc,
                               prior.param, dq.param.mcmc, rq.param.mcmc, fit.inla.ggamma,
                               n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-mcmc-w-inla.Rdata")
eta_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1,nrow(mcmc_w_inla_mod$eta))/nrow(mcmc_w_inla_mod$eta), n = 100, lims = c(-1,1,-1,1))
mcmc_w_inla_mod$eta_kern_joint = data.frame(expand.grid(x=eta_kern_mcmc$x, y=eta_kern_mcmc$y), z=as.vector(eta_kern_mcmc$z))
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                        kernel = "gaussian")[c(1,2)])
})
mcmc_w_inla_mod$mod = list(x=ImmunogG$Age,y = ImmunogG$IgG)
save(mcmc_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-mcmc-w-inla.Rdata")