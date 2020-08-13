library(spdep)
library(spData)
library(spatialreg)
library(parallel)
library(mvtnorm)
library(MASS)
library(coda)
library(rgdal)
library(sp)
library(RColorBrewer)

# LOad data
load("./data/load_data.RData")
library(INLA)


# Set hthis for parallel computing in INLA
inla.setOption(num.threads = 2)
# No Gaussian Error
zero.variance <- list(prec = list(initial = 15, fixed = TRUE))
# Number of grid points in each dimension
n.xy <-  c(40, 20) #rho, lambda

# indexing regions in map
turnout$idx <- 1:nrow(turnout)

# prior of autoregressive parameters
prior.rho.lambda <- function(x, log = TRUE) {
  sum(dunif(x, -1, 1, log = log))
}

# loading general functions
source("./sem/general_functions.R")

# function fitting conditional models with INLA
fit.inla <- function(data,eta){
  res <- sac.inla(TURNOUT01 ~ 1 + log(GDPCAP), 
                  d = as.data.frame(turnout), W.rho = W, W.lambda = W,
                  fhyper = list(prec = list(param = c(0.01, 0.01))),
                  rho = eta[1],
                  lambda = eta[2],
                  family = "gaussian", impacts = FALSE,
                  control.fixed = list(prec.intercept = 0.001),
                  control.family = list(hyper = zero.variance),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, cpo = TRUE),
                  control.inla = list(print.joint.hyper = TRUE, #), 
                                      strategy = "laplace", tolerance = 1e-10, h = 0.001),
                  # Start close to ML estimate
                  control.mode = list(theta = log(0.2), restart = TRUE),
                  improve = FALSE,
                  verbose = FALSE
  )
  logdet <- res$logdet
  res <- inla.rerun(res)
  res$mlik <- res$mlik + 0.5 * logdet
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           GDPCAP = res$marginals.fixed[[2]],
                           tau = res$marginals.hyperpar[[1]]),
              logdet = logdet))
}

# intial parameters of the proposal distribution in AMIS and IS
init = list(mu = c(0,0),cov = diag(2))

# proposal distribution of autoregressive parameters in AMIS and IS
dq.rho.lambda <- function(y, x, sigma = init$cov, log =TRUE) {
  dmvt(y,sigma = sigma, df=3, delta = x, type = "shifted",log=log)
}
rq.rho.lambda <- function(x, sigma = init$cov) {
  as.vector(rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted"))
}

# loading AMIS with INLA functions
source("./sem/amis_w_inla.R")
amis_w_inla_mod <- amis.w.inla(data = turnout, init = init, prior.rho.lambda,
                               dq.rho.lambda, rq.rho.lambda, fit.inla,
                               N_t = seq(25,50,1)*10, N_0 = 250)
save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")
# approximating densities with the weighted set of samples
eta_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
amis_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_amis$x, y=eta_kern_amis$y), z=as.vector(eta_kern_amis$z))
save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")

# loading IS with INLA functions
source("./sem/is_w_inla.R")
is_w_inla_mod <- is.w.inla(data = turnout, init = init, prior.rho.lambda,
                           dq.rho.lambda, rq.rho.lambda,fit.inla, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./sims/sem-is-w-inla.Rdata")
# approximating densities with the weighted set of samples
eta_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
is_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_is$x, y=eta_kern_is$y), z=as.vector(eta_kern_is$z))
save(is_w_inla_mod, file = "./sims/sem-is-w-inla.Rdata")


# initial state and proposal distribution for MCMC with INLA
init = list(mu = c(0,0),cov = 0.5*diag(2))

# proposal distribution autoregressive parameters in MCMC with INLA
dq.rho.lambda <- function(y, x, sigma = init$cov, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
}
rq.rho.lambda <- function(x, sigma = init$cov) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
}

# loading mcmc with inla functions
source("./sem/mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = turnout, init = init,
                               prior.rho.lambda, dq.rho.lambda, rq.rho.lambda, fit.inla,
                               n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_w_inla_mod, file = "./sims/sem-mcmc-w-inla.Rdata")
# approximating densities with samples
eta_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1,nrow(mcmc_w_inla_mod$eta))/nrow(mcmc_w_inla_mod$eta), n = 100, lims = c(-1,1,-1,1))
mcmc_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_mcmc$x, y=eta_kern_mcmc$y), z=as.vector(eta_kern_mcmc$z))
save(mcmc_w_inla_mod, file = "./sims/sem-mcmc-w-inla.Rdata")

# running mcmc simulation on the dataset using spatialreg package
mcmc_mod <- spBreg_sac(TURNOUT01 ~ log(GDPCAP), data = as.data.frame(turnout),
                       listw = it.lw, control = list(ndraw = 110000L, nomit = 10000L, thin = 10L,
                                                     prior = list(nu = 0.01, d0 = 0.01, a1 = 1, a2 = 1,
                                                                  rho = sac1$rho, lambda = sac1$lambda, Tbeta = diag(2) * 1000)))

mcmc_mod = list(samples = as.data.frame(mcmc_mod),
                times = attributes(mcmc_mod)$timings)
mcmc_mod$samples$tau = 1/mcmc_mod$samples$sige
eta_kern_mc = kde2d.weighted(x = mcmc_mod$samples$rho, y = mcmc_mod$samples$lambda, w = rep(1,nrow(mcmc_mod$samples))/nrow(mcmc_mod$samples), n = 100, lims = c(-1,1,-1,1))
mcmc_mod$eta_kern = data.frame(expand.grid(x=eta_kern_mc$x, y=eta_kern_mc$y), z=as.vector(eta_kern_mc$z))
mcmc_mod$margs = lapply(seq(ncol(mcmc_mod$samples)), function(x){
  as.data.frame(density(x = mcmc_mod$samples[,x],
                        weights = rep(1/nrow(mcmc_mod$samples),nrow(mcmc_mod$samples)),
                        kernel = "gaussian")[c(1,2)])
})
names(mcmc_mod$margs) = names(mcmc_mod$samples)
save(mcmc_mod,file = "./sims/sem/sem-mcmc.Rdata")


# BMA with INLA from original paper
source("./sem/bma_w_inla.R")
eta.grid = as.matrix(expand.grid(seq(-1,1,length.out=100),seq(-1,1,length.out=100)))
bma_w_inla_mod <- bma.w.inla(data=turnout,eta.grid, prior.rho.lambda,fit.inla)
eta_kern = kde2d.weighted(x = bma_w_inla_mod$eta[,1], y = bma_w_inla_mod$eta[,2], w = bma_w_inla_mod$weight, n = 100, lims = c(-1,1,-1,1))
bma_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern$x, y=eta_kern$y), z=as.vector(eta_kern$z))
save(mcmc_w_inla_mod, file = "./sims/sem-bma-w-inla.Rdata")

