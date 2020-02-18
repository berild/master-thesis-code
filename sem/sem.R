library(spdep)
library(spData)
library(parallel)
library(mvtnorm)
library(MASS)
library(coda)


library(rgdal)
library(sp)
library(RColorBrewer)
library(spdep)

# LOad data
load("./data/load_data.RData")

# Models with INLABMA
library(INLA)
library(INLABMA)


# Set hthis for parallel computing
inla.setOption(num.threads = 2)
# No Gaussian Error
zero.variance <- list(prec = list(initial = 15, fixed = TRUE))
# Number of grid points in each dimension
n.xy <-  c(40, 20) #rho, lambda)

turnout$idx <- 1:nrow(turnout)
form <- TURNOUT01 ~ 1 + log(GDPCAP)

init = list(mu = c(0,0),cov = diag(2))

source("./sem/general_functions.R")

source("./sem/amis_w_inla.R")
amis_w_inla_mod <- amis.w.inla(data = turnout, init = init, prior.rho.lambda, 
                               dq.rho.lambda, rq.rho.lambda, fit.inla, 
                               N_t = seq(25,50,1)*10, N_0 = 250)
save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")

#source("./sem/sem_is_w_inla.R")
#is_w_inla_mod <- is.w.inla(data = df, init = init, prior.rho, 
#                           dq.rho, rq.rho,fit.inla, N_0 = 800, N = 10000)
#save(is_w_inla_mod, file = "./sem/sims/sem-is-w-inla.Rdata")
#
#source("./sem/sem_mcmc_w_inla.R")
#mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = init$mu,
#                               prior.rho, dq.rho, rq.rho, fit.inla,
#                               n.samples = 100500, n.burnin = 500, n.thin = 10)
#save(mcmc_w_inla_mod, file = "./sem/sims/sem-mcmc-w-inla.Rdata")


