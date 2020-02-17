library(INLA)
library(spdep)
library(spData)
library(spatialreg)
library(INLABMA)
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

# Fit models
library(parallel)
options(mc.cores = 4)

# Set hthis for parallel computing
inla.setOption(num.threads = 2)

# Number of grid points in each dimension
n.xy <-  c(40, 20) #rho, lambda)


form <- TURNOUT01 ~ 1 + log(GDPCAP)

init = list()

source("./sem/general_functions.R")

source("./sem/amis_w_inla.R")
amis_w_inla_mod <- amis.w.inla(data = turnout, init = init, prior.rho, 
                               dq.rho, rq.rho, fit.inla, 
                               N_t = seq(25,50,1)*10, N_0 = 250)
save(amis_w_inla_mod, file = "./sem/sims/sem-amis-w-inla.Rdata")

source("./sem/sem_is_w_inla.R")
is_w_inla_mod <- is.w.inla(data = df, init = init, prior.rho, 
                           dq.rho, rq.rho,fit.inla, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./sem/sims/sem-is-w-inla.Rdata")

source("./sem/sem_mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = init$mu,
                               prior.rho, dq.rho, rq.rho, fit.inla,
                               n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_w_inla_mod, file = "./sem/sims/sem-mcmc-w-inla.Rdata")


ggplot(italy_map,aes(long,lat,group=group))+
  geom_polygon(color = "white") + 
  theme(legend.position = "none")
