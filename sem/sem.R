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

# source("./sem/amis_w_inla.R")
# amis_w_inla_mod <- amis.w.inla(data = turnout, init = init, prior.rho.lambda,
#                                dq.rho.lambda, rq.rho.lambda, fit.inla,
#                                N_t = seq(25,50,1)*10, N_0 = 250)
# save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")
# eta_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
# amis_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_amis$x, y=eta_kern_amis$y), z=as.vector(eta_kern_amis$z))
# save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")
# 
# source("./sem/is_w_inla.R")
# is_w_inla_mod <- is.w.inla(data = turnout, init = init, prior.rho.lambda,
#                            dq.rho.lambda, rq.rho.lambda,fit.inla, N_0 = 800, N = 10000)
# save(is_w_inla_mod, file = "./sims/sem-is-w-inla.Rdata")
# eta_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
# is_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_is$x, y=eta_kern_is$y), z=as.vector(eta_kern_is$z))
# save(is_w_inla_mod, file = "./sims/sem-is-w-inla.Rdata")

init$cov = 0.1*init$cov
source("./sem/sem_mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = turnout, init = init,
                               prior.rho.lambda, dq.rho.lambda, rq.rho.lambda, fit.inla,
                               n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_w_inla_mod, file = "./sims/sem-mcmc-w-inla.Rdata")
eta_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = mcmc_w_inla_mod$weight/(sum(mcmc_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
mcmc_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_mcmc$x, y=eta_kern_mcmc$y), z=as.vector(eta_kern_mcmc$z))
save(mcmc_w_inla_mod, file = "./sims/sem-mcmc-w-inla.Rdata")


# library(ggplot2)
# 
# p1 <- ggplot() +
#   geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,color="AMIS with INLA")) +
#   geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y,color="IS with INLA")) +
#   theme_bw() +
#   theme(legend.position="bottom")
# p1
# 
# 
# p2 <- ggplot() +
#   geom_line(data = amis_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="AMIS with INLA")) +
#   geom_line(data = is_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="IS with INLA")) +
#   theme_bw() +
#   theme(legend.position="bottom")
# p2
# 
# p3 <- ggplot() +
#   geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,color="AMIS with INLA")) +
#   geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,color="IS with INLA")) +
#   theme_bw() +
#   theme(legend.position="bottom")
# p3
# 
# p4 <- ggplot() +
#   geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
#   geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
#   labs(color = "",x=expression(rho),y=expression(lambda),linetype="") +
#   theme_bw() +
#   theme(legend.position="bottom")
# p4
# 
# spplot(turnout["TURNOUT01"],  col.regions= hcl.colors(16,palette =  hcl.pals("sequential")[21], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
# 
# spplot(turnout["GDPCAP"], col.regions= hcl.colors(16,palette =  hcl.pals("sequential")[21], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
