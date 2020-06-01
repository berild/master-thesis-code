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
#                                N_t = seq(25,50,1), N_0 = 25)
# save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")
# eta_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
# amis_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_amis$x, y=eta_kern_amis$y), z=as.vector(eta_kern_amis$z))
# save(amis_w_inla_mod, file = "./sims/sem-amis-w-inla.Rdata")

source("./sem/is_w_inla.R")
is_w_inla_mod <- is.w.inla(data = turnout, init = init, prior.rho.lambda,
                           dq.rho.lambda, rq.rho.lambda,fit.inla, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./sims/sem-is-w-inla.Rdata")
eta_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(-1,1,-1,1))
is_w_inla_mod$eta_kern = data.frame(expand.grid(x=eta_kern_is$x, y=eta_kern_is$y), z=as.vector(eta_kern_is$z))
save(is_w_inla_mod, file = "./sims/sem-is-w-inla.Rdata")

# p3 <- ggplot() + 
#   #geom_point(data = data.frame(x = 2, y = -2), aes(x = x, y = y), shape = 4,size = 3) + 
#   #geom_text(data = data.frame(x = 2, y = -2), aes(x = x, y = y), label="True Value", vjust=2) + 
#   #geom_line(data = data.frame(x=rep(1000,10),y = rep(1000,10),type = "INLA"), aes(x=x,y=y,linetype=type))+
#   geom_contour(data = eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) + 
#   #geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) + 
#   #geom_contour(data = mcmc_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) + 
#   labs(color = "",x=expression(rho),y=expression(lambda),linetype="") + 
#   #coord_cartesian(xlim = c(1.2,2.7),ylim=c(-2.7,-1.3))+
#   theme_bw() + 
#   theme(legend.position="bottom")
# p3
# #source("./sem/sem_is_w_inla.R")
# #is_w_inla_mod <- is.w.inla(data = df, init = init, prior.rho, 
# #                           dq.rho, rq.rho,fit.inla, N_0 = 800, N = 10000)
# #save(is_w_inla_mod, file = "./sem/sims/sem-is-w-inla.Rdata")
# #
# #source("./sem/sem_mcmc_w_inla.R")
# #mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = init$mu,
# #                               prior.rho, dq.rho, rq.rho, fit.inla,
# #                               n.samples = 100500, n.burnin = 500, n.thin = 10)
# #save(mcmc_w_inla_mod, file = "./sem/sims/sem-mcmc-w-inla.Rdata")
# 
# 
# spplot(turnout["TURNOUT01"],  col.regions= hcl.colors(16,palette =  hcl.pals("sequential")[21], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
# 
# spplot(turnout["GDPCAP"], col.regions= hcl.colors(16,palette =  hcl.pals("sequential")[21], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
