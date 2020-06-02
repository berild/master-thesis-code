library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(SemiPar)

source("./PQR/general_functions.R")
source("./PQR/amis_w_inla.R")
source("./PQR/is_w_inla.R")
load("./PQR/pqr-models.Rdata")


ggamma__model <- function(mod, n = 200){
  x = runif(n)
  params = DXX[DXX$mod==mod,-1]
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

gaussian_model <- function(mod, n = 200){
  x = runif(n)
  params = DXX[DXX$mod==mod,-1]
  mu = params[[1]] + params[[2]]*x
  tau = exp(params[[3]] + params[[4]]*x)
  y = rnorm(n,mean = mu, sd = 1/sqrt(tau))
  return(list(
    data = data.frame(x = x, y = y),
    params = params
  ))
}

init = list(mu = c(0,0),cov = 4*diag(2))

# D51 data
mod = gaussian_model("D51",n=500)
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.gaussian,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gaussian",kde = T)
amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
amis_w_inla_mod$mod = mod
save(amis_w_inla_mod, file = "./sims/pqr-gaussian1-amis-w-inla.Rdata")
is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
                          dq.param, rq.param, fit.inla.gaussian,
                          N_0 = 800, N = 10000, pqr = "gaussian", kde = T)
is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
is_w_inla_mod$mod = mod
save(is_w_inla_mod, file = "./sims/pqr-gaussian1-is-w-inla.Rdata")

mod = ggamma__model("D51",n=500)
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.ggamma,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gamma",kde = T)
amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
amis_w_inla_mod$mod = mod
save(amis_w_inla_mod, file = "./sims/pqr-gamma1-amis-w-inla.Rdata")
is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
                          dq.param, rq.param, fit.inla.ggamma,
                          N_0 = 800, N = 10000, pqr = "gamma", kde = T)
is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
is_w_inla_mod$mod = mod
save(is_w_inla_mod, file = "./sims/pqr-gamma1-is-w-inla.Rdata")


# # D52 data
# mod = gaussian_model("D52",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
#                               dq.param, rq.param, fit.inla.gaussian,
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/pqr-gaussian2-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
#                           dq.param, rq.param, fit.inla.gaussian,
#                           N_0 = 800, N = 10000, pqr = "gaussian", kde = T)
# is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(is_w_inla_mod, file = "./sims/pqr-gaussian2-is-w-inla.Rdata")
# 
# mod = ggamma__model("D52",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
#                               dq.param, rq.param, fit.inla.ggamma,
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/pqr-gamma2-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
#                           dq.param, rq.param, fit.inla.ggamma,
#                           N_0 = 800, N = 10000, pqr = "gamma", kde = T)
# is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(is_w_inla_mod, file = "./sims/pqr-gamma2-is-w-inla.Rdata")
# 
# 
# # D53 data
# mod = gaussian_model("D53",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
#                               dq.param, rq.param, fit.inla.gaussian,
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/pqr-gaussian3-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
#                           dq.param, rq.param, fit.inla.gaussian,
#                           N_0 = 800, N = 10000, pqr = "gaussian", kde = T)
# is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(is_w_inla_mod, file = "./sims/pqr-gaussian3-is-w-inla.Rdata")
# 
# mod = ggamma__model("D53",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
#                               dq.param, rq.param, fit.inla.ggamma,
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/pqr-gamma3-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
#                           dq.param, rq.param, fit.inla.ggamma,
#                           N_0 = 800, N = 10000, pqr = "gamma", kde = T)
# is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(is_w_inla_mod, file = "./sims/pqr-gamma3-is-w-inla.Rdata")
# 
# 
# # D54 data
# mod = gaussian_model("D54",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
#                               dq.param, rq.param, fit.inla.gaussian,
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/pqr-gaussian4-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
#                           dq.param, rq.param, fit.inla.gaussian,
#                           N_0 = 800, N = 10000, pqr = "gaussian", kde = T)
# is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(is_w_inla_mod, file = "./sims/pqr-gaussian4-is-w-inla.Rdata")
# 
# mod = ggamma__model("D54",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param,
#                               dq.param, rq.param, fit.inla.ggamma,
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/pqr-gamma4-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = mod$data, init = init, prior.param,
#                           dq.param, rq.param, fit.inla.ggamma,
#                           N_0 = 800, N = 10000, pqr = "gamma", kde = T)
# is_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(is_w_inla_mod, file = "./sims/pqr-gamma4-is-w-inla.Rdata")
# 
# 
# # Lidar data
# data(lidar)
# prior.param <- function(x, log = TRUE) {
#   sum(dnorm(x, 0, 100, log = log))
# }
# 
# init = list(mu = c(14,0),cov = diag(c(3,0.0001)))
# 
# amis_w_inla_mod = amis.w.inla(data = lidar, init = init, prior.param,
#                               dq.param, rq.param, fit.inla,
#                               N_t = seq(25,50,1)*10, N_0 = 250,kde = T)
# save(amis_w_inla_mod, file = "./sims/pqr-gaussian-lidar-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = lidar, init = init, prior.param, 
#                           dq.param, rq.param, fit.inla,
#                           N_0 = 800, N = 10000, kde = T)
# save(amis_w_inla_mod, file = "./sims/pqr-gaussian-lidar-is-w-inla.Rdata")
# 
# # ImmunogG data
# data("ImmunogG")
# amis_w_inla_mod = amis.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG),
#                               init = init, prior.param, dq.param, rq.param, fit.inla.ggamma, 
#                               N_t = seq(25,50,1)*10, N_0 = 250, pqr = "gamma",kde = T)
# save(amis_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-amis-w-inla.Rdata")
# is_w_inla_mod = is.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG), 
#                           init = init, prior.param, dq.param, rq.param, fit.inla.ggamma,
#                           N_0 = 800, N = 10000, pqr = "gamma", kde = T)
# save(is_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-is-w-inla.Rdata")
