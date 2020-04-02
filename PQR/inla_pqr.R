library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
source("./PQR/pqr_general_functions.R")
source("./PQR/pqr_amis_w_inla.R")
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

init = list(mu = c(0,0),cov = 0.005*diag(2))


# D51 data
# mod = gaussian_model("D51",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.gaussian, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/d51-gaussian-pqr-amis-w-inla.Rdata")
# 
# mod = ggamma__model("D51",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.ggamma, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/d51-gamma-pqr-amis-w-inla.Rdata")
# 
# 
# # D52 data
# mod = gaussian_model("D52",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.gaussian, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/d52-gaussian-pqr-amis-w-inla.Rdata")
# 
# mod = ggamma__model("D52",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.ggamma, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/d52-gamma-pqr-amis-w-inla.Rdata")


# D53 data
# mod = gaussian_model("D53",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.gaussian, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/d53-gaussian-pqr-amis-w-inla.Rdata")
# 
# mod = ggamma__model("D53",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.ggamma, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/d53-gamma-pqr-amis-w-inla.Rdata")
# 
# 
# # D54 data
# mod = gaussian_model("D54",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.gaussian, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gaussian",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gaussian")
# save(amis_w_inla_mod, file = "./sims/d54-gaussian-pqr-amis-w-inla.Rdata")

# mod = ggamma__model("D54",n=500)
# amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.ggamma, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# amis_w_inla_mod$pqr_truth = pqr_truth_lines(mod$data,mod$params,type="gamma")
# save(amis_w_inla_mod, file = "./sims/d54-gamma-pqr-amis-w-inla.Rdata")
# 
# 
# # ImmunogG data
data("ImmunogG")
# amis_w_inla_mod = amis.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG), 
#                               init = init, prior.param, 
#                               dq.param, rq.param, fit.inla.ggamma, 
#                               N_t = seq(25,50,1)*10, N_0 = 250,
#                               pqr = "gamma",kde = T)
# save(amis_w_inla_mod, file = "./sims/ImmunogG-gamma-pqr-amis-w-inla.Rdata")




data(engel)
data = data.frame(x=engel$income,y=engel$foodexp)

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.param,
                              dq.param, rq.param, fit.inla.gaussian,
                              N_t = seq(25,50,1)*10, N_0 = 250,
                              pqr = "gaussian",kde = T)
save(amis_w_inla_mod, file = "./sims/both-northern-monthly-gamma-pqr-amis-w-inla.Rdata")

ggplot(data, aes(x=x,y=y)) +geom_point()
