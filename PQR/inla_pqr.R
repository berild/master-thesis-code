library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
source("./PQR/pqr_general_functions.R")
source("./PQR/pqr_amis_w_inla.R")
source("./PQR/pqr_frq.R")
load("./PQR/pqr-models.Rdata")


GG_model <- function(mod, n = 200){
  x = runif(n)
  params = DXX[DXX$mod==mod,-1]
  mu = params[[1]] + params[[2]]*x
  sigma = exp(params[[3]] + params[[4]]*x)
  k = exp(params[[5]] + params[[6]]*x)
  theta = exp(mu)/(k^(sigma*sqrt(k)))
  beta = 1/(sigma*sqrt(k))
  scale = theta^beta
  y = (rgamma(n, shape = k, scale = scale))^(1/beta)
  return(list(
    data = data.frame(x = x, y = y),
    params = params,
    mod = data.frame(mu = mu, sigma = sigma, k = k, 
                     theta = theta,beta=beta)
  ))
}

gaussian_model <- function(mod, n = 200){
  x = runif(n)
  params = DXX[DXX$mod==mod,-1]
  mu = params[[1]] + params[[2]]*x
  tau = exp(params[[3]] + params[[4]]*x)
  y = rnorm(n,mean = mu, sd = sqrt(1/tau))
  return(list(
    data = data.frame(x = x, y = y),
    params = params,
    mod = data.frame(mu = mu, tau = tau)
  ))
  
}

mod = GG_model("D51",n=500)
init = list(mu = c(0,0),cov = 5*diag(2))
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
                              dq.param, rq.param, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 25,
                              quants = c(0.1,0.25,0.5,0.75,0.9))
amis_w_inla_mod$pqr_frq = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))

save(amis_w_inla_mod, file = "./sims/d51-pqr-amis-w-inla.Rdata")

mod = GG_model("D52",n=500)
init = list(mu = c(0,0),cov = 5*diag(2))
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
                              dq.param, rq.param, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 250)
amis_w_inla_mod$pqr_frq = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
save(amis_w_inla_mod, file = "./sims/d52-pqr-amis-w-inla.Rdata")

mod = GG_model("D53",n=500)
init = list(mu = c(0,0),cov = 5*diag(2))
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
                              dq.param, rq.param, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 250)
amis_w_inla_mod$pqr_frq = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
save(amis_w_inla_mod, file = "./sims/d53-pqr-amis-w-inla.Rdata")

mod = GG_model("D54",n=500)
init = list(mu = c(0,0),cov = 5*diag(2))
amis_w_inla_mod = amis.w.inla(data = mod$data, init = init, prior.param, 
                              dq.param, rq.param, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 250)
amis_w_inla_mod$pqr_frq = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
save(amis_w_inla_mod, file = "./sims/d54-pqr-amis-w-inla.Rdata")


data("ImmunogG")
init = list(mu = c(0,0),cov = 5*diag(2))
amis_w_inla_mod = amis.w.inla(data = data.frame(x = ImmunogG$Age,y = ImmunogG$IgG), init = init, prior.param, 
                              dq.param, rq.param, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 250)
amis_w_inla_mod$pqr_frq = PQR(ImmunogG$Age,ImmunogG$IgG, init =c(1, 0.05, -1,0, -3, 5),dom=c(0,7))
save(amis_w_inla_mod, file = "./sims/ImmunogG-pqr-amis-w-inla.Rdata")
