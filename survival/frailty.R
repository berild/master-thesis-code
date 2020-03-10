library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(survival)
source("./survival/frailty_general_functions.R")
source("./survival/frailty_amis_w_inla.R")

# variant  = 0
# n = 100
# n_class = 10
# frailty.param = 3
# u = rep(rgamma(n_class,shape = frailty.param, rate = frailty.param), each = n_class)
# 
# alpha = 1.1
# beta = 2.2
# x = c(scale(runif(n)))
# eta = 1+beta*x + log(u)
# lambda = exp(eta)
# 
# y = rweibull(n,
#              shape = alpha,
#              scale = lambda^(-1/alpha))          
# event = rep(1,n) 
# data = list(y=y, event=event, x=x, idx = rep(1:n_class,each = n/n_class))
# 
# init = list(mu = rep(1,n_class),cov = diag(n_class))
# 
# amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty, 
#                               dq.frailty, rq.frailty, fit.inla, 
#                               N_t = seq(25,50,1)*10, N_0 = 250, kde = T)
# amis_w_inla_mod$params = list(intercept = 1, beta = beta, alpha = frailty.param, params = unique(u))
# save(amis_w_inla_mod,file = "./sims/test1-frailty-amis-w-inla.Rdata")
# 

data(rats)
n_class = unique(rats$litter)
init = list(mu = rep(1,n_class),cov = diag(n_class))
amis_w_inla_mod = amis.w.inla(data = rats, init = init, prior.frailty, 
                              dq.frailty, rq.frailty, fit.inla.rats, 
                              N_t = seq(25,50,1)*10, N_0 = 250, kde = T)
save(amis_w_inla_mod,file = "./sims/rats-frailty-amis-w-inla.Rdata")