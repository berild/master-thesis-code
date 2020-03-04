library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(survival)
source("./survival/frailty_general_functions.R")
source("./survival/frailty_amis_w_inla.R")

variant  = 0
data(rats)
n = 100

u = rep(rgamma(10,shape = 1, rate = 1), each = 10)

alpha = 1.1
beta = 2.2
x = c(scale(runif(n)))
eta = 1+beta*x + log(u)
lambda = exp(eta)

y = rweibull(n,
             shape = alpha,
             scale = lambda^(-1/alpha))          
event = rep(1,n) 
n_class = 10
data = list(y=y, event=event, x=x, idx = rep(1:n_class,each = n/n_class))

init = list(mu = 1,cov = 1)

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty, 
                              dq.frailty, rq.frailty, fit.inla, 
                              N_t = seq(25,30,1)*10, N_0 = 250, kde = T)
amis_w_inla_mod$params = list(intercept = 1, beta = beta)
save(amis_w_inla_mod,file = "./sims/test-frailty-amis-w-inla.Rdata")

