library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
source("./survival/frailty_general_functions.R")
source("./survival/frailty_amis_w_inla.R")
variant  = 0

n = 100

u = rep(rgamma(10,shape = 2, scale = 2), each = 10)

alpha = 1.1
beta = 2.2
x = c(scale(runif(n)))
eta = 1+beta*x + u
lambda = exp(eta)

y = rweibull(n,
             shape = alpha,
             scale = lambda^(-1/alpha))          
event = rep(1,n) 
n_class = 10
data = list(y=y, event=event, x=x, idx = rep(1:n_class,each = n/n_class))

init = list(mu = rep(1,n_class),cov = diag(n_class))

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty, 
                              dq.frailty, rq.frailty, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 250, kde = T)
amis_w_inla_mod$params = list(intercept = alpha, beta = beta)
save(amis_w_inla_mod,file = "./sims/test-frailty-amis-w-inla.Rdata")