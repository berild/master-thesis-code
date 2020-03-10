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
n = 100
n_class = 10
frailty.param = 3
u = rep(rgamma(n_class,shape = frailty.param, rate = frailty.param), each = n_class)

alpha = 1.1
beta = 2.2
x = c(scale(runif(n)))
eta = 1+beta*x + log(u)
lambda = exp(eta)

y = rweibull(n,
             shape = alpha,
             scale = lambda^(-1/alpha))          
event = rep(1,n) 
data = list(y=y, event=event, x=x, idx = rep(1:n_class,each = n/n_class))

init = list(mu = rep(1,n_class),cov = diag(n_class))

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty, 
                              dq.frailty, rq.frailty, fit.inla, 
                              N_t = seq(25,50,1)*10, N_0 = 250, kde = T)
amis_w_inla_mod$params = list(intercept = 1, beta = beta, alpha = frailty.param, params = u)
save(amis_w_inla_mod,file = "./sims/test1-frailty-amis-w-inla.Rdata")


# data(rats)
# rats = rats[rats$sex == "f",]
# amis_w_inla_mod = amis.w.inla(data = data.frame(y= rats$time, idx = rats$litter, 
#                                                 event = rats$status, rx = factor(rats$rx)), 
#                               init = init, prior.frailty, 
#                               dq.frailty, rq.frailty, fit.inla.k, 
#                               N_t = seq(25,26,1), N_0 = 25, kde = T)
# amis_w_inla_mod$params = kidney$frail
# save(amis_w_inla_mod,file = "./sims/kidney-frailty-amis-w-inla.Rdata")