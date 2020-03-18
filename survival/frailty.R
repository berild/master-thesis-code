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
intercept = 1
alpha = 1.1
beta = 2.2

# frailty = 1
frailty.param = 1
u = rep(rgamma(n_class,shape = frailty.param, rate = frailty.param), each = n/n_class)

x = c(scale(runif(n)))
eta = intercept + beta*x + log(u)
lambda = exp(eta)

y = rweibull(n,
             shape = alpha,
             scale = lambda^(-1/alpha))
event = rep(1,n)
data = list(y=y, event=event, x=x, idx = rep(1:n_class,each = n/n_class))


formula = inla.surv(y,event)~ x + f(idx, model = "iid")
res_inla =inla(formula,
               family ="weibullsurv",
               data=data,
               control.family = list(list(variant = variant)))

# init = list(mu = rep(1,n_class),cov = 1*diag(n_class))
max.frail = max(res_inla$marginals.hyperpar$`Precision for idx`[res_inla$marginals.hyperpar$`Precision for idx`[,2]>1e-8,1])
init = list(mu = res_inla$summary.random$idx[,2],
            cov = diag(res_inla$summary.random$idx[,3]^2))

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.effect,
                              dq.frailty, rq.frailty, fit.inla,
                              N_t = seq(25,50,1)*10, N_0 = 250,frailty=T)
amis_w_inla_mod$params = list(intercept = intercept, beta = beta, alpha = alpha, frailty = frailty.param, params = unique(u))
save(amis_w_inla_mod,file = "./sims/test3-frailty-amis-w-inla.Rdata")

# data(kidney)
# n_class = length(unique(kidney$id))
# init = list(mu = rep(1,n_class),cov = diag(n_class))
# amis_w_inla_mod = amis.w.inla(data = kidney, init = init, prior.frailty,
#                               dq.frailty, rq.frailty, fit.inla.kidney,
#                               N_t = seq(25,26,1), N_0 = 25, kde = T)
# save(amis_w_inla_mod,file = "./sims/rats-frailty-amis-w-inla.Rdata")