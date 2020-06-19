library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(survival)
library(LaplacesDemon)
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
init = list(mu = c(res_inla$summary.random$idx[,2],res_inla$summary.hyperpar[2,1]),
            cov = diag(c(res_inla$summary.random$idx[,3],res_inla$summary.hyperpar[2,2])^2))

## Frailty model
prior.frailty <- function(x, log = TRUE) {
  if (log){
    sum(dgamma(exp(x[-length(x)]),shape = x[length(x)],rate = x[length(x)],log = T)) + dgamma(x[length(x)],shape = 1,rate = 0.01,log = log)
  }else{
    prod(dgamma(exp(x[-length(x)]),shape = x[length(x)],rate = x[length(x)],log = F))*dgamma(x[length(x)],shape = 1,rate = 0.01,log = F)
  }
}



dq.frailty <- function(y, x, sigma = init$cov, log =TRUE) {
  tmp = dst(y[length(y)],mu = x[length(x)],sigma=sqrt(diag(sigma)[length(x)]),nu=3)
  sum(dgamma(y[-length(y)],shape = x[length(x)],rate = x[length(x)],log=log)) + tmp
}

rq.frailty <- function(x, sigma = init$cov) {
  tmp  = abs(rst(1,mu = x[length(x)],sigma=sqrt(diag(sigma)[length(x)]),nu=3))
  c(rgamma(length(x)-1,shape = tmp,rate = tmp),tmp)
}

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty,
                              dq.frailty, rq.frailty, fit.inla,
                              N_t = seq(30,50,1)*10, N_0 = 4200,frailty=T)
amis_w_inla_mod$params = list(intercept = intercept, beta = beta, alpha = alpha, frailty = frailty.param, params = unique(u))
save(amis_w_inla_mod,file = "./sims/test-last-frailty-amis-w-inla.Rdata")

# data(kidney)
# n_class = length(unique(kidney$id))
# init = list(mu = rep(1,n_class),cov = diag(n_class))
# amis_w_inla_mod = amis.w.inla(data = kidney, init = init, prior.frailty,
#                               dq.frailty, rq.frailty, fit.inla.kidney,
#                               N_t = seq(25,26,1), N_0 = 25, kde = T)
# save(amis_w_inla_mod,file = "./sims/rats-frailty-amis-w-inla.Rdata")