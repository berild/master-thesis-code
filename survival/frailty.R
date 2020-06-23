library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(mvtnorm)
library(survival)
library(compositions)
source("./survival/frailty_general_functions.R")
source("./survival/frailty_amis_w_inla.R")

variant  = 0
n = 300
n_class = 100
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
init = list(mu = c(res_inla$summary.random$idx[,2],log(res_inla$summary.hyperpar[2,1])),
            cov = 1*diag(length(c(res_inla$summary.random$idx[,3],res_inla$summary.hyperpar[2,2])^2)))

## Frailty model
prior.frailty <- function(x, log = TRUE) {
  if (log){
    sum(dgamma(x[-length(x)],shape = x[length(x)],rate = x[length(x)],log = T)) + dgamma(x[length(x)],shape = 1,rate = 0.01,log = log)
  }else{
    prod(dgamma(x[-length(x)],shape = x[length(x)],rate = x[length(x)],log = F))*dgamma(x[length(x)],shape = 1,rate = 0.01,log = F)
  }
}

# dq.rho.lambda <- function(y, x, sigma = init$cov, log =TRUE) {
#   dmvt(y,sigma = sigma, df=3, delta = x, type = "shifted",log=log)
# }
# 
# rq.rho.lambda <- function(x, sigma = init$cov) {
#   as.vector(rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted"))
# }


dq.frailty <- function(y, x, sigma = init$cov, log =TRUE) {
  # dmvnorm(y, mean = x, sigma = sigma, log = log)
  #dmvt(y,sigma = sigma, df=3, delta = x, type = "shifted",log=log)
  # res <- dlnorm(y, meanlog = x, sdlog = sqrt(diag(sigma)), log = log)
  # if(log) {
  #   return(sum(res))
  # } else {
  #   return(prod(res))
  # }
  dlnorm.rplus(y,x,sigma)
}

rq.frailty <- function(x, sigma = init$cov) {
  #abs(as.vector(rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted")))
  # abs(as.vector(rmvnorm(1, mean = x, sigma = sigma)))
  # rlnorm(length(x), meanlog = x, sdlog = sqrt(diag(sigma)))
  as.vector(rlnorm.rplus(1,x,sigma))
}

amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty,
                              dq.frailty, rq.frailty, fit.inla,
                              N_t = rep(500,20), N_0 = 5000,frailty=T)
amis_w_inla_mod$params = list(intercept = intercept, beta = beta, alpha = alpha, frailty = frailty.param, params = unique(u))
save(amis_w_inla_mod,file = "./sims/frailty-100-amis-w-inla.Rdata")

# library(survival)
# data(kidney)
# n_class = length(unique(kidney$id))
# # formula = inla.surv(time,status)~ age + sex + disease + f(id, model = "iid")
# # res_inla =inla(formula,
# #                family ="weibullsurv",
# #                data=kidney,
# #                control.family = list(list(variant = variant)))
# init = list(mu = c(kidney$frail[seq(1,nrow(kidney)/2)*2],1), cov = diag(n_class+1))
# amis_w_inla_mod = amis.w.inla(data = kidney, init = init, prior.frailty,
#                               dq.frailty, rq.frailty, fit.inla.kidney,
#                               N_t = rep(500,20), N_0 = 5000,frailty=T)
# save(amis_w_inla_mod,file = "./sims/kidney-frailty-amis-w-inla.Rdata")