library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
library(mvtnorm)
library(survival)
library(compositions)
# loading general functions
source("./survival/frailty_general_functions.R")
# loading amis with inla functions
source("./survival/frailty_amis_w_inla.R")

# setting initial parameters for frailty model
variant  = 0
n = 300
n_class = 100
intercept = 1
alpha = 1.1
beta = 2.2

# frailty = 1
frailty.param = 1
# generating frailties
u = rep(rgamma(n_class,shape = frailty.param, rate = frailty.param), each = n/n_class)

# generating weibull data
x = c(scale(runif(n)))
eta = intercept + beta*x + log(u)
lambda = exp(eta)

y = rweibull(n,
             shape = alpha,
             scale = lambda^(-1/alpha))
event = rep(1,n)
data = list(y=y, event=event, x=x, idx = rep(1:n_class,each = n/n_class))

# fitting log normal model frailties with INLA
formula = inla.surv(y,event)~ x + f(idx, model = "iid")
res_inla =inla(formula,
               family ="weibullsurv",
               data=data,
               control.family = list(list(variant = variant)))

max.frail = max(res_inla$marginals.hyperpar$`Precision for idx`[res_inla$marginals.hyperpar$`Precision for idx`[,2]>1e-8,1])

# setting intial proposal distributoin for amis with inla using the approximation provided by inla
init = list(mu = c(res_inla$summary.random$idx[,2],log(res_inla$summary.hyperpar[2,1])),
            cov = 1*diag(length(c(res_inla$summary.random$idx[,3],res_inla$summary.hyperpar[2,2])^2)))

# prior distribution for frailties
prior.frailty <- function(x, log = TRUE) {
  if (log){
    sum(dgamma(x[-length(x)],shape = x[length(x)],rate = x[length(x)],log = T)) + dgamma(x[length(x)],shape = 1,rate = 0.01,log = log)
  }else{
    prod(dgamma(x[-length(x)],shape = x[length(x)],rate = x[length(x)],log = F))*dgamma(x[length(x)],shape = 1,rate = 0.01,log = F)
  }
}

# proposal distribution in amis with inla
dq.frailty <- function(y, x, sigma = init$cov, log =TRUE) {
  dlnorm.rplus(y,x,sigma)
}
rq.frailty <- function(x, sigma = init$cov) {
  as.vector(rlnorm.rplus(1,x,sigma))
}

# function fitting conditional models with inla
fit.inla <- function(data,eta){
  # data$oset = log(eta[data$idx])
  eta = eta[-length(eta)]
  data$oset = log(eta[data$idx])
  formula = inla.surv(y,event) ~ x + offset(oset)
  res=inla(formula,
           family ="weibullsurv",
           data=data,
           control.family = list(list(variant = variant)))
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           beta = res$marginals.fixed[[2]],
                           alpha = res$marginals.hyperpar[[1]])))
}


# running amis with inla algorithm
amis_w_inla_mod = amis.w.inla(data = data, init = init, prior.frailty,
                              dq.frailty, rq.frailty, fit.inla,
                              N_t = rep(500,20), N_0 = 5000,frailty=T)
amis_w_inla_mod$params = list(intercept = intercept, beta = beta, alpha = alpha, frailty = frailty.param, params = unique(u))
save(amis_w_inla_mod,file = "./sims/frailty-100-amis-w-inla.Rdata")

# save(amis_w_inla_mod,file = "./sims/kidney-frailty-amis-w-inla.Rdata")