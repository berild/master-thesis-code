# Lidar data PQR
library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)
library(spatstat)
source("./PQR/pqr_general_functions.R")
source("./PQR/pqr_amis_w_inla.R")
library(SemiPar)

data(lidar)

formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  T)

mod = inla(formula,  data = lidar,
           control.predictor = list(compute = T))


fit.inla <- function(data,eta){
  formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
  res = inla(formula, 
             data = data,
             control.predictor = list(compute = T),
             scale = exp(eta[1] + eta[2]*data$range), 
             family = "gaussian",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           prec.range = res$marginals.hyperpar[[1]])))
}

prior.param <- function(x, log = TRUE) {
  sum(dnorm(x, 0, 100, log = log))
}

init = list(mu = c(14,0),cov = diag(c(3,0.0001)))

amis_w_inla_mod = amis.w.inla(data = lidar, init = init, prior.param,
                              dq.param, rq.param, fit.inla,
                              N_t = seq(25,50,1)*10, N_0 = 250,kde = T)
save(amis_w_inla_mod, file = "./sims/lidar-gaussian-pqr-amis-w-inla.Rdata")


# formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
# res = inla(formula, 
#            data = lidar,
#            control.predictor = list(compute = T),
#            scale = exp(amis_w_inla_mod$theta$a.mu[8,1] + amis_w_inla_mod$theta$a.mu[8,2]*lidar$range), 
#            control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
#            verbose = FALSE)
# 
# 
# plot(lidar$range,lidar$logratio)
# lines(lidar$range,res$summary.linear.predictor[,1],col=2)
# #lines(lidar$range,mod$summary.linear.predictor[,1],col=4)
# lines(lidar$range,res$summary.linear.predictor[,3],col=3)
# #lines(lidar$range,mod$summary.linear.predictor[,3],col=5)
# lines(lidar$range,res$summary.linear.predictor[,5],col=3)
# #lines(lidar$range,mod$summary.linear.predictor[,5],col=5)
# 
# tmp = data.frame(x = NA, y=NA)
# for (i in seq(6,nrow(lidar), 5)){
#   tmp = rbind(tmp,c(mean(lidar$range[(i-5):i]),log(1/var(lidar$logratio[(i-5):i]))))
# }
# tmp = tmp[-1,]
# tmp
# tmp_res = lm(y~x, data = tmp)
# tmp_res$coefficients
# formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
# res_ml = inla(formula, 
#            data = lidar,
#            control.predictor = list(compute = T),
#            scale = exp(tmp_res$coefficients[[1]] + tmp_res$coefficients[[2]]*lidar$range), 
#            control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
#            verbose = FALSE)
# 
# res = inla(formula, 
#               data = lidar,
#               control.predictor = list(compute = T),
#               scale = exp(amis_w_inla_mod$theta$a.mu[8,1]+ amis_w_inla_mod$theta$a.mu[8,2]*lidar$range), 
#               control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
#               verbose = FALSE)
# 
# plot(lidar$range,lidar$logratio)
# lines(lidar$range,res$summary.linear.predictor[,1],col=2)
# lines(lidar$ranger,res_ml$summary.linear.predictor[,1],col=4)
# lines(lidar$range,res$summary.linear.predictor[,3],col=3)
# lines(lidar$range,res_ml$summary.linear.predictor[,3],col=5)
# lines(lidar$range,res$summary.linear.predictor[,5],col=3)
# lines(lidar$range,res_ml$summary.linear.predictor[,5],col=5)

