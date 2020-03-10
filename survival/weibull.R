# there are two different parametrizations in inla for the weibull
library(INLA)
variant  = 0

n = 100

u = rep(rnorm(10,mean = 0, sd = 3), each = 10)

alpha = 1.1
beta = 2.2
x = c(scale(runif(n)))
eta = 1+beta*x + u
lambda = exp(eta)

y = rweibull(n,
               shape = alpha,
               scale = lambda^(-1/alpha))
# this is the censoring, we set it to 1 so all 
#failure times are observed             
event = rep(1,n) 

data = list(y=y, event=event, x=x, idx = rep(1:10,each = 10))

## the function inla.surv() creates a new dataset that takes into account the 
## censoring, the likelihood contribution of the data point changes 
## according to if we have censoring or not

formula = inla.surv(y,event)~ x + f(idx, model = "iid")
r=inla(formula,
       family ="weibullsurv",
       data=data,
       control.family = list(list(variant = variant)))
plot(r)

fit.inla <- function(data,eta){
  oset = eta[data$idx]
  formula = inla.surv(y,event) ~ 1 + x + offset(log(oset))
  res=inla(formula,
           family ="weibullsurv",
           data=data,
           control.family = list(list(variant = variant)))
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           beta = res$marginals.fixed[[2]])))
}

plot_gamma <- function(shape, rate, n){
  require(ggplot2)
  return(ggplot() + 
           geom_line(data = data.frame(x = seq(0,100,length.out = n),y = dgamma(seq(0,100,length.out = n),rate =rate, shape = shape)), aes(x = x, y = y)))
}
plot_gamma(shape = 1, rate = 0.01, n =200)
