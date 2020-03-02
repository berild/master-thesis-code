# there are two different parametrizations in inla for the weibull
library(INLA)
variant  = 0

n = 1000

u = rep(rnorm(100,mean = 0, sd = 3), each = 10)

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

data = list(y=y, event=event, x=x, idx = rep(1:100,each = 10))

## the function inla.surv() creates a new dataset that takes into account the 
## censoring, the likelihood contribution of the data point changes 
## according to if we have censoring or not

formula = inla.surv(y,event)~ x + f(idx, model = "iid")
r=inla(formula,
       family ="weibullsurv",
       data=data,
       control.family = list(list(variant = variant)))
plot(r)

plot_gamma <- function(alpha, n){
  return(ggplot() + 
           geom_line(data = data.frame(x = seq(0,3*alpha*alpha,length.out = n),y = dgamma(seq(0,3*alpha*alpha,length.out = n),scale =alpha, shape = alpha)), aes(x = x, y = y)))
}
plot_gamma(alpha = 4, n =20)
