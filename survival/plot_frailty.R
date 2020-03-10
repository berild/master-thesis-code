library(ggplot2)

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$intercept)+ 
  labs(title = "Intercept") +
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta, aes(x=x,y=y)) +
  geom_vline(xintercept = amis_w_inla_mod$params$beta)+
  labs(title= "Beta") + 
  theme_bw()

plot_gamma <- function(shape, rate, n){
  return(ggplot() + 
           geom_line(data = data.frame(x = seq(0,100,length.out = n),y = dgamma(seq(0,100,length.out = n),rate =rate, shape = shape)), aes(x = x, y = y)))
}
plot_gamma(shape = 1, rate = 0.01, n =200)


ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[3]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[4]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[5]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[6]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[7]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[8]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[9]], aes(x=x,y=y)) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[10]], aes(x=x,y=y)) +
  geom_vline(xintercept = 1,color = "firebrick") + 
  labs(title= "Rate/Shape of frailty") + 
  theme_bw()
  

amis_w_inla_mod$frailty = calc.param(amis_w_inla_mod$mlik,amis_w_inla_mod$eta,amis_w_inla_mod$weight)

ggplot() + 
  geom_line(data = amis_w_inla_mod$frailty,aes(x=x,y=y))
