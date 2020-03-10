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
plot_gamma(shape = 3, rate = 3, n =200)


ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[1])+
  theme_bw()

ggplot() +   
  geom_line(data = amis_w_inla_mod$eta_kern[[2]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[2])+
  theme_bw()
  
ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[3]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[3])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[4]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[4])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[5]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[5])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[6]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[6])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[7]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[7])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[8]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[8])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[9]], aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$params[9])+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[10]], aes(x=x,y=y)) +
  geom_vline(xintercept = amis_w_inla_mod$params$params[10])+
  theme_bw()

amis_w_inla_mod$frailty = calc.param(amis_w_inla_mod$mlik,amis_w_inla_mod$eta,amis_w_inla_mod$weight)
amis_w_inla_mod$quants = kde.quantile(amis_w_inla_mod$eta_kern)

ggplot(amis_w_inla_mod$quants, aes(x=idx,y=q0.5)) + 
  geom_errorbar(aes(ymin=q0.025, ymax=q0.975), colour="black", width=.1) +
  geom_line() +
  geom_point(size=3) + 
  labs(y = "")+
  scale_x_discrete(name ="idx", 
                   limits=amis_w_inla_mod$quants$idx)+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$frailty,aes(x=x,y=y)) + 
  geom_vline(xintercept = 1,color = "firebrick") + 
  labs(title= "Rate/Shape of frailty") + 
  theme_bw()
