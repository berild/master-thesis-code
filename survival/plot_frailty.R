library(ggplot2)

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$intercept)

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta, aes(x=x,y=y)) +
  geom_vline(xintercept = amis_w_inla_mod$params$beta)

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]], aes(x=x,y=y,color = toString(1))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]], aes(x=x,y=y,color = toString(2))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[3]], aes(x=x,y=y,color = toString(3))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[4]], aes(x=x,y=y,color = toString(4))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[5]], aes(x=x,y=y,color = toString(5))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[6]], aes(x=x,y=y,color = toString(6))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[7]], aes(x=x,y=y,color = toString(7))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[8]], aes(x=x,y=y,color = toString(8))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[9]], aes(x=x,y=y,color = toString(9))) + 
  geom_line(data = amis_w_inla_mod$eta_kern[[10]], aes(x=x,y=y,color = toString(10))) +
  labs(color = "")

