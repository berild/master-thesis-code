library(ggplot2)

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$params$intercept)

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta, aes(x=x,y=y)) +
  geom_vline(xintercept = amis_w_inla_mod$params$beta)


ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern, aes(x=x,y=y)) + 
  labs(color = "")
  