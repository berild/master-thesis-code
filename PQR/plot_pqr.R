library(ggplot2)

ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$a), aes(x=x,y=y)) +
  geom_vline(xintercept = mod$params[[1]])

ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$b), aes(x=x,y=y)) + 
  geom_vline(xintercept = mod$params[[2]])

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]],aes(x=x,y=y)) + 
  geom_vline(xintercept = mod$params[[3]])

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]],aes(x=x,y=y)) + 
  geom_vline(xintercept = mod$params[[4]])

fig_pqr <- ggplot()+
  geom_line(data = amis_w_inla_mod$pqr_truth, aes(x=x,y=y,color = quants, linetype = "Truth")) + 
  geom_line(data = amis_w_inla_mod$pqr, aes(x=x,y=y,color = quants, linetype = "INLA")) + 
  #geom_point(data = amis_w_inla_mod$data, aes(x = x, y = y))+
  labs(title="ImmunogG",x="x",y="",color = "quantiles", linetype = "method") + 
  theme_bw() 
fig_pqr

