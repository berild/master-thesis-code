library(ggplot2)

ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$intercept), aes(x=x,y=y)) +
  geom_vline(xintercept = amis_w_inla_mod$pqr_frq$ml[1])

ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$beta), aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$pqr_frq$ml[2])

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern$f,aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$pqr_frq$ml[5])

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern$g,aes(x=x,y=y)) + 
  geom_vline(xintercept = amis_w_inla_mod$pqr_frq$ml[6])

fig_pqr <- ggplot()+
  geom_line(data = amis_w_inla_mod$pqr_frq$quantiles, aes(x=x,y=y,color = quants, linetype = "PQR")) + 
  geom_line(data = amis_w_inla_mod$pqr, aes(x=x,y=y,color = quants, linetype = "INLA")) + 
  labs(title="ImmunogG",x="x",y="",color = "quantiles", linetype = "method") + 
  theme_bw() 
fig_pqr

