library(ggplot2)
library(ggpubr)
load("./sims/sem/sem-amis-w-inla.Rdata")
load("./sims/sem/sem-is-w-inla.Rdata")
load("./sims/sem/sem-mcmc.Rdata")

p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_mod$margs$`(Intercept)`, aes(x=x,y=y,color = "MCMC")) + 
  labs(color = "",x=expression(beta[0]),y="") +
  theme_bw() +
  theme(legend.position="bottom")
p1


p2 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_mod$margs$`log(GDPCAP)`, aes(x=x,y=y,color = "MCMC")) + 
  labs(color = "",x=expression(beta[1]),y="") +
  theme_bw() +
  theme(legend.position="bottom")
p2

p3 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_mod$margs$tau, aes(x=x,y=y,color="MCMC")) +
  labs(color = "",x=expression(tau),y="") +
  theme_bw() +
  theme(legend.position="bottom")
p3

p4 <- ggplot() +
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  geom_contour(data = mcmc_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC"),bins = 6) +
  labs(color = "",x=expression(rho),y=expression(lambda),linetype="") +
  theme_bw() +
  theme(legend.position="bottom")
p4

ptot <- ggarrange(p1, p2, p3, p4,ncol=2, nrow=2, common.legend = T, legend="bottom")
ptot

sp1 <- spplot(turnout["TURNOUT01"],  col.regions= hcl.colors(40,palette = hcl.pals("sequential")[20], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
sp1
sp2 <- spplot(turnout["GDPCAP"], col.regions= hcl.colors(40,palette =  hcl.pals("sequential")[32], alpha = 1, rev = T),colorkey = list(space = "bottom"))
sp2
sptot <- ggarrange(sp1, sp2,ncol=2, nrow=1, common.legend = F)
sptot

