library(ggplot2)
load("./sims/disease-is-w-inla.Rdata")
load("./sims/disease-amis-w-inla.Rdata")

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$rf1, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$rf1, aes(x=x,y=y, color = "IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(alpha[1]),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$rf2, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$rf2, aes(x=x,y=y, color = "IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(alpha[2]),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$rf3, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$rf3, aes(x=x,y=y, color = "IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(alpha[3]),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau1, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau1, aes(x=x,y=y, color = "IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(tau[s]),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))


ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau2, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau2, aes(x=x,y=y, color = "IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(tau[v]),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

amis_w_inla_mod$eta_uni_kerns = lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
is_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
eta_joint_kern_amis1 = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(0,1,0,1))
eta_joint_kern_amis2 = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,3], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(0,1,0,1))
eta_joint_kern_amis3 = kde2d.weighted(x = amis_w_inla_mod$eta[,2], y = amis_w_inla_mod$eta[,3], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(0,1,0,1))
eta_joint_kern_is1 = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(0,1,0,1))
eta_joint_kern_is2 = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,3], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(0,1,0,1))
eta_joint_kern_is3 = kde2d.weighted(x = is_w_inla_mod$eta[,2], y = is_w_inla_mod$eta[,3], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(0,1,0,1))
amis_w_inla_mod$eta_joint_kern = list(eta12 = data.frame(expand.grid(x=eta_joint_kern_amis1$x, y=eta_joint_kern_amis1$y), z=as.vector(eta_joint_kern_amis1$z)),
                                      eta13 = data.frame(expand.grid(x=eta_joint_kern_amis2$x, y=eta_joint_kern_amis2$y), z=as.vector(eta_joint_kern_amis2$z)),
                                      eta23 = data.frame(expand.grid(x=eta_joint_kern_amis3$x, y=eta_joint_kern_amis3$y), z=as.vector(eta_joint_kern_amis3$z)))
is_w_inla_mod$eta_joint_kern = list(eta12 = data.frame(expand.grid(x=eta_joint_kern_is1$x, y=eta_joint_kern_is1$y), z=as.vector(eta_joint_kern_is1$z)),
                                    eta13 = data.frame(expand.grid(x=eta_joint_kern_is2$x, y=eta_joint_kern_is2$y), z=as.vector(eta_joint_kern_is2$z)),
                                    eta23 = data.frame(expand.grid(x=eta_joint_kern_is3$x, y=eta_joint_kern_is3$y), z=as.vector(eta_joint_kern_is3$z)))
ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(delta^(1)),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(delta^(2)),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[3]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[3]], aes(x=x,y=y,color="IS with INLA")) +
  labs(color = "", x = "",y="",title=expression(delta^(3)),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern$eta12, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern$eta12, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  labs(color = "",x=expression(delta^(1)),y=expression(delta^(2)),linetype="") +
  theme_bw() +
  coord_cartesian(xlim = c(0,1),ylim=c(0,1)) + 
  theme(legend.position="bottom")

ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern$eta13, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern$eta13, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  labs(color = "",x=expression(delta^(1)),y=expression(delta^(3)),linetype="") +
  theme_bw() +
  coord_cartesian(xlim = c(0,1),ylim=c(0,1)) + 
  theme(legend.position="bottom")

ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern$eta23, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern$eta23, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  labs(color = "",x=expression(delta^(2)),y=expression(delta^(3)),linetype="") +
  theme_bw() +
  coord_cartesian(xlim = c(0,1),ylim=c(0,1)) + 
  theme(legend.position="bottom")
