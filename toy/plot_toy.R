library(ggplot2)

load("./sims/toy-is-w-inla.Rdata")
load("./sims/toy-inla.Rdata")
load("./sims/toy-amis-w-inla.Rdata")
load("./sims/toy-mcmc-w-inla.Rdata")
source("./toy/general_functions.R")

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y, color = "MCMC with INLA")) +
  labs(color = "", x = "",y="",title=expression(alpha),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))


ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y, color = "IS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y, color = "MCMC with INLA")) +
  labs(color = "", x = "",y="",title=expression(tau),linetype = "")+
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
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                        kernel = "gaussian")[c(1,2)])
})

eta_joint_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1/length(mcmc_w_inla_mod$eta[,1]),length(mcmc_w_inla_mod$eta[,1])), n = 100, lims = c(0,2,-2,0))
amis_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))
is_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_is$x, y=eta_joint_kern_is$y), z=as.vector(eta_joint_kern_is$z))
mcmc_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_mcmc$x, y=eta_joint_kern_mcmc$y), z=as.vector(eta_joint_kern_mcmc$z))
ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "", x = "",y="",title=expression(delta^(1)),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "", x = "",y="",title=expression(delta^(2)),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) +
  labs(color = "",x=expression(delta^(1)),y=expression(delta^(2)),linetype="") +
  theme_bw() +
  coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) + 
  theme(legend.position="bottom")


ggplot() + 
  geom_path(data = data.frame(x = seq(length(mcmc_w_inla_mod$eta[,1])), y = mcmc_w_inla_mod$eta[,1]),aes(x=x,y=y)) + 
  geom_polygon(data = data.frame(x = inla_mod$marginals.fixed$x1[,1], 
                              y = 1500-inla_mod$marginals.fixed$x1[,2]/max(inla_mod$marginals.fixed$x1[,2])*1500), aes(x=y,y=x,fill = "target"),alpha = 0.5) +
  geom_path(data = data.frame(x = mcmc_w_inla_mod$eta_uni_kerns[[1]]$x,
                                 y = 1500 - mcmc_w_inla_mod$eta_uni_kerns[[1]]$y/max(mcmc_w_inla_mod$eta_uni_kerns[[1]]$y)*1500),aes(x=y,y=x,color = "estimate"))+
  scale_color_manual(values = "blue") + 
  scale_fill_manual(values = "red") + 
  labs(color = "",fill = "") + 
  theme_bw()

amis_adaptive = lapply(seq(13), function(x){
  tmp2 = c(0,125, seq(20,30)*5,5*31)
  tmp = rnorm(500,amis_w_inla_mod$theta$a.mu[x,1],amis_w_inla_mod$theta$a.cov[1,1,x])
  tmp = sort(tmp)
  y = dnorm(tmp,amis_w_inla_mod$theta$a.mu[x,1],amis_w_inla_mod$theta$a.cov[1,1,x])
  y = sum(tmp2[1:x]) + y/max(y)*tmp2[x+1]
  data.frame(x=tmp,y=y)
})
amis_adaptive2 = lapply(seq(13), function(x){
  tmp2 = c(0,125, seq(20,30)*5,5*31)
  tmp = inla_mod$marginals.fixed$x1[,1]
  y = inla_mod$marginals.fixed$x1[,2]
  y = sum(tmp2[1:x]) + y/max(y)*tmp2[x+1]
  data.frame(x=tmp,y=y)
})
ggplot() + 
  geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[8]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[9]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[10]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[11]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[12]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[13]],aes(x=y,y=x, color = "proposal")) + 
  geom_polygon(data = amis_adaptive2[[1]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[2]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[3]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[4]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[5]],aes(x=y,y=x,fill = "target"), alpha = 0.5) +
  geom_polygon(data = amis_adaptive2[[6]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[7]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[8]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[9]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[10]],aes(x=y,y=x,fill = "target"), alpha = 0.5) +
  geom_polygon(data = amis_adaptive2[[11]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[12]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  geom_polygon(data = amis_adaptive2[[13]],aes(x=y,y=x,fill = "target"), alpha = 0.5) + 
  scale_color_manual(values = "red") + 
  scale_fill_manual(values = "blue") + 
  labs(color="",fill="") +
  theme_bw()
  

