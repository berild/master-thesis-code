library(ggplot2)
library(ggpubr)
# simulated datasets
load("./sims/test1-frailty-amis-w-inla.Rdata") # frailty = 1
load("./sims/test2-frailty-amis-w-inla.Rdata") # frailty = 2
load("./sims/test3-frailty-amis-w-inla.Rdata") # frailty = 3

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

formula = inla.surv(y,event)~ x + f(idx, model = "iid")
res_inla =inla(formula,
       family ="weibullsurv",
       data=amis_w_inla_mod$data,
       control.family = list(list(variant = 0)))


p1 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y),color = gg_color_hue(1)) + 
  #geom_line(data = as.data.frame(res_inla$marginals.fixed[[1]]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$intercept,color = "firebrick")+ 
  labs(x = expression(beta[0]),y="") +
  theme_bw()
p1

p2 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta, aes(x=x,y=y),color = gg_color_hue(1)) +
  #geom_line(data = as.data.frame(res_inla$marginals.fixed[[2]]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$beta,color = "firebrick")+
  labs(x= expression(beta[1]),y="") + 
  theme_bw()
p2

p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$alpha, aes(x=x,y=y),color = gg_color_hue(1)) +
  #geom_line(data = as.data.frame(res_inla$marginals.hyperpar[[1]]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$alpha,color = "firebrick")+
  labs(x= expression(alpha),y="") + 
  theme_bw()
p3

#amis_w_inla_mod$quants = kde.quantile(amis_w_inla_mod$eta_kern)
#inla_frailty_idx = as.data.frame(res_inla$summary.random$idx[,4:6])
#names(inla_frailty_idx) = names(amis_w_inla_mod$frailty_idx[,2:4])
#inla_frailty_idx = data.frame(inla_frailty_idx,idx = amis_w_inla_mod$frailty_idx$idx[-11])

amis_w_inla_mod$frailty_idx = data.frame(amis_w_inla_mod$frailty_idx, rv = log(amis_w_inla_mod$params$params))


p5 <- ggplot() + 
  #geom_ribbon(data = inla_frailty_idx,aes(x =idx,ymin = q0.025, ymax=q0.975),fill = gg_color_hue(3)[2],color = gg_color_hue(3)[2] ,alpha = 0.3,show.legend = FALSE)  + 
  geom_ribbon(data = amis_w_inla_mod$frailty_idx,aes(x =idx, ymin=q0.025, ymax=q0.975),fill = gg_color_hue(3)[1],color = gg_color_hue(3)[1], alpha = 0.3,show.legend = FALSE) +
  geom_line(data = amis_w_inla_mod$frailty_idx,aes(x =idx,y=q0.5,color = "AMIS with INLA")) +
  #geom_line(data = inla_frailty_idx, aes(x =idx,y=q0.5,color = "INLA")) +
  geom_point(data = amis_w_inla_mod$frailty_idx,aes(x =idx,y=q0.5,color = "AMIS with INLA"),size=2) +
  #geom_point(data = inla_frailty_idx, aes(x =idx,y=q0.5,color = "INLA")) +
  geom_point(data = amis_w_inla_mod$frailty_idx,aes(x =idx,y = rv,color = "Real Value"),size=2) + 
  labs(y = "",color = "") +
  scale_color_manual(values = c(gg_color_hue(2)[1],"firebrick")) + 
  scale_x_discrete(name ="idx", 
                   limits=amis_w_inla_mod$frailty_idx$idx)+
  theme_bw() +
  theme(legend.position="bottom")
p5
amis_w_inla_mod$frailty.prec.kde= as.data.frame(density(x = amis_w_inla_mod$frailty.prec,
                                                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                                                        from = 0,
                                                        to = 10,
                                                        kernel = "gaussian")[c(1,2)])
p6 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$frailty.prec.kde,aes(x=x,y=y),color = gg_color_hue(1)) + 
  #geom_line(data = as.data.frame(res_inla$marginals.hyperpar$`Precision for idx`[-c(74,75),]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$frailty,color = "firebrick") + 
  labs(y="",x=expression(gamma))+
  coord_cartesian(xlim = c(0,10)) +
  theme_bw()
p6

ggarrange(p1,p2,p3,p6,nrow=2,ncol=2)
