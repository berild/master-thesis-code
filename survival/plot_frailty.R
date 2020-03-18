library(ggplot2)
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


ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,color = "AMIS w/ INLA")) + 
  geom_line(data = as.data.frame(res_inla$marginals.fixed[[1]]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$intercept,color = "firebrick")+ 
  labs(title = "Intercept") +
  labs(color = "")+
  theme_bw()


ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta, aes(x=x,y=y,color = "AMIS w/ INLA")) +
  geom_line(data = as.data.frame(res_inla$marginals.fixed[[2]]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$beta,color = "firebrick")+
  labs(title= "Beta") + 
  labs(color = "")+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$alpha, aes(x=x,y=y,color = "AMIS w/ INLA")) +
  geom_line(data = as.data.frame(res_inla$marginals.hyperpar[[1]]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$alpha,color = "firebrick")+
  labs(title= "Alpha") + 
  labs(color = "")+
  theme_bw()


#amis_w_inla_mod$frailty = calc.param(amis_w_inla_mod$mlik,amis_w_inla_mod$eta,amis_w_inla_mod$weight)
#amis_w_inla_mod$quants = kde.quantile(amis_w_inla_mod$eta_kern)
inla_frailty_idx = as.data.frame(res_inla$summary.random$idx[,4:6])
names(inla_frailty_idx) = names(amis_w_inla_mod$frailty_idx[,2:4])
inla_frailty_idx = data.frame(inla_frailty_idx,idx = amis_w_inla_mod$frailty_idx$idx)

amis_w_inla_mod$frailty_idx = data.frame(amis_w_inla_mod$frailty_idx, rv = amis_w_inla_mod$params$params)

ggplot() + 
  geom_ribbon(data = inla_frailty_idx,aes(x =idx,ymin = q0.025, ymax=q0.975),fill = gg_color_hue(3)[2],color = gg_color_hue(3)[2] ,alpha = 0.3,show.legend = FALSE)  + 
  geom_ribbon(data = amis_w_inla_mod$frailty_idx,aes(x =idx, ymin=q0.025, ymax=q0.975),fill = gg_color_hue(3)[1],color = gg_color_hue(3)[1], alpha = 0.3,show.legend = FALSE) +
  geom_line(data = amis_w_inla_mod$frailty_idx,aes(x =idx,y=q0.5,color = "AMIS w/INLA")) +
  geom_line(data = inla_frailty_idx, aes(x =idx,y=q0.5,color = "INLA")) +
  geom_point(data = amis_w_inla_mod$frailty_idx,aes(x =idx,y=q0.5,color = "AMIS w/INLA"),size=2) +
  geom_point(data = inla_frailty_idx, aes(x =idx,y=q0.5,color = "INLA")) +
  geom_point(data = amis_w_inla_mod$frailty_idx,aes(x =idx,y = rv,color = "Real Value"),size=2) + 
  labs(y = "",color = "") +
  scale_x_discrete(name ="idx", 
                   limits=amis_w_inla_mod$frailty_idx$idx)+
  theme_bw()

ggplot() + 
  geom_line(data = amis_w_inla_mod$frailty,aes(x=x,y=y,color = "AMIS w/ INLA")) + 
  geom_line(data = as.data.frame(res_inla$marginals.hyperpar$`Precision for idx`[-c(74,75),]),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$params$frailty,color = "firebrick") + 
  labs(title= "Rate/Shape of frailty") + 
  labs(color = "")+
  theme_bw()

