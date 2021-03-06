library(ggplot2)
library(INLA)
library(ggpubr)
source("./PQR/general_functions.R")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_temp = gg_color_hue(4)

p1 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$a), aes(x=x,y=y, color = "AMIS with INLA")) + 
  #geom_line(data = as.data.frame(is_w_inla_mod$margs$a),aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data = as.data.frame(mcmc_w_inla_mod$margs$a), aes(x=x,y=y, color = "MCMC with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$a) + 
  labs(color = "",x="a",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p1

p2 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$b), aes(x=x,y=y, color = "AMIS with INLA")) +
  geom_line(data = as.data.frame(mcmc_w_inla_mod$margs$b), aes(x=x,y=y, color = "MCMC with INLA")) + 
  #geom_line(data = as.data.frame(is_w_inla_mod$margs$b),aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$b) + 
  labs(color = "",x="b",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p2
p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]],aes(x=x,y=y,color="MCMC with INLA")) + 
  #geom_line(data = is_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$c) + 
  labs(color = "",x = "c",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]],aes(x=x,y=y,color="MCMC with INLA")) + 
  #geom_line(data = is_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$d) + 
  labs(color = "",x = "d",y ="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p4

find_post_stats_is <- function(mod){
  #tmp = running.ESS(mod$eta,mod$times,mod$weight,norm=T)
  data.frame(
    a = c(inla.zmarginal(mod$margs$a,silent=T)[[1]],inla.zmarginal(mod$margs$a,silent=T)[[2]]),
    b = c(inla.zmarginal(mod$margs$b,silent=T)[[1]],inla.zmarginal(mod$margs$b,silent=T)[[2]]),
    c = c(mod$theta$a.mu[28,1],sqrt(mod$theta$a.cov[1,1,28])),
    d= c(mod$theta$a.mu[28,2],sqrt(mod$theta$a.cov[2,2,28])),
    ess = c(tmp$ess[nrow(tmp)],tmp$time[nrow(tmp)])
  )
}

find_post_stats_is(amis_w_inla_mod)
create_pqr_text <- function(pqr,type=1){
  pqr_text = data.frame(x = c(),y = c(),text = c())
  for (x in unique(pqr$quants)){
    if (type == 1){
      tmp.idx = which.max(pqr$x[pqr$quants==x])
      tmpx = pqr$x[pqr$quants==x][tmp.idx]+0.03
      tmpy = pqr$y[pqr$quants==x][tmp.idx]
    }else{
      tmp.idx = which.min(pqr$x[pqr$quants==x])
      tmpx = pqr$x[pqr$quants==x][tmp.idx]-0.03
      tmpy = pqr$y[pqr$quants==x][tmp.idx]
    }
    pqr_text = rbind(pqr_text,data.frame(x=tmpx,y=tmpy,text = x)) 
  }
  return(pqr_text)
}

load("./sims/pqr/pqr-m1-gaussian-amis-w-inla.Rdata")
#load("./sims/pqr/pqr-m1-gaussian-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,1)
p1 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  #geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label="M1 Gaussian",x = 0.5, y = 10) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p1

load("./sims/pqr/pqr-m2-gaussian-amis-w-inla.Rdata")
#load("./sims/pqr/pqr-m2-gaussian-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,2)
p2 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  #geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label="M2 Gaussian",x = 0.5, y = 2.5) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p2

load("./sims/pqr/pqr-m1-gamma-amis-w-inla.Rdata")
#load("./sims/pqr/pqr-m1-gamma-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,2)
p3 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  #geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label="M1 Gamma",x=0.5,y=70) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw() + 
  coord_cartesian(ylim = c(0,80))
p3


load("./sims/pqr/pqr-m2-gamma-amis-w-inla.Rdata")
#load("./sims/pqr/pqr-m2-gamma-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,1)
p4 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  #geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label="M2 Gamma",x = 0.5, y = 0.47) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p4


ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = T,legend="bottom")

# Lidar
load("./sims/pqr/pqr-gaussian-lidar-is-w-inla.Rdata")
load("./sims/pqr/pqr-gaussian-lidar-mcmc-w-inla.Rdata")
load("./sims/pqr/pqr-gaussian-lidar-amis-w-inla.Rdata")
domain = matrix(c(11,17,-0.022,-0.008),nrow=2)
amis_w_inla_mod$theta$a.mu[,2]=amis_w_inla_mod$theta$a.mu[,2]/amis_w_inla_mod$scale[1]
amis_w_inla_mod$params = amis_w_inla_mod$theta$a.mu[nrow(amis_w_inla_mod$theta$a.mu),]
amis_w_inla_mod$eta[,2] = amis_w_inla_mod$eta[,2]/amis_w_inla_mod$scale[1]
amis_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        from = domain[1,x],
                        to = domain[2,x],
                        kernel = "gaussian")[c(1,2)])
})
amis_w_inla_mod$mod$range = amis_w_inla_mod$mod$range*amis_w_inla_mod$scale[1]

is_w_inla_mod$eta[,2] = is_w_inla_mod$eta[,2]/is_w_inla_mod$scale[1]
is_w_inla_mod$params = calc.post.mean(is_w_inla_mod$eta,is_w_inla_mod$weight)
is_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        from = domain[1,x],
                        to = domain[2,x],
                        kernel = "gaussian")[c(1,2)])
})
is_w_inla_mod$mod$range = is_w_inla_mod$mod$range*is_w_inla_mod$scale[1]

mcmc_w_inla_mod$eta[,2] = mcmc_w_inla_mod$eta[,2]/mcmc_w_inla_mod$scale[1]
mcmc_w_inla_mod$params = calc.post.mean(mcmc_w_inla_mod$eta)
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                        from = domain[1,x],
                        to = domain[2,x],
                        kernel = "gaussian")[c(1,2)])
})
mcmc_w_inla_mod$mod$range = amis_w_inla_mod$mod$range*amis_w_inla_mod$scale[1]


p1 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$intercept), aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = as.data.frame(mcmc_w_inla_mod$margs$intercept), aes(x=x,y=y,color="MCMC with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$intercept),aes(x=x,y=y,color="IS with INLA")) + 
  labs(y = "",x="a",color = "") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p1
p2 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$prec.range), aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$prec.range), aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data = as.data.frame(mcmc_w_inla_mod$margs$prec.range), aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(y = "",x=expression(tau)) + 
  theme_bw() + 
  coord_cartesian(xlim=c(0,0.5))
p2


p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[1]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]],aes(x=x,y=y,color="MCMC with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]],aes(x=x,y=y,color="IS with INLA")) +  
  labs(y="",x = "c",color = "") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[2]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]],aes(x=x,y=y,color="MCMC with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]],aes(x=x,y=y,color="IS with INLA")) +  
  labs(y= "",x = "d",color = "") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p4
eta_joint_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 200, lims = c(11,16,-0.018,-0.009))  
amis_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))

eta_joint_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 200, lims = c(11,16,-0.018,-0.009))  
is_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_is$x, y=eta_joint_kern_is$y), z=as.vector(eta_joint_kern_is$z))

eta_joint_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1/length(mcmc_w_inla_mod$eta[,1]),length(mcmc_w_inla_mod$eta[,1])), n = 200, lims = c(11,16,-0.018,-0.009))  
mcmc_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_mcmc$x, y=eta_joint_kern_mcmc$y), z=as.vector(eta_joint_kern_mcmc$z))

p5 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z,color = "AMIS with INLA"),bins = 6) + 
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,color = "MCMC with INLA"),bins = 6) + 
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,color = "IS with INLA"),bins = 6) + 
  theme_bw() +
  coord_cartesian(xlim = c(12.5,15.0), y = c(-0.0165,-0.012)) + 
  scale_color_manual(values = col_temp) + 
  labs(x="c",y="d",color ="")
p5
p6 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,color = "AMIS with INLA"),bins = 6) + 
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z,color = "MCMC with INLA"),bins = 6) + 
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,color = "IS with INLA"),bins = 6) + 
  theme_bw() +
  coord_cartesian(xlim = c(12.5,15.0), y = c(-0.0165,-0.012)) + 
  scale_color_manual(values = col_temp) + 
  labs(x="c",y="d",color ="")
p6
p7 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,color = "AMIS with INLA"),bins = 6) + 
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,color = "MCMC with INLA"),bins = 6) + 
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z,color = "IS with INLA"),bins = 6) + 
  theme_bw() +
  coord_cartesian(xlim = c(12.5,15.0), y = c(-0.0165,-0.012)) + 
  scale_color_manual(values = col_temp) + 
  labs(x="c",y="d",color ="")
p7

ggarrange(p5,p6,p7,ncol=1, nrow=3, common.legend = T,legend="bottom")

formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
res_amis = inla(formula,
           data = amis_w_inla_mod$mod,
           control.predictor = list(compute = T),
           scale = exp(amis_w_inla_mod$params[1]+ amis_w_inla_mod$params[2]*amis_w_inla_mod$mod$range),
           control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
           verbose = FALSE)
res_is = inla(formula,
                data = is_w_inla_mod$mod,
                control.predictor = list(compute = T),
                scale = exp(is_w_inla_mod$params[1]+ is_w_inla_mod$params[2]*is_w_inla_mod$mod$range),
                control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
                verbose = FALSE)
res_mcmc = inla(formula,
                data = amis_w_inla_mod$mod,
                control.predictor = list(compute = T),
                scale = exp(mcmc_w_inla_mod$params[1]+ mcmc_w_inla_mod$params[2]*amis_w_inla_mod$mod$range),
                control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
                verbose = FALSE)

tmp = data.frame(x = NA, y=NA)
for (i in seq(6,nrow(amis_w_inla_mod$mod), 5)){
  tmp = rbind(tmp,c(mean(amis_w_inla_mod$mod$range[(i-5):i]),log(1/var(amis_w_inla_mod$mod$logratio[(i-5):i]))))
}
pqr_inla_rw2 <- function(data,mod,params,type){
  quants = c(0.1,0.25,0.5,0.75,0.9)
  mu = mod$summary.linear.predictor[,1]
  tau = exp(params[1] + params[2]*data)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (i in seq(length(quants))){
    if (type == "gaussian"){
      tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
    }else if (type == "gamma"){
      tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau 
    }
    res = rbind(res,data.frame(x = data, y = tmpy, quants = rep(toString(quants[i]),length(data))))
  }
  return(res[-1,])
}

pqr_inla_rw2_more <- function(data,mod,params,type){
  quants = c(0.025,0.05,0.15,0.2,0.3,0.35,0.4,0.45,0.55,0.6,0.65,0.7,0.8,0.85,0.95,0.975)
  mu = mod$summary.linear.predictor[,1]
  tau = exp(params[1] + params[2]*data)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (i in seq(length(quants))){
    if (type == "gaussian"){
      tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
    }else if (type == "gamma"){
      tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau 
    }
    res = rbind(res,data.frame(x = data, y = tmpy, quants = rep(toString(quants[i]),length(data))))
  }
  return(res[-1,])
}

amis_w_inla_mod$pqr = pqr_inla_rw2(amis_w_inla_mod$mod$range,res_amis,amis_w_inla_mod$params,type="gaussian")
is_w_inla_mod$pqr = pqr_inla_rw2(amis_w_inla_mod$mod$range,res_is,is_w_inla_mod$params,type="gaussian")
mcmc_w_inla_mod$pqr = pqr_inla_rw2(amis_w_inla_mod$mod$range,res_mcmc,mcmc_w_inla_mod$params,type="gaussian")
amis_w_inla_mod$pqr_more = pqr_inla_rw2_more(amis_w_inla_mod$mod$range,res_amis,amis_w_inla_mod$params,type="gaussian")
is_w_inla_mod$pqr_more = pqr_inla_rw2_more(amis_w_inla_mod$mod$range,res_is,is_w_inla_mod$params,type="gaussian")
mcmc_w_inla_mod$pqr_more = pqr_inla_rw2_more(amis_w_inla_mod$mod$range,res_mcmc,mcmc_w_inla_mod$params,type="gaussian")
pqr_text = data.frame(x = c(),y = c(),text = c())
for (x in unique(amis_w_inla_mod$pqr$quants)){
  tmpx = amis_w_inla_mod$pqr$x[amis_w_inla_mod$pqr$quants==x] + 10
  tmpy = amis_w_inla_mod$pqr$y[amis_w_inla_mod$pqr$quants==x]
  tmptext = x
  pqr_text = rbind(pqr_text,data.frame(x=tmpx[length(tmpx)],y=tmpy[length(tmpy)],text = x))
}
p6 <- ggplot()+
  geom_line(data = amis_w_inla_mod$pqr, aes(x=x,y=y,linetype = quants)) + 
  geom_line(data = amis_w_inla_mod$pqr_more, aes(x=x,y=y,color = quants)) + 
  geom_point(data = amis_w_inla_mod$mod, aes(x = range, y = logratio),alpha = 0.4)+
  labs(x="range",y="logratio",color = "", linetype = "") + 
  scale_color_manual(values = rep("grey",18)) + 
  scale_linetype_manual(values = rep("solid",5)) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label= text)) + 
  guides(color=F,linetype=F) + 
  theme_bw() 
p6
p7 <- ggplot()+
  geom_line(data = mcmc_w_inla_mod$pqr, aes(x=x,y=y,linetype = quants)) + 
  geom_line(data = mcmc_w_inla_mod$pqr_more, aes(x=x,y=y,color = quants)) + 
  geom_point(data = amis_w_inla_mod$mod, aes(x = range, y = logratio),alpha = 0.4)+
  labs(x="range",y="logratio",color = "", linetype = "") + 
  scale_color_manual(values = rep("grey",18)) + 
  scale_linetype_manual(values = rep("solid",5)) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label= text)) + 
  guides(color=F,linetype=F) + 
  theme_bw() 
p7
p7 <- ggplot()+
  geom_line(data = is_w_inla_mod$pqr, aes(x=x,y=y,linetype = quants)) + 
  geom_line(data = is_w_inla_mod$pqr_more, aes(x=x,y=y,color = quants)) + 
  geom_point(data = is_w_inla_mod$mod, aes(x = range, y = logratio),alpha = 0.4)+
  labs(x="range",y="logratio",color = "", linetype = "") + 
  scale_color_manual(values = rep("grey",18)) + 
  scale_linetype_manual(values = rep("solid",5)) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label= text)) + 
  guides(color=F,linetype=F) + 
  theme_bw() 
p7

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight, norm = T)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight, norm = T)
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta, mcmc_w_inla_mod$times)
essp <- ggplot() + 
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,color = "AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,color = "MCMC with INLA")) +
  scale_x_continuous(labels = c("0 min", "5 min", "10 min", "15 min"),breaks=c(0,60*5,60*10,60*15)) + 
  labs(color = "",x="Runtime",y="ESS",color ="") +
  theme_bw() + 
  coord_cartesian(xlim = c(10,60*16)) + 
  scale_color_manual(values = col_temp) + 
  theme(legend.position="bottom")
essp

calc.post.sd(mcmc_w_inla_mod$eta)
ggarrange(p1,p3,p4,essp,ncol=2, nrow=2, common.legend = T,legend="bottom")

# ImmunogG
load("./sims/pqr/pqr-gamma-ImmunogG-amis-w-inla.Rdata")
load("./sims/pqr-gamma-ImmunogG-mcmc-w-inla.Rdata")
load("./sims/pqr-gamma-ImmunogG-is-w-inla.Rdata")
#load("./sims/pqr/pqr-gamma-ImmunogG-is-w-inla.Rdata")
amis_w_inla_mod$pqr = pqr_inla(amis_w_inla_mod$mod,amis_w_inla_mod$margs,amis_w_inla_mod$theta$a.mu[28,],type = "gamma",domain = c(0,7)) 
is_w_inla_mod$params = calc.post.mean(is_w_inla_mod$eta,is_w_inla_mod$weight)
is_w_inla_mod$pqr = pqr_inla(is_w_inla_mod$mod,is_w_inla_mod$margs,is_w_inla_mod$params,type = "gamma",domain = c(0,7)) 
mcmc_w_inla_mod$pqr = pqr_inla(mcmc_w_inla_mod$mod,mcmc_w_inla_mod$margs,colMeans(mcmc_w_inla_mod$eta),type = "gamma",domain = c(0,7)) 
amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight, norm = T)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight, norm = T)
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta, mcmc_w_inla_mod$times)
is_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})

p1 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$a), aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = as.data.frame(mcmc_w_inla_mod$margs$a), aes(x=x,y=y, color = "MCMC with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$a),aes(x=x,y=y,color="IS with INLA")) + 
  labs(color = "",x="a",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p1

p2 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$b), aes(x=x,y=y, color = "AMIS with INLA")) +
  geom_line(data = as.data.frame(mcmc_w_inla_mod$margs$b), aes(x=x,y=y, color = "MCMC with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$b),aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$b) + 
  labs(color = "",x="b",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p2
p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]],aes(x=x,y=y,color="MCMC with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]],aes(x=x,y=y,color="IS with INLA")) + 
  labs(color = "",x = "c",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]],aes(x=x,y=y,color="MCMC with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]],aes(x=x,y=y,color="IS with INLA")) + 
  labs(color = "",x = "d",y ="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p4

ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = T,legend="bottom")
calc.post.mean(mcmc_w_inla_mod$eta)
calc.post.sd(mcmc_w_inla_mod$eta)
inla.zmarginal(mcmc_w_inla_mod$margs$a)
inla.zmarginal(mcmc_w_inla_mod$margs$b)
calc.post.mean(is_w_inla_mod$eta,is_w_inla_mod$weight)
calc.post.sd(is_w_inla_mod$eta,is_w_inla_mod$weight)
inla.zmarginal(is_w_inla_mod$margs$a)
inla.zmarginal(is_w_inla_mod$margs$b)
calc.post.mean(amis_w_inla_mod$eta,amis_w_inla_mod$weight)
calc.post.sd(amis_w_inla_mod$eta,amis_w_inla_mod$weight)

pqr_inla_more <- function(data,margs,a.mu,type,domain=NA){
  tmpx = data$x
  quants = c(0.025,0.05,0.15,0.2,0.3,0.35,0.4,0.45,0.55,0.6,0.65,0.7,0.8,0.85,0.95,0.975)
  params = c(inla.zmarginal(margs$a,silent=T)[[1]],
             inla.zmarginal(margs$b,silent=T)[[1]],
             a.mu[1],
             a.mu[2])
  if (!anyNA(domain)){
    tmpx = seq(domain[1],domain[2],length.out = 500)
  }
  mu = params[1] + params[2]*tmpx
  tau = exp(params[3] + params[4]*tmpx)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (i in seq(length(quants))){
    if (type == "gaussian"){
      tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
    }else if (type == "gamma"){
      tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau 
    }
    res = rbind(res,data.frame(x = tmpx, y = tmpy, quants = rep(toString(quants[i]),length(tmpx))))
  }
  return(res[-1,])
}
amis_w_inla_mod$pqr = pqr_inla(amis_w_inla_mod$mod,amis_w_inla_mod$margs,amis_w_inla_mod$theta$a.mu[28,],type = "gamma",domain = c(0,7)) 
amis_w_inla_mod$pqr_more = pqr_inla_more(amis_w_inla_mod$mod,amis_w_inla_mod$margs,amis_w_inla_mod$theta$a.mu[28,],type="gamma",domain=c(0,7))
mcmc_w_inla_mod$pqr = pqr_inla(mcmc_w_inla_mod$mod,mcmc_w_inla_mod$margs,c(mean(mcmc_w_inla_mod$eta[,1]),mean(mcmc_w_inla_mod$eta[,2])),type="gamma",domain=c(0,7))
mcmc_w_inla_mod$pqr_more = pqr_inla_more(mcmc_w_inla_mod$mod,mcmc_w_inla_mod$margs,c(mean(mcmc_w_inla_mod$eta[,1]),mean(mcmc_w_inla_mod$eta[,2])),type="gamma",domain=c(0,7))
pqr_text = data.frame(x = c(),y = c(),text = c())
for (x in unique(amis_w_inla_mod$pqr$quants)){
  tmpx = amis_w_inla_mod$pqr$x[amis_w_inla_mod$pqr$quants==x] + 0.2
  tmpy = amis_w_inla_mod$pqr$y[amis_w_inla_mod$pqr$quants==x]
  tmptext = x
  pqr_text = rbind(pqr_text,data.frame(x=tmpx[length(tmpx)],y=tmpy[length(tmpy)],text = x))
}

p5 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants)) + 
  geom_line(data= amis_w_inla_mod$pqr_more,aes(x=x,y=y,color=quants)) + 
  #geom_line(data= mcmc_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants)) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod),aes(x=x,y=y),alpha=0.3) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label= text)) + 
  scale_color_manual(values = rep("grey",18)) + 
  scale_linetype_manual(values = rep("solid",5)) + 
  labs(color="",linetype="",x="Age (years)",y="IgG (g/L)") +
  guides(linetype=F,color = F) + 
  theme_bw()
p5

eta_joint_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 200, lims = c(1,2.5,-0.2,0.35))  
amis_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))

p6 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z),bins = 6,color = "black") + 
  theme_bw() +
  labs(x="c",y="d")
p6

essp <- ggplot() + 
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess)) +
  #geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,color = "IS with INLA")) +
  scale_x_continuous(labels = c("0 min", "5 min", "10 min", "15 min"),breaks=c(0,60*5,60*10,60*15)) + 
  labs(color = "",x="Runtime",y="ESS") +
  theme_bw() + 
  coord_cartesian(xlim = c(0,60*15)) + 
  scale_color_manual(values = col_temp[c(1,3)])
essp

ggarrange(p1,p2,p3,p4,p5,essp,ncol=2, nrow=3, common.legend = T,legend="bottom")


library(LaplacesDemon)
T_s = c(1,2,5,10,15,20,28)
amis_adaptive = lapply(seq(7), function(x){
  tmp = rst(500,mu = amis_w_inla_mod$theta$a.mu[T_s[x],1],sigma=sqrt(amis_w_inla_mod$theta$a.cov[1,1,T_s[x]]),nu=3)
  tmp = sort(tmp)
  y = dst(tmp,mu = amis_w_inla_mod$theta$a.mu[T_s[x],1],sigma=sqrt(amis_w_inla_mod$theta$a.cov[1,1,T_s[x]]),nu=3)
  y = y/max(y) + x
  data.frame(x=tmp,y=y)
})
amis_adaptive2 = lapply(seq(7), function(x){
  tmp = amis_w_inla_mod$eta_kern[[1]]$x
  y = amis_w_inla_mod$eta_kern[[1]]$y
  y = y/max(y) +x
  data.frame(x=tmp,y=y)
})

p1amis <- ggplot() + 
  geom_polygon(data = amis_adaptive2[[1]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[2]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[3]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[4]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[5]],aes(x=y,y=x,fill = "target")) +
  geom_polygon(data = amis_adaptive2[[6]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[7]],aes(x=y,y=x,fill = "target")) + 
  geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) + 
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  scale_x_continuous(label = T_s, breaks = seq(7)) + 
  labs(x = "T",y = expression(beta[1]),color="",fill="") +
  theme_bw() + 
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))+
  coord_cartesian(ylim = c(0,15))
p1amis


amis_adaptive21 = lapply(seq(7), function(x){
  tmp = rst(500,mu = amis_w_inla_mod$theta$a.mu[T_s[x],2],sigma=sqrt(amis_w_inla_mod$theta$a.cov[2,2,T_s[x]])/amis_w_inla_mod$scale[1],nu=3)
  tmp = sort(tmp)
  y = dst(tmp,mu = amis_w_inla_mod$theta$a.mu[T_s[x],2],sigma=sqrt(amis_w_inla_mod$theta$a.cov[2,2,T_s[x]])/amis_w_inla_mod$scale[1],nu=3)
  y = y/max(y) + x
  data.frame(x=tmp,y=y)
})
amis_adaptive22 = lapply(seq(7), function(x){
  tmp = amis_w_inla_mod$eta_kern[[2]]$x
  y = amis_w_inla_mod$eta_kern[[2]]$y
  y = y/max(y) +x
  data.frame(x=tmp,y=y)
})

p2amis <-ggplot() + 
  geom_polygon(data = amis_adaptive22[[1]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[2]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[3]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[4]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[5]],aes(x=y,y=x,fill = "target")) +
  geom_polygon(data = amis_adaptive22[[6]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[7]],aes(x=y,y=x,fill = "target")) + 
  geom_path(data = amis_adaptive21[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[2]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[3]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[4]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[5]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[6]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[7]],aes(x=y,y=x, color = "proposal")) +
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  scale_x_continuous(label = T_s, breaks = seq(7)) + 
  labs(x = "T",y = expression(beta[2]),color="",fill="") +
  theme_bw() + 
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01)) 
p2amis

ggarrange(p1amis,p2amis,common.legend = T,legend="bottom")

ggplot() + 
  geom_path(data= data.frame(x =seq(nrow(mcmc_w_inla_mod$eta)), y=mcmc_w_inla_mod$eta[,1]),aes(x=x,y=y)) 
