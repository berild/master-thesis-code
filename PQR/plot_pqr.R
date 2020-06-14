library(ggplot2)
library(INLA)
source("./PQR/general_functions.R")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_temp = gg_color_hue(4)

p1 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$a), aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$a),aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$a) + 
  labs(color = "",x="a",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p1

p2 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$b), aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$b),aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$b) + 
  labs(color = "",x="b",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p2
p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$c) + 
  labs(color = "",x = "c",y="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$d) + 
  labs(color = "",x = "d",y ="") + 
  scale_color_manual(values = col_temp) + 
  theme_bw()
p4

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
load("./sims/pqr/pqr-m1-gaussian-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,1)
p1 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label=paste("M1 Gaussian \n","a =", amis_w_inla_mod$mod$params$a,"   b =", amis_w_inla_mod$mod$params$b,"   c =",amis_w_inla_mod$mod$params$c,
                              "   d =",amis_w_inla_mod$mod$params$d),x = 0.5, y = 12) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p1

load("./sims/pqr/pqr-m2-gaussian-amis-w-inla.Rdata")
load("./sims/pqr/pqr-m2-gaussian-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,2)
p2 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label=paste("M2 Gaussian \n","a =", amis_w_inla_mod$mod$params$a,"   b =", amis_w_inla_mod$mod$params$b,"   c =",amis_w_inla_mod$mod$params$c,
                              "   d =",amis_w_inla_mod$mod$params$d),x = 0.5, y = 2.5) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p2

load("./sims/pqr/pqr-m1-gamma-amis-w-inla.Rdata")
load("./sims/pqr/pqr-m1-gamma-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,2)
p3 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label=paste("M1 Gamma \n","a =", amis_w_inla_mod$mod$params$a,"   b =", amis_w_inla_mod$mod$params$b,"   c =",amis_w_inla_mod$mod$params$c,
                              "   d =",amis_w_inla_mod$mod$params$d),x = 0.5, y = 80) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p3


load("./sims/pqr/pqr-m2-gamma-amis-w-inla.Rdata")
load("./sims/pqr/pqr-m2-gamma-is-w-inla.Rdata")
pqr_text <- create_pqr_text(amis_w_inla_mod$pqr,1)
p4 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  #geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_line(data = amis_w_inla_mod$pqr_truth,aes(x=x,y=y,linetype = quants,color="Truth")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod$data),aes(x=x,y=y),alpha=0.2) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label=text)) + 
  annotate("text",label=paste("M2 Gamma \n","a =", amis_w_inla_mod$mod$params$a,"   b =", amis_w_inla_mod$mod$params$b,"   c =",amis_w_inla_mod$mod$params$c,
                              "   d =",amis_w_inla_mod$mod$params$d),x = 0.5, y = 0.45) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="",y="") + 
  scale_color_manual(values = col_temp) + 
  guides(linetype = F) + 
  theme_bw()
p4

library(ggpubr)

ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = T,legend="bottom")

# Lidar
load("./sims/pqr/pqr-gaussian-lidar-is-w-inla.Rdata")
load("./sims/pqr/pqr-gaussian-lidar-amis-w-inla.Rdata")


p1 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$intercept), aes(x=x,y=y)) + 
  #geom_line(data = as.data.frame(is_w_inla_mod$margs$intercept),aes(x=x,y=y,color="IS with INLA")) + 
  labs(y = "",x="a") + 
  theme_bw()
p1
p2 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$prec.range), aes(x=x,y=y)) + 
  labs(y = "",x=expression(tau)) + 
  theme_bw() + 
  coord_cartesian(xlim=c(0,0.5))
p2


p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]],aes(x=x,y=y)) + 
  #geom_line(data = is_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="IS with INLA")) +  
  labs(y="",x = "c") + 
  theme_bw()
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]],aes(x=x,y=y)) + 
  #geom_line(data = is_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="IS with INLA")) +  
  labs(y= "",x = "d") + 
  theme_bw()
p4
eta_joint_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 200, lims = c(11,16,-0.018,-0.009))  
amis_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))

p5 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z),bins = 6,color = "black") + 
  theme_bw() +
  labs(x="c",y="d")
p5
# amis_w_inla_mod$theta$a.mu[,2]=amis_w_inla_mod$theta$a.mu[,2]/amis_w_inla_mod$scale[1]
# params = amis_w_inla_mod$theta$a.mu[nrow(amis_w_inla_mod$theta$a.mu),]
# amis_w_inla_mod$eta[,2] = amis_w_inla_mod$eta[,2]/amis_w_inla_mod$scale[1]
# amis_w_inla_mod$eta_kern[[2]]$x = amis_w_inla_mod$eta_kern[[2]]$x/amis_w_inla_mod$scale[1]
# amis_w_inla_mod$mod$range = amis_w_inla_mod$mod$range*amis_w_inla_mod$scale[1]

formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
res = inla(formula,
           data = amis_w_inla_mod$mod,
           control.predictor = list(compute = T),
           scale = exp(params[1]+ params[2]*amis_w_inla_mod$mod$range),
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

amis_w_inla_mod$pqr = pqr_inla_rw2(amis_w_inla_mod$mod$range,res,params,type="gaussian")
amis_w_inla_mod$pqr_more = pqr_inla_rw2_more(amis_w_inla_mod$mod$range,res,params,type="gaussian")
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
  #geom_line(data = is_w_inla_mod$pqr, aes(x=x,y=y,color = quants, linetype = "IS with INLA")) + 
  geom_point(data = amis_w_inla_mod$mod, aes(x = range, y = logratio),alpha = 0.5)+
  labs(x="range",y="logratio",color = "", linetype = "") + 
  scale_color_manual(values = rep("grey",18)) + 
  scale_linetype_manual(values = rep("solid",5)) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label= text)) + 
  guides(color=F,linetype=F) + 
  theme_bw() 
p6

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight, norm = F)
essp <- ggplot() + 
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,color = "AMIS with INLA")) +
  scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) + 
  labs(color = "",x="Runtime",y="ESS") +
  theme_bw() + 
  coord_cartesian(xlim = c(10,60*70)) + 
  scale_color_manual(values = col_temp) + 
  theme(legend.position="bottom")
essp


# ImmunogG
load("./sims/pqr/pqr-gamma-ImmunogG-amis-w-inla.Rdata")
load("./sims/pqr/pqr-gamma-ImmunogG-is-w-inla.Rdata")
amis_w_inla_mod$pqr = pqr_inla(amis_w_inla_mod$mod,amis_w_inla_mod$margs,amis_w_inla_mod$eta_kern,type = "gamma",domain = c(0,7)) 
is_w_inla_mod$pqr = pqr_inla(is_w_inla_mod$mod,is_w_inla_mod$margs,is_w_inla_mod$eta_kern,type = "gamma",domain = c(0,7)) 

p1 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$a), aes(x=x,y=y,color = "AMIS with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$a), aes(x=x,y=y,color = "IS with INLA")) + 
  labs(color = "",x="a",y="") + 
  scale_color_manual(values = col_temp[c(1,3)]) + 
  theme_bw()
p1
p2 <- ggplot() + 
  geom_line(data = as.data.frame(amis_w_inla_mod$margs$b), aes(x=x,y=y,color = "AMIS with INLA")) + 
  geom_line(data = as.data.frame(is_w_inla_mod$margs$b), aes(x=x,y=y,color = "IS with INLA")) + 
  labs(color = "",x="b",y="") + 
  scale_color_manual(values = col_temp[c(1,3)]) + 
  theme_bw()
p2
p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_kern[[1]],aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$c) + 
  labs(color = "",x = "c",y="") +
  scale_color_manual(values = col_temp[c(1,3)]) + 
  theme_bw()
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$eta_kern[[2]],aes(x=x,y=y,color="IS with INLA")) + 
  geom_vline(xintercept = amis_w_inla_mod$mod$params$d) + 
  scale_color_manual(values = col_temp[c(1,3)]) + 
  labs(color = "",x = "d",y="") + 
  theme_bw()
p4
pqr_text = data.frame(x = c(),y = c(),text = c())
for (x in unique(amis_w_inla_mod$pqr$quants)){
  tmpx = amis_w_inla_mod$pqr$x[amis_w_inla_mod$pqr$quants==x]
  tmpy = amis_w_inla_mod$pqr$y[amis_w_inla_mod$pqr$quants==x]
  tmptext = x
  pqr_text = rbind(pqr_text,data.frame(x=tmpx[length(tmpx)],y=tmpy[length(tmpy)],text = x))
}

p5 <- ggplot() + 
  geom_line(data= amis_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="AMIS with INLA")) + 
  geom_line(data= is_w_inla_mod$pqr,aes(x=x,y=y,linetype=quants,color="IS with INLA")) + 
  geom_point(data = as.data.frame(amis_w_inla_mod$mod),aes(x=x,y=y),alpha=0.3) + 
  geom_text(data = pqr_text, aes(x=x,y=y,label= text)) + 
  scale_linetype_manual(values = c("solid","solid","solid","solid","solid")) + 
  labs(color="",linetype="",x="Age",y="immunogG") + 
  scale_color_manual(values = col_temp[c(1,3)]) + 
  guides(linetype=F) + 
  theme_bw()
p5

essp <- ggplot() + 
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,color = "AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,color = "IS with INLA")) +
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
  labs(x = "T",y = expression(beta[1]),color="",fill="",title="AMIS with INLA") +
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
  labs(x = "T",y = expression(beta[2]),color="",fill="",title="AMIS with INLA") +
  theme_bw() + 
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01)) 
p2amis

ggarrange(p1amis,p2amis,common.legend = T,legend="bottom")