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

params = amis_w_inla_mod$theta$a.mu[nrow(amis_w_inla_mod$theta$a.mu),]
res = inla(formula,
           data = lidar,
           control.predictor = list(compute = T),
           scale = exp(params[1]+ params[2]*lidar$range),
           control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
           verbose = FALSE)

tmp = data.frame(x = NA, y=NA)
for (i in seq(6,nrow(lidar), 5)){
  tmp = rbind(tmp,c(mean(lidar$range[(i-5):i]),log(1/var(lidar$logratio[(i-5):i]))))
}
tmp = tmp[-1,]
tmp_res = lm(y~x, data = tmp)
formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
res_ml = inla(formula,
           data = lidar,
           control.predictor = list(compute = T),
           scale = exp(tmp_res$coefficients[[1]] + tmp_res$coefficients[[2]]*lidar$range),
           control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
           verbose = FALSE)

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

amis_w_inla_mod$pqr = pqr_inla_rw2(lidar$range,res,params,type="gaussian")
amis_w_inla_mod$pqr_ml = pqr_inla_rw2(lidar$range,res_ml,tmp_res$coefficients,type="gaussian")
fig_pqr <- ggplot()+
  geom_line(data = amis_w_inla_mod$pqr_ml, aes(x=x,y=y,color = quants, linetype = "ML")) + 
  geom_line(data = amis_w_inla_mod$pqr, aes(x=x,y=y,color = quants, linetype = "INLA")) + 
  geom_point(data = amis_w_inla_mod$data, aes(x = range, y = logratio))+
  labs(title="Lidar",x="range",y="logratio",color = "Quantiles", linetype = "Method") + 
  theme_bw() 
fig_pqr

