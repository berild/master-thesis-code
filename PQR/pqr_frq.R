library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)

# DXX <- data.frame(mod = c("D31","D32","D33","D34",
#                           "D41","D42","D43","D44",
#                           "D51","D52","D53","D54",
#                           "D61","D62","D63","D64"),
#                   a = c(-1.6,4.5,-1.0,0.5,
#                         1.5,5.0,0.2,-1.0,
#                         6.0,3.0,1.0,-1.5,
#                         0.5,-2.0,2.0,3.0),
#                   b = c(0,0,0,0,
#                         0.5,-0.1,7.0,-2.0,
#                         -0.4,0.5,-0.1,-6.0,
#                         5.00,-0.75,-1.50,0.10),
#                   c = c(-0.693,1.030,-1.050,0.000,
#                         -0.299,0.000,0.405,0.693,
#                         -2.5,-0.1,-1.5,0.5,
#                         1.00,-0.50,0.25,-5.00),
#                   d = c(0,0,0,0,
#                         0,0,0,0,
#                         0.8,-0.5,-2.0,1.5,
#                         -0.85,-4.00,0.10,1.50),
#                   f = c(-0.693,0.693,1.609,1.099,
#                         -1.609,0.405,1.609,1.386,
#                         0.405,1.099,-0.288,1.386,
#                         4.0,-0.2,-3.0,2.0),
#                   g = c(0,0,0,0,
#                         0,0,0,0,
#                         0,0,0,0,
#                         1.5,-1.0,5.0,-4.0))


fix_params <- function(x){
  if (length(x)==3){
    x = c(x[[1]],0,x[[2]],0,x[[3]],0)
  }else if (length(x)==4){
    x = c(x[[1]],x[[2]],x[[3]],0,x[[4]],0)
  }else if (length(x)==5){
    x = c(x,0)
  }
  return(x)
}



ml_gg <- function(x,covar,res){
  x = fix_params(x)
  mu = x[[1]] + x[[2]]*covar
  sigma = exp(x[[3]] + x[[4]]*covar)
  k = exp(x[[5]] + x[[6]]*covar)
  w = (res-mu)/sigma
  return(sum(-log(sigma) + (k-1/2)*log(k) - lgamma(k) + w*sqrt(k) - k*exp(w/sqrt(k))))
}

ml_optim <- function(covar,response,init){
  fr = function(x){
    ml_gg(x=x, covar=covar,res = log(response))
  }
  tmpres = optim(init,fr,control=list(fnscale=-1))
  return(tmpres$par)
}

q_reg_GG <- function(params,covar,quants){
  params = fix_params(params)
  mu = params[[1]] + params[[2]]*covar
  sigma = exp(params[[3]] + params[[4]]*covar)
  k = exp(params[[5]] + params[[6]]*covar)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (quant in quants){
    tmpr = qgamma(quant,shape = k,scale = 1)
    tmpquant = exp(mu)*(tmpr/k)^(sigma*sqrt(k)) 
    res = rbind(res,data.frame(x = covar, y = tmpquant, quants = rep(toString(quant),length(covar))))
  }
  return(res[-1,])
}


PQR <- function(x,y,init,dom){
  quants = c(0.1,0.25,0.5,0.75,0.9)
  ml_params = ml_optim(x,y,init)
  quantiles = q_reg_GG(ml_params,seq(dom[1],dom[2],length.out = 200),quants)
  return(list(ml = ml_params,
              quantiles = quantiles))
}

  