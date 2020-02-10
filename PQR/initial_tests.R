library(parallel)


DXX <- data.frame(mod = c("D31","D32","D33","D34",
                          "D41","D42","D43","D44",
                          "D51","D52","D53","D54",
                          "D61","D62","D63","D64"),
                  a = c(-1.6,4.5,-1.0,0.5,
                        1.5,5.0,0.2,-1.0,
                        6.0,3.0,1.0,-1.5,
                        0.5,-2.0,2.0,3.0),
                  b = c(0,0,0,0,
                        0.5,-0.1,7.0,-2.0,
                        -0.4,0.5,-0.1,-6.0,
                        5.00,-0.75,-1.50,0.10),
                  c = c(-0.693,1.030,-1.050,0.000,
                        -0.299,0.000,0.405,0.693,
                        -2.5,-0.1,-1.5,0.5,
                        1.00,-0.50,0.25,-5.00),
                  d = c(0,0,0,0,
                        0,0,0,0,
                        0.8,-0.5,-2.0,1.5,
                        -0.85,-4.00,0.10,1.50),
                  f = c(-0.693,0.693,1.609,1.099,
                        -1.609,0.405,1.609,1.386,
                        0.405,1.099,-0.288,1.386,
                        4.0,-0.2,-3.0,2.0),
                  g = c(0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        1.5,-1.0,5.0,-4.0))

GG_model <- function(mod, n = 200){
  x = runif(n)
  params = DXX[DXX$mod==mod,-1]
  mu = params[[1]] + params[[2]]*x
  sigma = exp(params[[3]] + params[[4]]*x)
  k = exp(params[[5]] + params[[6]]*x)
  theta = exp(mu)/(k^(sigma*sqrt(k)))
  beta = 1/(sigma*sqrt(k))
  y = (rgamma(n = n, shape = k, scale = theta^beta))^(1/beta)
  return(list(
    data = data.frame(x = x, y = y),
    params = params,
    mod = data.frame(mu = mu, sigma = sigma, k = k, 
                     theta = theta,beta=beta)
  ))
}

qGG <- function(params,quant,x){
  if (length(params)==3){
    params = c(params[1],0,params[2],0,params[3],0)
  }else if (length(params==4)){
    params = c(params[1:3],0,params[4],0)
  }else if (length(params)==5){
    params = c(params,0)
  }
  exp(params[[1]] + params[[2]]*x)*
    (qgamma(quant,scale = 1, shape = exp(params[[5]]+params[[6]]*x))/
       exp(params[[5]]+params[[6]]*x))^(exp(params[[3]] + params[[4]]*x)/
                                          sqrt(exp(params[[5]] + params[[6]]*x)))
}

optim_quant <- function(quant,data,init){
  res = rep(NA,nrow(data))
  for (i in seq(length(res))){
    fr = function(x){
      qGG(params = x,quant=quant, x = data$x[i])
    }
    res[i] = optim(init,fr)
  }
  return(res)
}

PQR <- function(mod,n = 200){
  res = GG_model(mod,n)
  quants = c(0.1,0.25,0.5,0.75,0.9)
  init = rep(0,sum(res$params!=0))
  for (quant in quants[1]){
    res2 = optim_quant(quant,res$data,init = init)
  }
  # quant_optim = mclapply(quant,function(x){
  #   optim_quant(quant,data,init)
  # }, mc.set.seed = TRUE, mc.cores = ncores)
  return(res2)
}
