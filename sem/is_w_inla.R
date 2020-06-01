require(INLA)
require(spdep)
require(INLABMA)
require(parallel)
require(mvtnorm)

dq.rho.lambda <- function(y, x, sigma = init$cov, log =TRUE) {
  dmvt(y,sigma = sigma, df=3, delta = x, type = "shifted",log=log)
  #dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.rho.lambda <- function(x, sigma = init$cov) {
  as.vector(rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted"))
  #as.vector(rmvnorm(1, mean = x, sigma = sigma))
}

par.is <- function(x, data, theta, t, prior, d.prop, r.prop, fit.inla){
  INLA_crash = T
  while(INLA_crash){
    tryCatch({
      eta = r.prop(theta$a.mu[t,], sigma = theta$a.cov[,,t])
      mod = fit.inla(data = data ,eta = eta)
      INLA_crash = F 
    },error=function(e){
    },finally={})
  }
  weight = mod$mlik + prior(eta) - d.prop(y = eta, x = theta$a.mu[t,], sigma = theta$a.cov[,,t])
  return(list(mlik = mod$mlik, dists = mod$dists, eta = eta, weight = weight, times = Sys.time()))
}


is.w.inla <- function(data, init, prior, d.prop, r.prop, fit.inla, N_0 = 200, N = 400){
  if (detectCores()>10){
    ncores = 10
  }else{
    ncores = detectCores()
  }
  eta = matrix(NA, ncol = length(init$mu), nrow = N_0)
  weight = numeric(N_0)
  theta = list(a.mu = matrix(NA, ncol = length(init$mu), nrow = 2),
               a.cov = array(NA, dim = c(length(init$mu), length(init$mu), 2)))
  theta$a.mu[1,] = init$mu
  theta$a.cov[,,1] = init$cov
  starttime = Sys.time()
  times = numeric(N)
  pb <- txtProgressBar(min = 0, max = N+N_0, style = 3)
  is.list = mclapply(seq(N_0), function(x){
    par.is(x, data, theta, 1, prior, d.prop,r.prop, fit.inla)
  }, mc.set.seed = TRUE, mc.cores = ncores)
  for (i in seq(length(is.list))){
    setTxtProgressBar(pb, i)
    eta[i,]= is.list[[i]]$eta
    weight[i] = is.list[[i]]$weight
  }
  theta = calc.theta(theta,weight,eta,length(eta),2)
  margs = NA
  eta = matrix(NA, ncol = length(init$mu), nrow = N_0)
  weight = numeric(N)
  mlik = numeric(N)
  is.list = mclapply(seq(N), function(x){
    par.is(x, data, theta, 2, prior, d.prop,r.prop, fit.inla)
  }, mc.set.seed = TRUE, mc.cores = ncores)
  for (i in seq(length(is.list))){
    setTxtProgressBar(pb, i+N_0)
    eta[i,] = is.list[[i]]$eta
    margs = store.post(is.list[[i]]$dists,margs,i,N)
    weight[i] = is.list[[i]]$weight
    mlik[i] = is.list[[i]]$mlik
    times[i] = as.numeric(difftime(is.list[[i]]$times,starttime,units = "secs"))
  }
  weight = exp(weight - max(weight))
  return(list(eta = eta,
              theta = theta,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight,
              times = times))
}