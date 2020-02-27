require(INLA)
require(INLABMA)
require(coda)
require(spdep)
require(mvtnorm)
require(MASS)

prior.param <- function(x, log = TRUE) {
  sum(dunif(x, -100, 100, log = log))
}


fit.inla <- function(data,eta){
  res = inla(y ~ 1 + x, 
             data = data,
             scale = exp(eta[1] + eta[2]*data$x), 
             family = "gamma",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE,
             quantiles=c(0.1, 0.25, 0.5, 0.75, 0.9))
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           beta = res$marginals.fixed[[2]]),
              quants = list(a = res$summary.fixed[1,3:7],
                            b = res$summary.fixed[2,3:7])))
}


calc.theta <- function(theta,weight,eta,i_tot,i_cur){
  weight[1:i_tot] = exp(weight[1:i_tot] - max(weight[1:i_tot]))
  for (i in seq(ncol(eta))){
    theta$a.mu[i_cur,i] = sum(eta[1:i_tot,i]*weight[1:i_tot])/sum(weight[1:i_tot])
  }
  for (i in seq(ncol(eta))){
    for (j in seq(i,ncol(eta))){
      theta$a.cov[i,j,i_cur] = theta$a.cov[j,i,i_cur] = sum(weight[1:i_tot]*(eta[1:i_tot,i]-theta$a.mu[i_cur,i])*
                                                              (eta[1:i_tot,j]-theta$a.mu[i_cur,j]))/(sum(weight[1:i_tot]))
    }
  }
  return(theta)
}

calc.quants <- function(quants,weight){
  weight = weight/sum(weight)
  new_quants = quants
  for (i in seq(length(quants))){
    new_quants[[i]] = rep(0,ncol(quants[[i]]))
    names(new_quants[[i]]) = colnames(quants[[i]])
    for (j in seq(nrow(quants[[i]]))){
      new_quants[[i]] = new_quants[[i]] + weight[j]*quants[[i]][j,]
    }
  }
  return(new_quants)
}

calc.stats <- function(stats,weight){
  for (i in seq(length(stats))){
    new.stat = c(0,0)
    new.stat[1] =  sum(stats[[i]]*weight)/sum(weight)
    new.stat[2] = sum((stats[[i]] - new.stat[1])*(stats[[i]] - new.stat[1])*weight)/sum(weight)
    stats[[i]] = new.stat
  }
  return(stats)
}

store.quants <- function(quant,quants,j,n.prop){
  if (anyNA(quants)){
    quants = quant
    for (i in seq(length(quant))){
      quants[[i]] = matrix(NA, nrow = n.prop, ncol = length(quant[[i]]))
      colnames(quants[[i]]) = names(quant[[i]])
      quants[[i]][j,] = as.numeric(quant[[i]])
    }
    return(quants)
  }else{
    for (i in seq(length(quant))){
      quants[[i]][j,] = as.numeric(quant[[i]])
    }
    return(quants)
  }
}

store.stats <- function(stat,stats,j,n.prop){
  if (anyNA(stats)){
    stats = stat
    for (i in seq(length(stat))){
      stats[[i]] = rep(NA, n.prop)
      stats[[i]][j] = stat[[i]]
    }
    return(stats)
  }else{
    for (i in seq(length(stat))){
      stats[[i]][j] = stat[[i]]
    }
    return(stats)
  }
}

store.post <- function(marg,margs,j,n.prop){
  if (anyNA(margs)){
    margs = marg
    for (i in seq(length(marg))){
      margs[[i]] = list(x = matrix(NA, nrow = n.prop, ncol = length(marg[[i]][,1])),
                        y = matrix(NA, nrow = n.prop, ncol = length(marg[[i]][,2])))
      margs[[i]]$x[j,] = marg[[i]][,1]
      margs[[i]]$y[j,] = marg[[i]][,2]
    }
    return(margs)
  }else{
    for (i in seq(length(marg))){
      margs[[i]]$x[j,] = marg[[i]][,1]
      margs[[i]]$y[j,] = marg[[i]][,2]
    }
    return(margs)
  }
}

amis_pqr <- function(eta,weight,quants,data){
  require(spatstat)
  amis_kerns = lapply(seq(ncol(eta)), function(x){
    dens = density(x = eta[,x],
                   weights = weight/sum(weight), 
                   kernel = "gaussian")
    quants = quantile(dens,c(0.1,0.25,0.5,0.75,0.9))
    return(list(dens = as.data.frame(dens[c(1,2)]), quants = quants))
  })
  quants = list(a = quants$a, b = quants$b, f = amis_kerns[[1]]$quants, g = amis_kerns[[2]]$quants)
  pqr = pqr_inla(x = data$x, a = quants$a, b=quants$b, f = quants$f, g = quants$g)
  return(list(eta_kern=list(f = amis_kerns[[1]]$dens,g = amis_kerns[[2]]$dens),pqr = pqr, quants=quants))
}

pqr_inla <- function(x, a, b, f, g){
  quants = c(0.1,0.25,0.5,0.75,0.9)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (i in seq(length(a))){
    tmpquant = exp(a[i] + b[i]*x)*qgamma(quants[i],shape = exp(f[i] + g[i]*x),scale = 1)/exp(f[i] + g[i]*x)
    res = rbind(res,data.frame(x = x, y = tmpquant, quants = rep(toString(quants[i]),length(x))))
  }
  return(res[-1,])
}

running.ESS <- function(eta, times, ws = NA, norm = TRUE,step = 100){
  if (anyNA(ws)){
    require(coda)
    ess = sapply(lapply(seq(2,nrow(eta)),function(x){
      effectiveSize(eta[1:x,])
    }),min)
    times = times[-1]
  }else{
    if (norm){
      ws = ws/sum(ws)
    }
    ess = sapply(seq(length(ws)),function(x){
      sum(ws[1:x])^2/(sum(ws[1:x]^2))
    })
    rm.ess = !is.na(ess)
    times = times[rm.ess]
    ess = ess[rm.ess]
  }
  ess.df = data.frame(time = c(times[1],times[rev(seq(length(times),100,-step))]),
                      ess = c(ess[1],ess[rev(seq(length(ess),100,-step))]))
  return(ess.df)
}

fit.marginals <- function(ws,margs,len = 400){
  ws = ws/sum(ws)
  xmin <- quantile(apply(margs[[1]],1,function(X){min(X)}),0.25)
  xmax <- quantile(apply(margs[[1]],1,function(X){max(X)}),0.75)
  xx <- seq(xmin, xmax, len = len)
  marg = numeric(len)
  for (i in seq(nrow(margs[[1]]))){
    marg = marg + ws[i]*inla.dmarginal(xx, list(x = margs[[1]][i,], y = margs[[2]][i,]))
  }
  data.frame(x = xx, y = marg)
}

kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n) 
  gy <- seq(lims[3], lims[4], length = n) 
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] 
  ay <- outer(gy, y, "-")/h[2]
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) 
  return(list(x = gx, y = gy, z = z))
}

