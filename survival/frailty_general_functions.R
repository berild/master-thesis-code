require(INLA)
require(INLABMA)
require(coda)
require(spdep)
require(mvtnorm)
require(MASS)

prior.frailty <- function(x, log = TRUE) {
  sum(dgamma(x,shape = 1,rate = 0.01,log = log))
}

fit.inla.kidney <- function(data,eta){
  data$oset = log(eta[data$id])
  formula = inla.surv(time,status) ~ age + sex + disease + offset(oset)
  res=inla(formula,
           family ="weibullsurv",
           data=data,
           control.family = list(list(variant = variant)))
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = arginals.fixed[[1]],
                           age = res$marginals.fixed[[2]],
                           sex = res$marginals.fixed[[3]],
                           disease = res$marginals.fixed[[4]])))
}

fit.inla <- function(data,eta){
  data$oset = log(eta[data$idx])
  formula = inla.surv(y,event) ~ x + offset(oset)
  res=inla(formula,
           family ="weibullsurv",
           data=data,
           control.family = list(list(variant = variant)))
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           beta = res$marginals.fixed[[2]],
                           alpha = res$marginals.hyperpar[[1]])))
}

dq.param <- function(param, eta, mlik){
  weight = numeric(length(param))
  for (j in seq(length(param))){
    for (i in seq(nrow(eta))){
      weight[j] = weight[j] + 
        exp(sum(dgamma(eta[i,],rate = param[j],shape = param[j],log=TRUE)) + 
              mlik[i] + prior.frailty(param[j]))
    }
  }
  weight = weight/(sum(weight)*(param[2]-param[1]))
  return(weight)
}

calc.param <- function(mlik,eta,weight){
  frailty.seq = seq(0,5,length.out = 101)[-1]
  frailty.weight = dq.param(frailty.seq,eta,mlik)
  return(data.frame(x=frailty.seq,y = frailty.weight))
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


calc.stats <- function(stats,weight){
  for (i in seq(length(stats))){
    new.stat = c(0,0)
    new.stat[1] =  sum(stats[[i]]*weight)/sum(weight)
    new.stat[2] = sum((stats[[i]] - new.stat[1])*(stats[[i]] - new.stat[1])*weight)/sum(weight)
    stats[[i]] = new.stat
  }
  return(stats)
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

amis_kde <- function(eta,weight){
  if (is.null(ncol(eta))){
    return(as.data.frame(density(x = eta,
                          weights = weight/sum(weight), 
                          kernel = "gaussian")[c(1,2)]))
  }else{
    return(lapply(seq(ncol(eta)), function(x){
      as.data.frame(density(x = eta[,x],
                            weights = weight/sum(weight), 
                            kernel = "gaussian")[c(1,2)])
    }))
  }
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

kde.quantile <- function(kde){
  n_kd = length(kde)
  quants = data.frame(idx = seq(n_kd),q0.025 = numeric(n_kd),q0.5 = numeric(n_kd),q0.975 = numeric(n_kd))
  for (i in seq(n_kd)){
    tmp = 0
    q25 = T
    q5 = T
    step.size = kde[[i]]$x[2]- kde[[i]]$x[1]
    for (j in seq(nrow(kde[[i]]))){
      tmp = tmp + kde[[i]]$y[j]*step.size
      if ((tmp > 0.025)&(q25)){
        quants$q0.025[i] = kde[[i]]$x[j]
        q25 = F
      }else if ((tmp > 0.5)&(q5)){
        quants$q0.5[i] = kde[[i]]$x[j]
        q5 = F
      }else if ((tmp > 0.975)){
        quants$q0.975[i] = kde[[i]]$x[j]
        break
      }
    }
  }
  return(quants)
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

