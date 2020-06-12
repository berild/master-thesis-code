library(sp)
library(spdep)
library(INLA)

#Load data
load("./data/dismap_sim_data.RData")
source("./disease/general_functions.R")

#Rename data
spdf <- OyE.sim
rm(OyE.sim)
#Check names
names(spdf)
#[1] "obs.OCavity"   "esp.OCavity"   "obs.Esophagus" "esp.Esophagus"
#[5] "obs.Stomach"   "esp.Stomach"  

#Set shorter names
names(spdf) <- c("Obs.Cb", "Esp.Cb", "Obs.Eso", "Esp.Eso", 
                 "Obs.Est", "Esp.Est","SMR.Cb","SMR.Eso","SMR.Est")

spplot(spdf["SMR.Cb"],col.regions= hcl.colors(20,palette = hcl.pals("sequential")[20], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
spplot(spdf["SMR.Eso"],col.regions= hcl.colors(20,palette = hcl.pals("sequential")[20], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
spplot(spdf["SMR.Est"],col.regions= hcl.colors(40,palette = hcl.pals("sequential")[20], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
#Create dataset for INLA (n x 3)
n <- nrow(spdf)
d <- list(OBS = matrix(NA, nrow = n*3, ncol = 3))
#Add observed
d$OBS[1:n, 1] <- spdf$Obs.Cb#Bucal cancer
d$OBS[n + 1:n, 2] <- spdf$Obs.Eso#Esophagous cancer
d$OBS[2*n + 1:n, 3] <- spdf$Obs.Est#Stomach cancer

#Expected cases
d$EXP <- c(spdf$Esp.Cb, spdf$Esp.Eso, spdf$Esp.Est)

#Area ID
d$AREAID <- rep(1:n, 3)

#Adjacency matrix
nb.spain <- poly2nb(spdf)
W <- nb2mat(nb.spain, style = "B")

#Define groups in data
d$r <- rep(1:3, each = n)
d$rf <- as.factor(d$r)


fit.inla <- function(data, m.delta) {
  data$delta <- m.delta[data$r]
  form.joint <- OBS ~ -1 + rf + 
    f(AREAID, model = "besag", graph = W, replicate = r,
      hyper = list(prec = list(param = c(0.01, 0.01)))) + 
    f(AREAID.com, delta, model = "besag", graph = W,
      hyper = list(prec = list(param = c(0.01, 0.01)))) 
  
  res <- inla(form.joint, data = data, family = rep("poisson", 3), E = data$EXP)
  return(list(mlik = res$mlik[[1]],
              dists = list(rf1 = res$marginals.fixed[[1]],
                           rf2 = res$marginals.fixed[[2]],
                           rf3 = res$marginals.fixed[[3]],
                           tau1 = res$marginals.hyperpar[[1]],
                           tau2 = res$marginals.hyperpar[[2]])))
}

dq.delta <- function(x, y, sigma, log = TRUE) {
  #loc = log(x) - 1/2*(1+log(diag(sigma))-2*log(x))
  #scale = sqrt(1 + log(diag(sigma)) - 2*log(x))
  res <- dlnorm(y, meanlog = x, sdlog = sqrt(diag(sigma)), log = log)	
  if(log) {
    return(sum(res)) 
  } else {
    return(prod(res))
  }
}
#Generate random values
rq.delta <- function(x, sigma) {
  #loc = log(x) - 1/2*(1+log(diag(sigma))-2*log(x))
  #scale = sqrt(1 + log(diag(sigma)) - 2*log(x))
  rlnorm(length(x), meanlog = x, sdlog = sqrt(diag(sigma)))
}

#Prior for delta (log-normal)
prior.delta <- function(x, mulog = 0, sigmalog = sqrt(1/.1), log = TRUE) {
  res <- dlnorm(x, meanlog = mulog, sdlog = sigmalog, log = log)
  
  if(log) {
    return(sum(res)) 
  } else {
    return(prod(res))
  }
}


#init = list(mu = exp(log(c(0.2522, 0.2753, 0.0558)) + 1/2), cov = exp(2*c(0.2522, 0.2753, 0.0558))*diag(3))
init = list(mu = log(c(0.2522, 0.2753, 0.0558)), cov = 2*diag(3))
# source("./disease/amis_w_inla.R")
# amis_w_inla_mod <- amis.w.inla(data = d, init = init, prior.delta,
#                                dq.delta, rq.delta, fit.inla,
#                                N_t = seq(25,50,1)*10, N_0 = 250)
# save(amis_w_inla_mod, file = "./sims/disease-amis-w-inla.Rdata")

source("./disease/is_w_inla.R")
is_w_inla_mod <- is.w.inla(data = d, init = init, prior.delta,
                           dq.delta, rq.delta, fit.inla, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./sims/disease-is-w-inla.Rdata")
