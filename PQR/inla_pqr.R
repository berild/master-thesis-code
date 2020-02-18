library(parallel)
library(ggplot2)
library(ggpubr)
library(Brq)
library(INLA)

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
                  f = c(1,1,1,1,
                        1,1,1,1,
                        1,1,1,1,
                        1,1,1,1),
                  g = c(0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0))

GG_model <- function(mod, n = 200){
  x = runif(n)
  params = DXX[DXX$mod==mod,-1]
  mu = params[[1]] + params[[2]]*x
  sigma = exp(params[[3]] + params[[4]]*x)
  k = exp(params[[5]] + params[[6]]*x)
  theta = exp(mu)/(k^(sigma*sqrt(k)))
  beta = 1/(sigma*sqrt(k))
  scale = theta^beta
  y = (rgamma(n = n, shape = k, scale = scale))^(1/beta)
  return(list(
    data = data.frame(x = x, y = y),
    params = params,
    mod = data.frame(mu = mu, sigma = sigma, k = k, 
                     theta = theta,beta=beta)
  ))
}



eta = mod$params$a + mod$params$b*mod$data$x
mu = exp(eta)
  prec.scale = mod$params$c + mod$params$d*mod$data$x
  prec.par = 1
  a = prec.par * prec.scale
b = mu / (prec.par * prec.scale)
y = rgamma(n, shape = a, scale = b)
r = inla(y ~ 1 + x, data = mod$data,
         scale = prec.scale, family = "gamma")


mod = GG_model("D51",n=500)
res = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
fig1 <- ggplot(res$quantiles, aes(x=x,y=y,color = quant))+
  geom_line() + 
  labs(title="D51",x="x",y="",color = "quantiles") + 
  theme_bw()
fig1
mod = GG_model("D52",n=500)
res = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
fig2 <- ggplot(res$quantiles, aes(x=x,y=y,color = quant))+
  geom_line() + 
  labs(title="D52",x="x",y="",color = "quantiles") + 
  theme_bw()
fig2
mod = GG_model("D53",n=500)
res = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
fig3 <- ggplot(res$quantiles, aes(x=x,y=y,color = quant))+
  geom_line() + 
  labs(title="D53",x="x",y="",color = "quantiles") + 
  theme_bw()
fig3
mod = GG_model("D54",n=500)
res = PQR(mod$data$x,mod$data$y,init = rep(0,sum(mod$params!=0)),dom = c(0,1))
fig4 <- ggplot(res$quantiles, aes(x=x,y=y,color = quant))+
  geom_line() + 
  labs(title="D54",x="x",y="",color = "quantiles") + 
  theme_bw()
fig4


