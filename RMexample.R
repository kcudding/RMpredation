#using Sim.DiffProc 
rm(list=ls())
library("Sim.DiffProc")
library(rootSolve)
# the model : RM predation
sg=c(0, 0.001)
tt=c("RM predation without noise (stable focus)", "RM predation with noise (stable focus)")
par(mfrow=c(1,2))

for (i in seq_along(sg)) {
#parameters all scaled for non-D
alpha = .9 
mu= .1
K = 1 
sigma=sg[i] #noise

#approximate eq'm for starting conditions
s1=mu/(1-mu)
s2=(1+s1)*(1-s1)
s1=s1+rnorm(1,0,0.005)
s2=s2+rnorm(1,0,0.005)

rate <- function(x, a=.9, m=.08) {
  F1= x[1]*(1-(x[1])) - ((x[1]*x[2])/(1+x[1]))
  F2= (a*x[1]*x[2])/(1+x[1]) - m*x[2]
  c(F1 = F1, F2 = F2)
}
  ss <- multiroot(f = rate, start = c(s1, s2))
  
  
  fx <- expression(x*(1-(x)) - ((x*y)/(1+x)),(alpha*x*y)/(1+x) - mu*y)
  gx <- expression(sigma*y,sigma*x)
  
  mod2d2 <- snssde2d(drift=fx,diffusion=gx,M=1,T=1000, N=10000,x0=c(x0=s1,y0=s2),method='rk1')
  mod2d2
  ts.plot(mod2d2$X, mod2d2$Y, main=tt[i], col=c(1,2), lwd=2)


}
  legend("left", c("prey", "predator"), col=c(1,2), lwd=2, bty="n")
 plot2d(mod2d2)
 
 df=data.frame(prey=mod2d2$X, predator=mod2d2$X)

 write.csv(df, "rmexample.csv")
  