setwd('GPRPanel')
source('GPRPanel.R')
library(rstan)
set.seed(99)
x1=rnorm(15*15)
x2=rnorm(15*15)
y=1+.5*x1+-.95*x2+rnorm(15*15,0,.25)
data=data.frame(y=y, x1=x1, x2=x2, group=rep(1:15, each=15), time=rep(1:15, 15))
tt = GPRPanel(y~x1+x2, 'group', 'time', data)
summary(tt)$summary
tt2 = GPRPanel(y~x1+x2, 'group', 'time', data, loglik=T, iter=500)
summary(tt2)$summary


#### prediction ####

source('GPRPanelPred.R')
set.seed(91)
x1=rnorm(15*15)
x2=rnorm(15*15)
y=1+.5*x1+-.95*x2+rnorm(15*15,0,.25)
data=data.frame(y=y, x1=x1, x2=x2, group=rep(1:15, each=15), time=rep(1:15, 15))
Z = cbind(rnorm(15), rnorm(15), diag(1,15,15))
z_time = 1.788854
Z_corr = cbind(Z, z_time)
tt3 = GPRPanelPred(y~x1+x2, 'group', 'time', data, Z=Z, Z_corr=Z_corr, iter=500, seed=7867)
cor(1+.5*Z[,1]-.95*Z[,2],
    summary(tt3)$summary[paste0('z[',1:15,']'),'mean'])
