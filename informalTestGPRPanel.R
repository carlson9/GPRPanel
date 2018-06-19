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
