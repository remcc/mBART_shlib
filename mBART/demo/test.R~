## load mbart dynamic library and monbart.R and rpmonbart.R
source("load-monbart.R")

## simulate simple data with one x
set.seed(14) 
sigma = .1
n=200
x = matrix(sort(-1+2*runif(n)),ncol=1)
y = x^3 + sigma*rnorm(n)
x=matrix(x,ncol=1)

## run monotonic bart
set.seed(99)
bfmc = monbart(x,y)

## plot results
plot(x,y)
lines(x,bfmc$yhat.train.mean,col="blue",lwd=3,xlab="x",ylab="posterior mean of f(x)")




