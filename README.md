# projSGLMM
This package contains functions to fit the projection-based SGLMM using Monte Carlo EM and Laplace approximation EM algorithms
The method paper "Fast expectation-maximization algorithms for spatial generalized linear mixed models" can be found here https://arxiv.org/abs/1909.05440

library(devtools)
install_github("yawenguan/projSGLMM")
library(projSGLMM)
library(fields)

# simulate poisson
set.seed(2021)
family = "poisson"
n = 1000
grid = as.matrix(expand.grid(seq(0,1,l=10),seq(0,1,l=10)))
coords = rbind(matrix(runif(n*2),ncol=2), grid)
X = matrix(rnorm(nrow(coords)*2),nc=2)
K = exp(-rdist(coords)/0.2)
r.e = t(chol(K))%*%rnorm(nrow(K))
expo = exp(X[,1]+X[,2]+r.e)
Z = rpois(length(expo),expo)
df = data.frame(s.x=coords[,1],s.y=coords[,2],Z=Z,r.e = r.e,X=X)
dftrain = df[1:n,]
dftest = df[-c(1:n),]

Z = dftrain$Z 
X = as.matrix(dftrain[,-c(1:4)])
coords = cbind(dftrain$s.x,dftrain$s.y) 
X.pred = as.matrix(dftest[,-c(1:4)])
coords.pred = cbind(dftest$s.x,dftest$s.y)
MCN=40 # number of EM iteration

# fit proj. SGLMM using mcem
mcem = sparse.sglmmGP.mcem(Z,X,coords,nu=1.5)
mcem.pred = sparse.sglmmGP.mcem.pred(mcem,X.pred,coords.pred)
quilt.plot(coords.pred,rowMeans(mcem.pred$pred.mean))
