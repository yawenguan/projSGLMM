library(devtools)
install_github("yawenguan/projSGLMM")
library(projSGLMM)
library(fields)

# simulate poisson
set.seed(2021)
effrange = 0.2 ###EFF###
family = "poisson"

# use the fields package to compute range from effective range
rhofromfields <- Matern.cor.to.range(effrange, nu=nu, cor.target=.05)
rho <- rhofromfields*sqrt(2*nu)
n = 300
grid = as.matrix(expand.grid(seq(0,1,l=10),seq(0,1,l=10)))
X = coords = rbind(matrix(runif(n*2),ncol=2), grid)
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