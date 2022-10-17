# projSGLMM
This package contains functions to fit the projection-based SGLMM using Monte Carlo EM and Laplace approximation EM algorithms
The method paper "Fast expectation-maximization algorithms for spatial generalized linear mixed models" can be found here https://arxiv.org/abs/1909.05440

To install our package in R. 

```R
library(devtools)
install_github("yawenguan/projSGLMM")
```

An example of using R function sparse.sglmmGP.mcem to run MCEM algorithm is provided below

```R
library(projSGLMM)
library(fields)

# simulate poisson data
set.seed(2021)
family = "poisson"
n = 1000
grid = as.matrix(expand.grid(seq(0,1,l=20),seq(0,1,l=20)))
coords = rbind(matrix(runif(n*2),ncol=2), grid)
X = matrix(rnorm(nrow(coords)*2),nc=2) #simulate from iidN(0,1).no confounding

# library(RandomFields) # RandomFields R package is not available for R 4.2.1
# rfgp = RFsimulate(RMmatern(nu=1.5, scale=0.18),coords) # if RandomFields package is loaded
# r.e  = rfgp$variable1

# form covariance matrix
covfn <- covfndef(nu=1.5) # nu= 0.5, 1.5, 2.5 or 10
K = covfn(rdist(coords),phi=0.1) # phi = 0.1 corresponds to effective range ~ 0.27; as covfn(0.27,phi=0.1)~=0.05 

# simulate data
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

# # initial values (optional)
# g.lm <-glm(Z~X-1,family=family) 
# beta.init = g.lm$coefficients
# w.init<-g.lm$residuals # residuals taken as random effects
# tau.init = var(w.init)
# offset = NULL
# init = list(tau.init = tau.init,beta.init=beta.init,phi.init = 0.3)

# fit proj. SGLMM using mcem to the simulated data
mcem = sparse.sglmmGP.mcem(Z,X,coords,nu=1.5,family="poisson") #no rank is given, so rank will be selected based on BIC
mcem.pred = sparse.sglmmGP.mcem.pred(mcem,X.pred,coords.pred)

# plot results
set.panel(1,2)
re.est  = mcem$wmean.update[mcem$stopiter+1,]
quilt.plot(coords,re.est,zlim = range(dftrain$r.e),main="Estimated r.e")
quilt.plot(coords,dftrain$r.e,zlim= range(dftrain$r.e),main="True r.e")

linear.mean.est  = (re.est+mcem$X%*%mcem$beta.est[mcem$stopiter+1,])
linear.mean     = rowSums(dftrain[,4:6])
quilt.plot(coords,linear.mean.est,zlim = range(linear.mean),main="Estimated xb+r.e")
quilt.plot(coords,linear.mean,zlim= range(linear.mean),main="True xb+r.e")

image.plot(matrix(rowMeans(mcem.pred$pred.re),20),zlim=range(dftest$r.e),main="Predicted r.e")
image.plot(matrix(dftest$r.e,20),zlim=range(dftest$r.e),main="True r.e")

linear.mean = rowSums(dftest[,4:6])
image.plot(matrix(log(rowMeans(mcem.pred$pred.mean)),20),zlim=range(linear.mean),main="Predicted xb+r.e")
image.plot(matrix(linear.mean,20),zlim=range(linear.mean),main="True xb+r.e")
```
