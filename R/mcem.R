#' Fit MCEM to SGLMM xxxxx
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param xx Path to the input file
#' @param yy Path to the input file
#' @return Aa matrix of the infile
#' @return Ab matrix of the infile
#' @export
sparse.sglmmGP.mcem <- function(Z,X,family,offset = NULL,q,MCN,size,nu = nu,covfn=covfn,mc.cores=1,coords = coords,rdist=rdist,
                                 mul=2,init=NULL,ceil = ceil, epsilon=1e-3, zalpha = zalpha, ifbreak = ifbreak,RSR=RSR,
                                 RP=RP,tune=0.1,track=T, ALL){
  ptm = proc.time()
  n = length(Z)
  Nbinom = rep(1,n)
  size <- 1e2*q # initial size
  stopiter <- stoptime <- Mstop<-MMstop<-Ustop<-dstop<-samp <- NULL;  flag <- T; #for stopping EM
  
  if(is.null(offset)) offset = rep(1,n)
  p <- ncol(X) # define the number of fixed effects
  PPERP <- diag(n)- X%*%chol2inv(chol(crossprod(X,X)))%*%t(X)
  dist <- rdist(coords,coords)
  #initials for MCEMG --------------
  rhomax <- Matern.cor.to.range(1.5*max(dist), nu=nu, cor.target=.05)
  rhomin <- Matern.cor.to.range(0.1*max(dist), nu=nu, cor.target=.05)
  phigrid = seq(rhomin*sqrt(2*nu),rhomax*sqrt(2*nu),by=0.01)
  
  if(is.null(init$phi.init)){
    # if no initial; starts from effrange = the half maxdist 
    phi.init =  Matern.cor.to.range(0.5*max(dist), nu=nu, cor.target=.05) * sqrt(2*nu)
  }else{phi.init = init$phi.init}
  phiidx = which.min(abs(phigrid-phi.init));phi.init = phigrid[phiidx]
  
  K = covfn(dist,phi.init)
  if(!RP){
    K = eigen(K)
    U = K$vectors[,1:q]
    d = K$values[1:q]
  }else{
    # use random projection
    K = rpK(n,mul*q,K)
    U = K$u[,1:q]
    d = K$d[1:q]
  }
  MM <- U%*%diag(sqrt(d))
  
  g.lm <-glm(Z~X-1 + offset(log(offset)),family=family) 
  
  if(is.null(init$beta.init)){
    beta.init<- g.lm$coefficients
    w.init<-g.lm$residuals # residuals taken as random effects
    tau.init = var(w.init)
    deltalast<-delta.init<-rep(0,q) # t(MM)%*%w.init
    bin1 <- likfit(list(coords=coords,data=w.init), ini = c(tau.init,0.2), fix.nugget = T, fix.kappa = T,kappa=0.5, lik.method = "ML", cov.model="matern") 
    phi.init = bin1$phi
    tau.init = bin1$sigmasq
    phiidx = which.min(abs(phigrid-phi.init));phi.init = phigrid[phiidx]
  }else{
    beta.init<- init$beta.init
    deltalast<- delta.init<-rep(0,q)
    w.init<-rep(0,n)
    tau.init = init$tau.init
    phi.init = init$phi.init
  }
  
  # Define link function  
  if(family=="binomial"){
    # define the link function
    link = function(x){
      ExpZ = 1/(1+exp(-x)) # YG: extend for binomial
      VarZ = ExpZ*(1-ExpZ)
      return(list(ExpZ=ExpZ,VarZ=VarZ))
    }
    
    delta.lf <-function(M,delta,tau,xbeta){
      r.e <- M%*%delta
      eta <- xbeta+r.e
      foo <- -sum(log(1+exp((1-2*Z)*eta))) # Binary likelihood YG: extend to binomial
      TWQW <- crossprod(delta,delta)
      foo <- foo - q/2*log(tau) - 1/(2*tau)*TWQW # foo- 1/(2*tau)*TWQW 
      return(list(ln = foo,TWQW=TWQW))
    }
  }
  
  if(family=="poisson"){
    # define the link function
    link = function(x){
      VarZ = ExpZ = exp(x+log(offset)) # or exp(x)/(1+exp(x))
      return(list(ExpZ=ExpZ,VarZ=VarZ))
    } 
    
    delta.lf <-function(M,delta,tau,xbeta){ # M, xbeta, offset, Z,
      r.e<- M%*%delta
      expo<-exp(xbeta+r.e+log(offset))
      foo<-sum(dpois(Z,lambda = expo,log=T)) # Poisson likelihood
      TWQW <- crossprod(delta,delta)
      foo <- foo - q/2*log(tau) - 1/(2*tau)*TWQW # foo- 1/(2*tau)*TWQW
      #         return(list(ln = foo,TWQW=TWQW,xbeta = xbeta,M=M))
      return(list(ln = foo,TWQW=TWQW))
    }
  }
  
  # matrices to store some output from the algorithm
  tau.update<-matrix(NA,nrow=MCN+1,ncol=1)
  phi.update<-matrix(NA,nrow=MCN+1,ncol=1)
  beta.update<-matrix(NA,nrow=MCN+1,ncol=p)
  deltamean.update<-matrix(NA,nrow=MCN+1,ncol=q) # matrix to store estimates of random effect 
  wmean.update<-matrix(NA,nrow=MCN+1,ncol=n)
  log.list<- vector("list", MCN+1) 
  adaptsize = c()
  tau.update[1,]<-tau.init
  phi.update[1,]<-phi.init
  beta.update[1,]<-beta.init
  deltamean.update[1,]<-0
  wmean.update[1,]<-0
  # store some acceptance rate
  accep = c()
  checksize = setsize = 100*q # minimum MC sample size
  proposal = diag(tune,q)
  
  deltaloop <- function(U,d,tau,beta){
    xbeta <- X%*%beta
    tic = proc.time()
    MMnew <- U%*%diag(sqrt(d))
    signdiag = sign(diag(t(MMnew)%*%MM)) # find the sign change, MM is loaded from global environment
    signdiag = as.logical(1-signdiag)
    MMnew[,signdiag] = -MMnew[,signdiag]
    U[,signdiag] = -U[,signdiag]
    if(RSR) M <- PPERP%*%MMnew else M <- MMnew
    
    delta.update <-matrix(NA,nrow=size,ncol=q)
    TWQW.update  <-log.l <- matrix(NA,nrow=size,ncol=1)
    if(ALL) a.delta=0 else a.delta = rep(0,q)
    # obtain an MC sample of random effects
    deltastar <- delta <- delta.update[1,] <-deltalast
    lrcur <- delta.lf(M,delta,tau,xbeta)
    log.l[1] <- lrcur$ln
    TWQW.update[1] <-lrcur$TWQW
    EZ_approx <- rep(0,n)
    
    for (m in 1:(size-1)){ # number of MC iteration
      # 1.sample delta;2.
      
      if(ALL){# update random effects all-at-once, M-H algo to generate MC samples from f(delta|Y,beta,tau)
        deltastar <- proposal%*%rnorm(q)+delta
        lu <- log(runif(1))
        lrstar <- delta.lf(M,deltastar,tau,xbeta)
        lr <- lrstar$ln - lrcur$ln
        
        if(lu < lr){
          delta   <- deltastar
          a.delta <- a.delta+1 
          lrcur   <- lrstar
        }}else{
          for(re.j in 1:q){
            # update random effects all-at-once, M-H algo to generate MC samples from f(delta|Y,beta,tau)
            deltastar[re.j] <- rnorm(1,delta[re.j], sd= proposal[re.j,re.j])
            lu <- log(runif(1))
            lrstar <- delta.lf(M,deltastar,tau,xbeta)
            lr <- lrstar$ln - lrcur$ln
            
            if(lu < lr){
              delta   <- deltastar
              a.delta[re.j] <- a.delta[re.j]+1 
              lrcur   <- lrstar
            }
          }
        }
      
      delta.update[m+1,] <- delta
      log.l[m+1]         <- lrcur$ln
      TWQW.update[m+1]   <- lrcur$TWQW
      
      if((m+1)%%setsize == 0 & i == 1){ # check ESS for the first EM iter
          if(min(ess(delta.update[1:(m+1),])) >= 2*q){break}
      }
      if((m+1)%%(setsize/2) == 0 & i > 1){ # check adapt MC size
        foo <- link(c(xbeta) + M%*%t(delta.update[(m+2-setsize/2):(m+1),])) # list of two "ExpZ" "VarZ"
        EZ_approx <- (rowMeans(foo[[1]]) * setsize/2 + EZ_approx * (m+1-setsize/2))/(m+1) # weighted mean
        
        if((m+1)>=max(setsize,checksize)){ # allow MCMC to stop when condition met, also make sure the mc size always increases
          if(family =="binomial"){VZ_approx<-EZ_approx*(1-EZ_approx)}else{VZ_approx<-EZ_approx}
          TWQW_approx <- mean(TWQW.update[1:(m+1)])
          
          # update parameter
          betanew<- beta + solve(t(X)%*%(VZ_approx*X),t(X)%*%(Z-EZ_approx))
          taunew <- tau - (q/(2*tau^2)-TWQW_approx/tau^3)^(-1)*(-q/(2*tau) + TWQW_approx/(2*tau^2))
          
          # if(taunew<=0) cat("Update Tau",taunew,"\n")
          KK  = 1
          while(taunew <=0){
            taunew <- tau - 1/(KK*2)*(q/(2*tau^2)-TWQW_approx/tau^3)^(-1)*(-q/(2*tau) + TWQW_approx/(2*tau^2))
            KK = KK+1
          }
          if(KK>1){
            betanew<- beta + 1/(KK*2)*solve(t(X)%*%(VZ_approx*X),t(X)%*%(Z-EZ_approx))
          }
          
          r.e <- M%*%t(delta.update[seq(1,(m+1),by=1e1),])
          EZ <- link(c(X%*%betanew) + r.e)$ExpZ
          if(family =="binomial"){foo<-dbinom(Z,Nbinom,prob = EZ,log=T)}else{foo<-dpois(Z,lambda = EZ,log=T)} # Poisson likelihood
          foo <- colSums(foo)
          lfnew <- foo - q/2*log(taunew) - 1/(2*taunew)*TWQW.update[seq(1,(m+1),by=1e1)] 
          lf0 <- lfnew - log.l[seq(1,(m+1),by=1e1)]
          bmlf<- bm(lf0)
          # cat("taunew",taunew,"\n")
          #           if(!is.na(bmlf$est)&!is.na(bmlf$se)) if((bmlf$est - zalpha*bmlf$se >= 0)) break # if < 0 then simulate more
          # break condition: increase with high probablity or near
          if(bmlf$est - 2*bmlf$se >= 0 | KK > 1) break # if < 0 then simulate more
        }
      }
    }
    
    if(i==1){
      r.e <- M%*%t(delta.update[seq(1,(m+1),len=1e3),])
      foo <- link(c(xbeta) + r.e) # list of two "ExpZ" "VarZ"
      EZ_approx <- rowMeans(foo[[1]]) 
      if(family =="binomial"){VZ_approx<-EZ_approx*(1-EZ_approx)}else{VZ_approx<-EZ_approx}
      TWQW_approx = mean(TWQW.update[seq(1,(m+1),len=1e3)])
      betanew <- beta +solve(t(X)%*%(VZ_approx*X),t(X)%*%(Z-EZ_approx))
      taunew <- 1/q * TWQW_approx  
      
      # if(taunew<=0) cat("Update Tau",taunew,"\n")
      KK  = 1
      while(taunew <=0){taunew <- tau - 1/(KK*2)*(q/(2*tau^2)-TWQW_approx/tau^3)^(-1)*(-q/(2*tau) + TWQW_approx/(2*tau^2));KK = KK+1}
      if(KK>1){
        betanew<- beta + 1/(KK*2)*solve(t(X)%*%(VZ_approx*X),t(X)%*%(Z-EZ_approx))
      }
      
      EZ <- link(c(X%*%betanew) + r.e)$ExpZ
      
      if(family =="binomial"){foo<-dbinom(Z,Nbinom,prob = EZ,log=T)}else{foo<-dpois(Z,lambda = EZ,log=T)} # Poisson likelihood
      foo<-colSums(foo)
      lfnew <- foo - q/2*log(taunew) - 1/(2*taunew)*TWQW.update[seq(1,(m+1),len=3e3)]
      lf0 <- lfnew - log.l[seq(1,(m+1),len=3e3)]
      bmlf <- bm(lf0)
    }
    
    meanlnc = mean(log.l[1:(m+1)])
    tac = proc.time() - tic
    return(list(deltasample = delta.update[1:(m+1),],meanlnc=meanlnc, tac=tac,
                log.l= log.l[1:(m+1)],log.ldiff= lf0, M=M,MM=MMnew,U=U,d=d,
                betanew = betanew,taunew=taunew, deltalast = deltalast,
                accept = a.delta/(m+1),mcsize=(m+1)))
  }
  
  phiupdate <- function(phix,tau,delta.update,Uold,dold){ # xbeta MM tau delta.update can pulled outside from here
    tic = proc.time()
    if(phix<=0|phix>=length(phigrid)) return(list(meanlnc=NA))
    phi = phigrid[phix]
    K = covfn(dist,phi)
    
    if(!RP){
      K = eigen(K)
      Unew = K$vectors[,1:q]
      dnew = K$values[1:q]
    }else{
      K = rpK(n,2*q,K)
      Unew = K$u[,1:q]
      dnew = K$d[1:q]
    }
    
    MM = Uold%*%diag(sqrt(dold))
    # compute the difference in likelihood function due to phi
    # -1/(2sigma^2) w'(UDinvU_new - UDinvU_old)W
    foo1 = Unew%*%((1/dnew)*t(Unew)) - Uold%*%((1/dold)*t(Uold))
    r.e <- MM%*%t(delta.update) #old MM
    TWQdiffW = -1/(2*tau)*colSums(r.e*(foo1%*%r.e)) # === apply(r.e,2,t(r.e)%*%foo1%*%r.e)
    Detdiff = -1/2*( sum(log(dnew)) -  sum(log(dold)) )
    lf0 <- Detdiff + TWQdiffW
    bmlf<- bm(lf0)
    tac = proc.time()-tic
    list(log.ldiff = lf0,meanlnc = bmlf$est, Unew = Unew,dnew=dnew,tac=tac)    
  }
  
  # Start MCEMG 
  for (i in 1:MCN)
  {
    xbeta <- X%*%beta.update[i,]
    update <- deltaloop(U,d,tau.update[i,],beta.update[i,])
    bmlf <- bm(update$log.ldiff)
    
    tic = proc.time()
    out <- mclapply(c(phiidx-5:1,phiidx+1:5),
                    function(phix) phiupdate(phix,update$taunew,update$deltasample[seq(1,update$mcsize,len=3e3),],U,d),mc.cores=mc.cores)
    tac = proc.time()-tic
    
    loc <- which.max(sapply(out,"[[","meanlnc"))
    if(out[[loc]]$meanlnc > 0){
      phiidx <- c(phiidx-5:1,phiidx+1:5)[loc]
      d <- out[[loc]]$dnew 
      U <- out[[loc]]$Unew  
    } 
   
    mcsize = update$mcsize
    delta.update = update$deltasample
    deltalast = delta.update[mcsize,]
    
    # different adaptation 6/26
    sampvar = var(delta.update)
    proposal = 0.95*2.38^2/q*sampvar + 0.05*diag(0.1^2/q,q)
    proposal = diag(diag(chol(proposal)))
    
    M = update$M
    MM =update$MM
    log.list[[i]] = update$log.l
    accep = rbind(accep,update$accept)
    adaptsize[i]=mcsize
    
    if(i>1) checksize =setsize= mcsize
    beta.update[i+1,]<- update$betanew
    tau.update[i+1,]<- update$taunew
    phi.update[i+1,]<- phigrid[phiidx]
    deltamean.update[i+1,]<- colMeans(delta.update)
    wmean.update[i+1,]<- MM %*% deltamean.update[i+1,]
    if(i == 1) size <- ceil*q # increase mcsize to maximum allowed
    if(track) {
      # quilt.plot(coords,wmean.update[i+1,],main=paste0("RandomEffects Est",i))
      cat("---------------------------------------------- \n")
      cat("Beta Estimates",beta.update[i+1,],"\n")
      cat("Tau Estimates",tau.update[i+1,],"\n")
      cat("phi Estimates",phi.update[i+1,],"\n")
      cat("tau/phi",tau.update[i+1,]/phi.update[i+1,],"\n")
      cat("Acc.rate",range(update$accept),"\n")
      cat("mcsize",mcsize,"\n")
      cat("----END ITERATION---",i,"--------------------- \n")
    }
    
    if(i>=3 & flag){ # save output for stoping iteration
      if(bmlf$est + 2*bmlf$se < epsilon*abs(mean(update$log.l))&
         all(abs(beta.update[i+1,])-abs(beta.update[i,])< epsilon*abs(beta.update[i,]))){
        stoptime <- proc.time() - ptm 
        stopiter <- i
        flag  <- FALSE
        samp <- update$deltasample
        Ustop <- update$U
        dstop <- update$d
        Mstop <- update$M
        MMstop <- update$MM
        if(ifbreak) break
      }
    } 
  }
  
  runtime=proc.time() - ptm
  return(list(log.l=log.list, wmean.update =wmean.update,coords=coords, Z=Z, X=X,covfn = covfn,
              tau.est = tau.update,beta.est = beta.update, phi.est = phi.update, re.samples = delta.update, 
              Ustop=Ustop,dstop=dstop, Mstop=Mstop,MMstop=MMstop,family = family,Ufinal=update$U,dfinal = update$d,
              adaptsize = adaptsize, accep = accep, deltamean = deltamean.update, M = M, MM = MM, q = q,runtime=runtime,
              stoptime = stoptime,stopiter=stopiter,samp = samp,
              betainit = beta.init, tauinit = tau.init, deltainit= delta.init,ceil = ceil, epsilon= epsilon,zalpha=zalpha))
}

#' @export
sparse.sglmmGP.mcem.pred <- function(fit,pred.X,pred.coords,m=1e3,burn=0.5,rdist = rdist){
  require(fields)
  p = nrow(pred.coords)
  if(!is.null(fit$stopiter)){
    stopiter = fit$stopiter+1
    Ustop = fit$Ustop
    dstop = fit$dstop
  }else{
    stopiter = length(fit$adaptsize)
    Ustop = fit$Ufinal
    dstop = fit$dfinal
    fit$samp = fit$re.samples
  }
  
  phi = fit$phi.est[stopiter,]
  betahat = fit$beta.est[stopiter,]
  sigma2 = fit$tau.est[stopiter,]
  covfn = fit$covfn
  
  predidx = seq(nrow(fit$samp)*burn+1,nrow(fit$samp),length=m)
  K = covfn(rdist(fit$coords, pred.coords),phi)
  K11 = covfn(rdist(fit$coords),phi)
  K11inv = solve(K11)
  
  MM =  Ustop%*%diag(sqrt(dstop))
  CondMean = (t(K)%*%K11inv%*%(MM%*%t(fit$samp[predidx,])))
  CondVar = sigma2*(covfn(rdist(pred.coords),phi)-t(K)%*%K11inv%*%K)
  
  # cholCV = chol(CondVar)
  # pred.re = t(matrix(rnorm(m*p), ncol=p)%*%cholCV) + CondMean
  pred.re = t(matrix(rnorm(m*p), ncol=p)%*%chol(diag(diag(CondVar)))) + CondMean
  
  if(fit$family == "poisson"){
    pred.pi =  exp(c(pred.X%*%betahat) + pred.re)
    pred.mean=  exp(c(pred.X%*%betahat) + CondMean)
    pred.Z = apply(pred.pi,2,function(x) rpois(p,lambda = x))
  }else{
    expo =  exp(c(pred.X%*%betahat) + pred.re)
    pred.pi = expo/(1+expo)
    pred.mean = exp(c(pred.X%*%betahat) + CondMean)/(1+exp(c(pred.X%*%betahat) + CondMean))
    pred.Z = apply(pred.pi,2,function(x) rbinom(p,1,x))
  }
  return(list(pred.pi = pred.pi, pred.re =pred.re, pred.Z = pred.Z, pred.X=pred.X, pred.coords = pred.coords,pred.mean=pred.mean))
}

rpK<-function(n,r,K){ # r is demension selected, K,n is loaded from global environment
  omega <- matrix(rnorm(n*r,mean=0,sd=1/sqrt(r)),nrow=n,ncol=r)
  phi <- svd(K%*%omega)$u
  Kphi <-K%*%phi
  K1 <-  t(phi)%*% Kphi
  #   Bt<-chol(K1)
  #   C= K%*%phi%*%solve((Bt))
  svdk1  = svd(K1)
  C= Kphi%*%svdk1$v%*%diag(sqrt(1/svdk1$d))
  foo<-svd(C)
  return(list(u = foo$u,d = foo$d^2))
}

# n = 1000
# coords = matrix(runif(n*2),ncol=2)
# dist = rdist(coords)
# covfn <- covfndef(0.5)
# K = covfn(dist,1.5)
# set.seed(1)
# M0=rpKtest(n,r*2,K)
# set.seed(1)
# M1=rpK(1000,100,K)
# max(abs(M0$d-M1$v))

#' @export
# define covariance function
covfndef <- function(nu){
  # exponential 
  if(nu == 0.5) covfn <- function(dist,phi) { 
    K = exp(-1/phi*dist) 
    return(K)
  }
  # matern 1.5
  if(nu == 1.5) covfn <- function(dist,phi) { 
    K = (1+sqrt(3)/phi*dist)*exp(-sqrt(3)/phi*dist)
    return(K)
  }
  # matern 1.5
  if(nu == 2.5) covfn <- function(dist,phi) { 
    K = (1+sqrt(5)/phi*dist+ 5/(3*phi^2)*dist^2)*exp(-sqrt(5)/phi*dist)
    return(K)
  }
  # square exponential
  if(nu ==10) covfn <- function(dist,phi) { 
    K = exp(-1/(2*phi^2)*dist^2) 
    return(K)
  }
  return(covfn)
}

tr <- function(x) sum(diag(x))

# generate multivariate normal
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

SE<-function(Z,X,offset=NULL,beta,tau,delta.update,M,family,A=NULL){
  # Observed infromation matrix
  # second derivative
  q = ncol(delta.update)
  n = length(Z)
  if(is.null(offset)) offset = rep(1,n)
  QQ<-diag(c(A%*%matrix(1,n,1)))-A # precision matrix as calculated in Hughes and Haran,
  Q <- t(M)%*%QQ%*%M
  p <-ncol(X) # define the number of fixed effects
  I = matrix(0,p+1,p+1)
  xbeta = X%*%beta
  
  if(family=="binomial"){ #add offset
    # define the link function
    link = function(x){
      ExpZ = 1/(1+exp(-x)) # or exp(x)/(1+exp(x))
      VarZ = ExpZ*(1-ExpZ) 
      return(list(ExpZ=ExpZ,VarZ=VarZ))
    }
  }
  
  if(family=="poisson"){
    # define the link function
    link = function(x){
      VarZ = ExpZ = exp(x+log(offset)) # or exp(x)/(1+exp(x))
      return(list(ExpZ=ExpZ,VarZ=VarZ))
    } 
  }
  
  
  r.e <- M%*%t(delta.update)
  foo <- link(c(xbeta) + r.e) # list of two "ExpZ" "VarZ"
  EZ_approx <- rowMeans(foo[[1]]) 
  if(family =="binomial"){VZ_approx<-EZ_approx*(1-EZ_approx)}else{VZ_approx<-EZ_approx}
  
  EZ_conditional <- foo[[1]]
  VZ_conditional <- foo[[2]]
  
  #   foo = apply(delta.update,1, function(delta) link(xbeta + M%*%delta))
  #   foo = do.call("rbind",foo)
  #   EZ_conditional <- do.call("cbind",foo[,1])
  #   EZ_approx = apply(EZ_conditional,1,mean)
  #   VZ_conditional = do.call("cbind",foo[,2])
  #   VZ_approx = apply(VZ_conditional,1,mean)
  
  TWQW_conditional =  apply(delta.update,1,function(x) t(x)%*%Q%*%x)
  TWQW_approx = mean(TWQW_conditional)
  
  I[1:p,1:p] = t(X)%*%(VZ_approx*X)
  I[p+1,p+1] = 1/(2*tau^2/q) 
  
  # first derivative
  S_beta = apply(EZ_conditional,2,function(x) (t(X)%*%(Z-x)))
  S_tau = (q/(2*tau) - 1/2*TWQW_conditional)
  S = rbind(S_beta,S_tau)
  StS = 1/ncol(S)*S%*%t(S)
  
  # following quantity should be around zero.
  # foo should be zero if converged 
  foo= c(t(X)%*%(Z-EZ_approx),(q/(2*tau) - 1/2*TWQW_approx))%*%t(c(t(X)%*%(Z-EZ_approx),(q/(2*tau) - 1/2*TWQW_approx)))
  
  # take inverse of observed infromation matrix
  Obsinfm = solve(I-StS+foo)
  #   Obsinfmfoo = solve(I+cov(t(S)))
  
  # below is the standard error of each point estimates. 
  sqrt(diag(Obsinfm))
  Paramest = c(beta,tau)
  uci = Paramest + 2*sqrt(diag(Obsinfm))
  lci = Paramest - 2*sqrt(diag(Obsinfm))
  return(list(uci = uci,lci=lci,foo=foo,Obsinfm=Obsinfm,I=I,StS=StS))
}

SE_MCEM_GP<-function(Z,X,offset=NULL,beta,tau,phi,delta.update,M,family){
  # Observed infromation matrix
  # second derivative
  q = ncol(delta.update)
  n = length(Z)
  if(is.null(offset)) offset = rep(1,n)
  #   QQ<-diag(c(A%*%matrix(1,n,1)))-A # precision matrix as calculated in Hughes and Haran,
  #   Q <- t(M)%*%QQ%*%M
  p <-ncol(X) # define the number of fixed effects
  I = matrix(0,p+1,p+1)
  xbeta = X%*%beta
  
  if(family=="binomial"){ #add offset
    # define the link function
    link = function(x){
      ExpZ = 1/(1+exp(-x)) # or exp(x)/(1+exp(x))
      VarZ = ExpZ*(1-ExpZ) 
      return(list(ExpZ=ExpZ,VarZ=VarZ))
    }
  }
  
  if(family=="poisson"){
    # define the link function
    link = function(x){
      VarZ = ExpZ = exp(x+log(offset)) # or exp(x)/(1+exp(x))
      return(list(ExpZ=ExpZ,VarZ=VarZ))
    } 
  }
  #
  r.e <- M%*%t(delta.update)
  foo <- link(c(xbeta) + r.e) # list of two "ExpZ" "VarZ"
  EZ_approx <- rowMeans(foo[[1]]) 
  if(family =="binomial"){VZ_approx<-EZ_approx*(1-EZ_approx)}else{VZ_approx<-EZ_approx}
  
  EZ_conditional <- foo[[1]]
  VZ_conditional <- foo[[2]]
  TWQW_conditional =  apply(delta.update,1,function(x) t(x)%*%x)
  TWQW_approx = mean(TWQW_conditional)
  
  # observed information -E(second derive)
  I[1:p,1:p] = t(X)%*%(VZ_approx*X)
  I[p+1,p+1] = -(q/(2*tau^2)-TWQW_approx/tau^3)
  # first derivative
  S_beta = apply(EZ_conditional,2,function(x) (t(X)%*%(Z-x)))
  S_tau = (-q/(2*tau) + TWQW_conditional/(2*tau^2))
  ###
  S = rbind(S_beta,S_tau)
  StS = 1/ncol(S)*S%*%t(S)
  
  # following quantity should be around zero.
  # foo should be zero if converged 
  # YG: why is this not zero
  foo= c(t(X)%*%(Z-EZ_approx),(-q/(2*tau) + TWQW_approx/(2*tau^2)))%*%t(c(t(X)%*%(Z-EZ_approx),(-q/(2*tau) + TWQW_approx/(2*tau^2))))
  
  # take inverse of observed infromation matrix
  Obsinfm = solve(I-StS+foo)
  
  # below is the standard error of each point estimates. 
  sqrt(diag(Obsinfm))
  1/sqrt(diag(I-StS+foo))
  
  Paramest = c(beta,tau)
  uci = Paramest + 2*sqrt(diag(Obsinfm))
  lci = Paramest - 2*sqrt(diag(Obsinfm))
  return(list(uci = uci,lci=lci,foo=foo,Obsinfm=Obsinfm))
}

# compute batchmeans, function from R packages batchmeans
bm <- function (x, size = "sqroot", warn = FALSE){
  n = length(x)
  if (n < 1000) {
    if (warn) 
      warning("too few samples (less than 1,000)")
    if (n < 10) 
      return(NA)
  }
  if (size == "sqroot") {
    b = floor(sqrt(n))
    a = floor(n/b)
  } else if (size == "cuberoot") {
    b = floor(n^(1/3))
    a = floor(n/b)
  } else {
    if (!is.numeric(size) || size <= 1 || size == Inf) 
      stop("'size' must be a finite numeric quantity larger than 1.")
    b = floor(size)
    a = floor(n/b)
  }
  y = sapply(1:a, function(k) return(mean(x[((k - 1) * b + 1):(k * b)])))
  mu.hat = mean(y)
  var.hat = b * sum((y - mu.hat)^2)/(a - 1)
  se = sqrt(var.hat/n)
  (list(est = mu.hat, se = se, var.hat = var.hat))
}

bmmat <- function(x) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("'x' must be a matrix or data frame.")
  num = ncol(x)
  bmvals = matrix(NA, num, 2)
  colnames(bmvals) = c("est", "se")
  rownames(bmvals) = colnames(x)
  bmres = apply(x, 2, bm)
  for (i in 1:num) bmvals[i, ] = c(bmres[[i]]$est, bmres[[i]]$se)
  bmvals
}

AdjustInf <- function(MM,beta,re.samples,X){
  AP <- solve(t(X)%*%(X))%*%t(X)
  w = MM%*%t(re.samples)
  Ajust <- AP%*%w
  Abeta <- beta - apply(Ajust,1,post.mode)
  Abetadistri <- beta - Ajust
  return(list(Abeta = Abeta,Abetadistri=Abetadistri, Ajust = rowMeans(Ajust)))
}

logit.inv <-function(x) log(x) - log(1-x)

post.mode <- function(s){
  d <- density(s)
  d$x[which.max(d$y)]
}
#' @export
quants <- function(x){
  CI = quantile(x, prob = c(0.025,0.975))
  return(c(mean = mean(x), CI[1],CI[2]))
}
