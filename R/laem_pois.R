#' Fit LAEM to projection-based SGLMM
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
sparse.sglmmGP.laem <- function(Z,X,family="poisson",offset = NULL,q,MCN,size,nu = nu,covfn=covfn,rdist=rdist,
                                mc.cores=1,coords = coords,mul=2,init=NULL,ceil = ceil, epsilon=1e-3,
                                zalpha = zalpha, ifbreak = ifbreak,RSR=RSR,RP=RP,tune=0.1,track=T,
                                tol = 1e-4, order= 1,optiter = 10){
  library(fields)
  ptm = proc.time()
  n   = length(Z)
  p     <- ncol(X) # define the number of fixed effects
  PPERP <- diag(n)- X%*%chol2inv(chol(crossprod(X,X)))%*%t(X)
  dist  <- rdist(coords,coords)

  flag = T
  stopiter=Ustop =dstop = deltastop = stoptime= NULL

  # set up initials for MCEMG --------------
  rhomax <- Matern.cor.to.range(1.5*max(dist), nu=nu, cor.target=.05)
  rhomin <- Matern.cor.to.range(0.1*max(dist), nu=nu, cor.target=.05)
  phigrid = seq(rhomin*sqrt(2*nu),rhomax*sqrt(2*nu),by=0.01)
  cat("min range", round(min(phigrid),3), "; max range",round(max(phigrid),3),"\n")
  # phigrid = seq(0.1,1.5,by=0.01)
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
    K = rpK(n,2*q,K)
    U = K$u[,1:q]
    d = K$d[1:q]
    # call c function for computing eigen
    # K = .Call("rp",phi.init,0,coords,as.integer(n) ,as.integer(rk),nu,as.integer(core))
    # K.rp = list(d = K[[1]],u = K[[2]][,1:q])
    # d <- (K.rp$d[1:q])^2
    # U1<- U <- K.rp$u[,1:q]
    # U <- .Call("mmC",PPERP,U,as.integer(n),as.integer(q),as.integer(core))
  }
  MM <- U%*%diag(sqrt(d))

  g.lm <-glm(Z~X-1,family=family) # YG: add weights for binomial
  if(is.null(init$beta.init)){
    beta.init<- g.lm$coefficients
    w.init   <-g.lm$residuals # residuals taken as random effects
    tau.init = var(w.init)
    deltalast<-delta.init<-rep(0,q) # t(MM)%*%w.init
    bin1     <- likfit(list(coords=coords,data=w.init), ini = c(tau.init,0.2), fix.nugget = T, fix.kappa = T,kappa=0.5, lik.method = "ML", cov.model="matern")
    phi.init = bin1$phi
    tau.init = bin1$sigmasq
    phiidx = which.min(abs(phigrid-phi.init));phi.init = phigrid[phiidx]
  }else{
    beta.init<- init$beta.init
    deltalast<-delta.init<-rep(0,q)
    w.init   <-rep(0,n)
    tau.init = init$tau.init
    phi.init = init$phi.init
    phiidx = which.min(abs(phigrid-phi.init));phi.init = phigrid[phiidx]
  }

  # matrices to store some output from the algorithm -------
  tau.update <-matrix(NA,nrow=MCN+1,ncol=1)
  phi.update <-matrix(NA,nrow=MCN+1,ncol=1)
  beta.update<-matrix(NA,nrow=MCN+1,ncol=p)
  deltamean.update<-matrix(NA,nrow=MCN+1,ncol=q) # matrix to store estimates of random effect
  wmean.update    <-matrix(NA,nrow=MCN+1,ncol=n)
  wmean.update    <-matrix(NA,nrow=MCN+1,ncol=n)
  wmean.update    <-matrix(NA,nrow=MCN+1,ncol=n)
  U.update   <- array(NA, dim = c(n,q,MCN+1))
  d.update   <- matrix(NA,nrow=q,ncol=MCN+1)

  log.lf<-c()
  log.lf[1] <- NA
  tau.update[1,]<-tau.init
  phi.update[1,]<-phi.init
  beta.update[1,]<-beta.init
  deltamean.update[1,]<-0
  wmean.update[1,]<-0

  # 1/17 added Gaussian approximation
  LAfn <-function(M,tau,xbeta,q){ # M is the projection matrix, either PperpUD or UD
    v=rep(0,q)
    i = 0
    diff = 10
    while(diff>tol & i < optiter){
      # poisson different than binomial
      d = D = c(exp(xbeta+M%*%v))
      DM = D*M
      Q = crossprod(M,DM)  # t(M)%*%D%*%M + Sigmainv
      diag(Q) = diag(Q) + rep(1/tau,q)
      b = crossprod(M, Z - d + DM%*%v)
      vnew = solve(Q,b) # this will converge
      diff = crossprod(v-vnew)
      v = vnew
      i = i+1
    }
    return(list(v=vnew,M=M,Q=Q,b=b))
  }

  # 1/23 added Gaussian approximation for w.
  # LAfnw <-function(tau,N=1,xbeta,n,U,dvec){ # M is the projection matrix, either PperpUD or UD
  #   v = rep(0,n)
  #   i = 0
  #   cat("here",n,"\n")
  #   diff = 10
  #   while(diff>tol){
  #     D = c(N*exp(xbeta+v)/(1+exp(xbeta+v))^2)
  #     d = N*exp(xbeta+v)/(1+exp(xbeta+v))
  #     Q = U%*%diag(1/dvec)%*%t(U)/tau # t(H)%*%D%*%H + Sigmainv,for this H is Identity.
  #     diag(Q) = diag(Q) + D
  #     b = Z - d + D*v
  #     vnew = solve(Q,b) # this will converge
  #     diff = crossprod(v-vnew)
  #     v = vnew
  #     i = i+1
  #   }
  #   cat("Took itertaion",i,"to converage","\n")
  #   return(list(v=vnew,Q=Q,b=b))
  # }
  # LAapprox <- LAfn(U[,1:30]%*%diag(sqrt(d[1:30])),tau.update[i],N=1,xbeta,30) #return v=vnew,M=M,Q=Q,diff=diff, b=b
  # v1 = U[,1:30]%*%diag(sqrt(d[1:30]))%*%LAapprox$v
  # LAapproxw <- LAfnw(U[,1:30]%*%diag(sqrt(d[1:30])),tau.update[i],N=1,xbeta,n,U,d) #return v=vnew,M=M,Q=Q,diff=diff, b=b
  # plot(v1,LAapproxw$v)

  # gaussian approx to conditional mean
  LACondM <- function(xbeta,M,delta,Q){
    m <- rep(NA,n)
    expo <- exp(xbeta+ M%*%delta)
    # poisson different than binomial
    D = expo
    # deriv2 = t(M)%*%diag(D)%*%M
    for(i in 1:n){ # rewrite the for loop later
      # poisson different than binomial
      m[i]<-expo[i]+ 1/2*D[i]*(M[i,]%*%solve(Q,M[i,]))
    }
    return(m)
  }

  # gaussian approx to conditional var
  # poisson different than binomial
  # LACondV is the same as LACondM
  # LACondV <- function(xbeta,M,delta,Q){
  #   v <- rep(NA,n)
  #   expo <- exp(xbeta+ M%*%delta)
  #   D = expo*(expo^2-4*expo+1)/(expo+1)^4
  #   # deriv2 = t(M)%*%diag(D)%*%M
  #   for(i in 1:n){ # rewrite the for loop later
  #     v[i]<-expo[i]/(expo[1]+1)^2 + 1/2*D[i]*(M[i,]%*%solve(Q,M[i,]))
  #   }
  #   return(v)
  # }

  phiupdate0 <- function(phix,xbeta,tau,delta,Q,MMold,q){ # xbeta MM tau delta.update can pulled outside from here
    if(phix<=0|phix>=length(phigrid)) return(list(log.ldiff=NA))
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

    MMnew = Unew%*%diag(sqrt(dnew))
    signdiag = sign(diag(t(MMnew)%*%MMold)) # find the sign change, MM is loaded from global environment
    signdiag = as.logical(1-signdiag)
    MMnew[,signdiag] = -MMnew[,signdiag]
    #change sign of Unew then below
    M  = PPERP%*%MMnew
    ##################################################
    # compute the difference in likelihood function due to phi
    # -1/(2sigma^2) w'(UDinvU_new - UDinvU_old)W
    ##################################################
    y = xbeta+M%*%delta
    prob=exp(y)/(1+exp(y))
    lf = t(Z)%*%y - sum(log(1+exp(y))) - 1/(2*tau)*t(delta)%*%delta - q/2*log(tau) +
      1/2*t(prob*(1-prob))%*%diag(M%*%Q%*%t(M)) - 1/2*sum(diag(Q))*tau
    return(list(lf=lf,MM=MMnew))
  }

  phiupdate <- function(phix,tau,delta,deltaQ,Uold,dold){ # xbeta MM tau delta.update can pulled outside from here
    tic = proc.time()
    if(phix<=0|phix>=length(phigrid)) return(list(log.ldiff=NA))
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
    ##################################################
    # compute the difference in likelihood function due to phi
    # -1/(2sigma^2) w'(UDinvU_new - UDinvU_old)W
    ##################################################
    foo1 = Unew%*%((1/dnew)*t(Unew)) - Uold%*%((1/dold)*t(Uold))
    r.e <- MM%*%delta #old MM
    TWQdiffW = -1/(2*tau)*sum(diag((solve(deltaQ) + delta%*%t(delta))%*%(t(MM)%*%foo1%*%MM)))
    Detdiff = -1/2*(sum(log(dnew)) -  sum(log(dold)))
    lf0 <- Detdiff + TWQdiffW
    tac = proc.time()-tic
    list(log.ldiff = lf0,Unew = Unew,dnew=dnew,tac=tac)
  }

  tau = tau.update[1]


  # Start MCEMG
  for (i in 1:MCN)
  {
    if(RSR) M <- PPERP%*%MM else M <- MM
    xbeta <- X%*%beta.update[i,]
    # 1/17 added Gaussian approximation
    LAapprox <- LAfn(M,tau.update[i],xbeta,q) #return v=vnew,M=M,Q=Q,diff=diff, b=b
    delta <- LAapprox$v
    deltaQ <- LAapprox$Q
    # poisson different than binomial
    VZ_approx <- EZ_approx <- LACondM(xbeta,M,delta,deltaQ)
    TWQW_approx <- sum(diag(solve(deltaQ) + delta%*%t(delta))) # this is E delta'delta

    # update beta and tau
    betanew <- beta.update[i,] + solve(t(X)%*%(VZ_approx*X),t(X)%*%(Z-EZ_approx))

    if(order == 1 |i ==1 ){
      # first order update
      taunew <- 1/q * TWQW_approx #***double check this***
    }else{# # second order update
      taunew <- (tau -(q/(2*tau^2)-TWQW_approx/tau^3)^(-1)*
                   (-q/(2*tau) + TWQW_approx/(2*tau^2))
      )

      # if(taunew<=0) cat("Update Tau",taunew,"\n")
      KK  = 1
      while(taunew <=0){
        taunew <- (tau - 1/(KK*2)*(q/(2*tau^2)-TWQW_approx/tau^3)^(-1)*
                     (-q/(2*tau) + TWQW_approx/(2*tau^2)));KK = KK+1
      }
      if(KK>1){
        beta.update[i+1,]<- beta.update[i,] + 1/((KK-1)*2)*solve(t(X)%*%(VZ_approx*X),
                                                                 t(X)%*%(Z-EZ_approx))
      }
    }

    # # # update phi YG: continue here
    out <- mclapply(c(phiidx-5:1,phiidx+1:5),
                    function(phix) phiupdate(phix,tau.update[i],delta,deltaQ,U,d),mc.cores=mc.cores)
    # phiupdate <- function(phix,tau,delta,deltaQ,Uold,dold)
    # loc <- which.max(sapply(out,"[[","lf"))
    # if(out[[loc]]$lf - lf0$lf > 0){
    loc <- which.max(sapply(out,"[[","log.ldiff"))
    if(out[[loc]]$log.ldiff > 0){
      phiidx <- c(phiidx-5:1,phiidx+1:5)[loc]
      # MM <- out[[loc]]$MM
      U <- out[[loc]]$Unew
      d <- out[[loc]]$dnew
      MM   <- U%*%diag(sqrt(d))
    }

    beta.update[i+1,]<- betanew
    tau<-tau.update[i+1,] <- taunew
    phi.update[i+1,] <- phigrid[phiidx]
    deltamean.update[i+1,]<- delta
    wmean.update[i+1,]<- MM %*% delta
    U.update[,,i+1]   <- U
    d.update[,i+1]   <- d
    # poisson different than binomial
    xbeta <- X%*%beta.update[i,]
    log.lf[i+1] <- Z%*%(xbeta+M %*% delta)-sum(exp(xbeta+M%*%delta)) -
      q/2*log(tau.update[i,]) - 1/(2*tau.update[i,])*crossprod(delta,delta) -
      1/2*sum(diag(solve(deltaQ,t(M)%*%(c(exp(xbeta+M %*% delta))*M) + diag(1/tau.update[i,],q))))

    if(i>=3 & flag){ # save output for stoping iteration
      if((log.lf[i+1]-log.lf[i] < epsilon*abs(log.lf[i]))&
         all(abs(abs(beta.update[i+1,])-abs(beta.update[i,]))< epsilon*abs(beta.update[i,]))){
        stoptime <- proc.time() - ptm
        stopiter <- i+1
        flag <- FALSE
        Ustop <- U
        dstop <- d
        deltastop <- delta
        if(ifbreak) break
      }
    }

    if(track){
      cat("---------------------------------------------- \n")
      cat("Beta Estimates",beta.update[i+1,],"\n")
      cat("Tau Estimates",tau.update[i+1,],"\n")
      cat("phi Estimates",phi.update[i+1,],"\n")
      cat("tau/phi",tau.update[i+1,]/phi.update[i+1,],"\n")
      cat("----END ITERATION---",i,"--------------------- \n")
    }

  }

  # I = matrix(0,p+1,p+1) # E(-d^2 Lc | Z) = -Q''
  # # complete information -E(second derive)
  # I[1:p,1:p] = t(X)%*%(VZ_approx*X)
  # I[p+1,p+1] = -(q/(2*tau^2)-TWQW_approx/tau^3)
  # # missing information delta|y = E(SS'|Z)
  # txz = t(X)%*%Z
  # y = X%*%betanew + M%*%delta
  # prob=exp(y)/(1+exp(y))
  # S_beta = txz%*%t(txz)-txz%*%t(EZ_approx)%*%X -
  #   t(txz%*%t(EZ_approx)%*%X) + t(X)%*%(prob%*%t(prob))%*%X # the last term is a plug in estimate
  # # S_tau = (-q/(2*tau) + TWQW_conditional/(2*tau^2))
  # ###
  # # S = rbind(S_beta,S_tau)
  # # StS = 1/ncol(S)*S%*%t(S)
  #
  # # following quantity should be around zero.
  # # foo should be zero if converged
  # # YG: why is this not zero
  # # foo= c(t(X)%*%(Z-EZ_approx),(-q/(2*tau) + TWQW_approx/(2*tau^2)))%*%t(c(t(X)%*%(Z-EZ_approx),(-q/(2*tau) + TWQW_approx/(2*tau^2))))
  #
  # # take inverse of observed infromation matrix
  # Obsinfm = solve(I[1:p,1:p]-S_beta)
  #
  # # below is the standard error of each point estimates.
  # SE = sqrt(diag(Obsinfm))

  runtime=proc.time() - ptm
  return(list(time=runtime,log.l=log.lf,family=family, M = M, MM = MM, U = U,d = d, delta=delta,deltaQ = deltaQ,
              wmean.update =wmean.update, coords=coords, Z=Z, X=X,covfn = covfn, stopiter= stopiter, Ustop = Ustop,
              dstop = dstop,stoptime = stoptime,deltastop=deltastop,
              tau.est = tau.update,beta.est = beta.update, phi.est = phi.update,delta.update = deltamean.update,
              U.update = U.update,d.update=d.update))
}

#' @export
sparse.sglmmGP.laem.pred <- function(fit,pred.X,pred.coords,m,rdist=rdist){
  require(fields)
  p = nrow(pred.coords)
  if(!is.null(fit$stopiter)){
    stopiter = fit$stopiter
    Ustop = fit$U.update[,,stopiter]
    dstop = fit$d.update[,stopiter]
    delta = fit$delta.update[stopiter,]
  }else{
    stopiter = length(fit$phi.est)
    Ustop = fit$U
    dstop = fit$d
    delta = fit$delta
  }

  phi = fit$phi.est[stopiter,]
  betahat = fit$beta.est[stopiter,]
  sigma2 = fit$tau.est[stopiter,]
  covfn = fit$covfn

  K = covfn(rdist(fit$coords, pred.coords),phi)
  K22 = covfn(rdist(pred.coords),phi)

  # lowrank prediction MS
  CondMean = t(K)%*%(Ustop%*%diag(sqrt(1/dstop))%*%delta)
  CondVar = sigma2*(K22-t(K)%*%Ustop%*%diag(1/dstop)%*%t(Ustop)%*%K)

  # fullrank prediction
  # K11 = covfn(rdist(fit$coords),phi)
  # K11inv = solve(K11)
  # MM =  Ustop%*%diag(sqrt(dstop))
  # CondMean = (t(K)%*%K11inv%*%(MM%*%delta))
  # CondVar = sigma2*(K22-t(K)%*%K11inv%*%K)

  cholCV = chol(CondVar)
  # simulate mulrivar
  pred.re = t(matrix(rnorm(m*p), ncol=p)%*%cholCV) + matrix(CondMean, ncol = m,nrow =p)
  # simulate individually
  # pred.re = t(matrix(rnorm(m*p), ncol=p)%*%chol(diag(diag(CondVar)))) + matrix(CondMean, ncol = m,nrow =p)

  if(fit$family == "poisson"){
    pred.pi =  exp(c(pred.X%*%betahat) + pred.re)
    pred.mean =  exp(c(pred.X%*%betahat) + CondMean)
    pred.Z = apply(pred.pi,2,function(x) rpois(p,lambda = x))
  }else{
    expo =  exp(c(pred.X%*%betahat) + pred.re)
    pred.pi = expo/(1+expo)
    pred.mean =   exp(c(pred.X%*%betahat) + CondMean)/(1+ exp(c(pred.X%*%betahat) + CondMean))
    pred.Z = apply(pred.pi,2,function(x) rbinom(p,1,x))
  }
  return(list(pred.pi = pred.pi, pred.re =pred.re, pred.Z = pred.Z, pred.X=pred.X, pred.coords = pred.coords,pred.mean=pred.mean))
}


SE_LAEM_GP<-function(Z,X,offset=NULL,beta,tau,phi,delta.update,M,family){
  # Observed infromation matrix
  # second derivative
  q = length(delta.update)
  n = length(Z)
  if(is.null(offset)) offset = rep(1,n)
  p <- ncol(X) # define the number of fixed effects
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
  r.e <- M%*%delta.update
  foo <- link(c(xbeta) + r.e) # list of two "ExpZ" "VarZ"
  EZ_approx <- c(foo[[1]])
  if(family =="binomial"){VZ_approx<-EZ_approx*(1-EZ_approx)}else{VZ_approx<-EZ_approx}

  EZ_conditional <- foo[[1]]
  VZ_conditional <- foo[[2]]
  TWQW_conditional =  t(delta.update)%*%delta.update
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
