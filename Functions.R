############################################################################
#####    codes for survival Analysis via spatial- and temporal-smoothing
##############################################################################

estimate.beta <- function(dat,  p, q, beta0, lambda1, lambda2, weight, maxiter, epsilon){
  betap <- array(0,dim=c(p,q,m))
  betaold <- beta0
  betanew <- betaold
  diff <- 1
  iter <- 0
  while (diff > epsilon){
    log_like<-rep(NA,p*q)
    iter <- iter + 1
    print(paste("iter =" , iter))
    likelihood <- array(NA,dim=c(p*q))
    for (i in 1:p){
      for (j in 1:q){
        datij <- dat[[i]][[j]]
        yij <- as.numeric(datij$y)
        deltaij <- as.numeric(datij$failed)
        Xij <- as.matrix(subset(datij, select = -c(y, failed)))
        
        ord <- order(yij)
        revord <- order(ord)
        ordyij <- yij[ord] 
        ordXij <- Xij[ord,]
        orddeltaij <- deltaij[ord]
      
        NResult<- NRupdate(ordyij, ordXij, orddeltaij, betaold, weight,
                           lambda1, lambda2, i, j, alpha = 0.1)
        betaijord <- NResult[[1]]
        log_like[which(is.na(log_like))[1]] <-NResult[[2]]
        betanew[i,j,] = betaijord
      }
    }
    
    
    for (a in 1:m) {
      betap[,,a]<-sum((t(diff(t(betanew[,,a]))))^2)
    }
    p1<-sum(betap[1,1,])/2*lambda1
    
    p2 <- array(NA,dim=c(p*(p-1)/2))
    for (c in 1:(p-1)) {
      for (d in (c+1):p) {
        p2[which(is.na(p2))[1]]<-sum(weight[c,d]*(betanew[c,,]-betanew[d,,])^2)
      }
    }
    
    p2<-lambda2*sum(p2)/2    
    print(paste(" objfunc=" , (sum(log_like,na.rm = TRUE)-p1-p2)))
    betaoldv <- apply(betaold, 3, as.vector)
    betanewv <- apply(betanew, 3, as.vector)
    diff <- norm(betanewv- betaoldv, "F")
    betaold <- betanew
    if (iter > maxiter){
      print("The maximum iterations is reached")
      break
    }
  }
  return(betanew)
}


cox.logLik <- function(time, status, X, beta){
  if (any(diff(time)<0)){
    cat("the times must in ascending")
    stop
  }
  D.set <- unique(time)
  n <- length(time)
  m <- ncol(X)
  loglik <- 0
  grad <- rep(0,m)
  hessian <- matrix(0,m,m)
  eta <- X%*%beta
  w <- as.vector(exp(eta))
  W <- rev(cumsum(rev(w)))
  loglik <- sum(status*(eta -log (W)))
  Wx <- apply(w*X, 2, function(x){rev(cumsum(rev(x)))})
  grad <- colSums(status*(X - Wx/W))
  ##hessian
  Wx2.vec <- w*t(apply(X, 1, function(x){crossprod(t(x))}))
  Wx2 <- apply(Wx2.vec, 2, function(x){rev(cumsum(rev(x)))})
  W2x <- t(apply(Wx/W, 1,function(x){crossprod(t(x))}))
  hessian.vec <- colSums (-status*(Wx2/W- W2x))
  hessian <- matrix(hessian.vec, m, m)
  res <- list(logLik= as.numeric(loglik), gradient = as.vector(grad), hessian = hessian)
  return(res)
}



NRupdate <- function(yij, Xij, deltaij, betaold, weight, lambda1, lambda2, i, j, alpha){
  n <- length(yij)
  eta <- Xij %*% betaold[i,j,]
  #cox.log <- coxpen::cox_log(y  = yij, delta  = deltaij, X  = Xij, beta  = betaold[i,j,]) 
  #if you run the previous line, you must pre-install package "coxpen" that was uploaded and named "coxpen_1.0.tar.gz". This package works the same as the function "cox.logLik", but it can speed up the computation.
  
  cox.log <- cox.logLik(time = yij, status = deltaij, X = Xij, beta = betaold[i,j,])
  likelihood<- cox.log$logLik
  mu <- cox.log$gradient
  H <- cox.log$hessian
  diagsH <- diag(x = diag(H),m,m)
  if(any(diag(H) ==0) | any(is.na(H))){
    up = diag(1, m, m)
    mu = rep(0,m)
    diagsH  = diag(1,m,m)
  }
  if (j ==1){
    up2 <- lambda1 + lambda2*(sum(weight[i,]) - weight[i,i])
    down2 <- lambda1*betaold[i,j+1,]+ rowSums(sapply((1:p)[-i], function(a){weight[a,i]*betaold[a,j,]}))
  }else if (j == q){
    up2 <- lambda1 + lambda2*(sum(weight[i,]) - weight[i,i])
    down2 <- lambda1*betaold[i,j-1,]+ rowSums(sapply((1:p)[-i], function(a){weight[a,i]*betaold[a,j,]}))
  }else{
    up2 <- 2*lambda1 + lambda2*(sum(weight[i,]) - weight[i,i])
    down2 <- lambda1*(betaold[i,j-1,] +  betaold[i, j+1,])+ rowSums(sapply((1:p)[-i], function(a){weight[a,i]*betaold[a,j,]}))
  }
  
  up <- diagsH - diag(x = up2 ,m,m)
  down <- as.numeric(diagsH %*% betaold[i,j,]) - mu - down2
  new <- solve(up) %*% matrix(down, nrow = m)
  betaijnew <- betaold[i,j,] + alpha*(new  -  betaold[i,j,])
  return(NRupdate<-list(betaijnew,likelihood))
}


############################################################################
#####    codes for the separate method
############################################################################

estimate.beta.cox <- function(dat,  p, q, m){
  
  beta <- array(NA, dim=c(p,q,m))
  for (i in 1:p){
    for (j in 1:q){
      datij <- dat[[i]][[j]]
      yij <- as.numeric(datij$y)
      deltaij <- as.numeric(datij$failed)
      Xij <- as.matrix(subset(datij, select = -c(y, failed)))
      ord <- order(yij)
      revord <- order(ord)
      ordyij <- yij[ord] 
      ordXij <- Xij[ord,]
      orddeltaij <- deltaij[ord]
      
      fit = survival::coxph(Surv(ordyij, orddeltaij) ~ ordXij)
      beta[i,j, ] = fit$coefficients
      
    }
  }
  for (h in 1:m) {
    for (i in 1:p) {
      beta[i,which(is.na(beta[i,,h])),h]=0
    }
    
  }
  return(beta)
}


##############################################################################################################################
#####    codes for data generating
##############################################################################################################################


###################Case 1
gdat1=function(n,m,beta_true){
  #1. The generation of covariates
  
  Z <- rmvnorm(n, mean = rep(0,m), sigma = diag(x= 1, m,m))
  
  
  
  #2. The generation of failure time
  u=runif(n)
  v=(-3*log(1-u)/exp(Z%*%beta_true))-0.5^3
  TT=rep(NA,n)
  for(i in 1:n)
  {
    if(v[i]>0)
      TT[i]=v[i]^(1/3)+0.5
    else
      TT[i]=-(-v[i])^(1/3)+0.5
  }
  
  
  
  #3. The generation of censoring time
  C=rexp(n,0.123)
  #C=pmin(C,ptau)
  #4. The generation of observed time and delta.
  Y=pmin(TT,C)
  Delta=(TT<=C)
  
  ######censor rate###
  cr=1-sum(Delta)/n
  #  q=quantile(TT)
  X=as.data.frame(Z)
  part1<-data.frame(y=Y,X=Z,failed=Delta)
  dat=list(part1,cr)
  return(dat)
}

###################Case 2
gdat2=function(n,m,beta_true)
{
  #1. The generation of covariates
  
  Z1 <- rmvnorm(n, mean = rep(0,m), sigma = diag(x= 1, m,m))
  Z2 <- matrix(NA,n,m)
  for (v in 1:m) {
    u=runif(1,0,1)
    Z2[,v]=rbinom(n,1,u)
  }
  
  Z<- cbind(Z1[,c(1:5)],Z2[,c(6:10)])
  
  #2. The generation of failure time
  u=runif(n)
  v=(-3*log(1-u)/exp(Z%*%beta_true))-0.5^3
  TT=rep(NA,n)
  for(i in 1:n)
  {
    if(v[i]>0)
      TT[i]=v[i]^(1/3)+0.5
    else
      TT[i]=-(-v[i])^(1/3)+0.5
  }
  
  
  
  #3. The generation of censoring time
  C=rexp(n,0.12)
  #C=pmin(C,ptau)
  #4. The generation of observed time and delta.
  Y=pmin(TT,C)
  Delta=(TT<=C)
  
  ######censor rate###
  cr=1-sum(Delta)/n
  #  q=quantile(TT)
  X=as.data.frame(Z)
  part1<-data.frame(y=Y,X=Z,failed=Delta)
  dat=list(part1,cr)
  return(dat)
}

###################Case 3######################
gdat3=function(n,m,beta_true)
{
  #1. The generation of covariates
  
  Z <- matrix(NA,n,m)
  for (v in 1:m) {
    u=runif(1,0,1)
    Z[,v]=rbinom(n,1,u)
  }
  
  
  #2. The generation of failure time
  u=runif(n)
  v=(-3*log(1-u)/exp(Z%*%beta_true))-0.5^3
  TT=rep(NA,n)
  for(i in 1:n)
  {
    if(v[i]>0)
      TT[i]=v[i]^(1/3)+0.5
    else
      TT[i]=-(-v[i])^(1/3)+0.5
  }
  
  
  
  #3. The generation of censoring time
  C=rexp(n,0.123)
  #C=pmin(C,ptau)
  #4. The generation of observed time and delta.
  Y=pmin(TT,C)
  Delta=(TT<=C)
  
  ######censor rate###
  cr=1-sum(Delta)/n
  #  q=quantile(TT)
  X=as.data.frame(Z)
  part1<-data.frame(y=Y,X=Z,failed=Delta)
  dat=list(part1,cr)
  return(dat)
}









