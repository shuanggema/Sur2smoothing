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


##############################################################################################################################
#####    codes for giving count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##############################################################################################################################


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


##############################################################################################################################
#####    codes for arrangement of real data
##############################################################################################################################

arrange_brain <- function(p, q, data){
  realdat2 <- list()
  realdat1<- list()
  for (i in 1:p) {
    for (j in 1:q) {
      a=brain[[i]][[j]]
      colnames(a)[18] <- "y"
      error<-runif(nrow(a),0.0001,0.001)
      a[,18] <- a[,18]+error
      colnames(a)[19] <- "failed"
      realdat1[[j]] <- a[,-c(1,9)]
    }
    realdat2[[i]] <- realdat1
  }
  return(realdat2)
}


##############################################################################################################################
#####    codes for generating weights in real data
##############################################################################################################################

sites <- function(longitude,latitude){
  locations<-data.frame(-longitude,latitude)
  locations<-as.matrix(locations)
  wei = as.matrix(sqrt(2)-dist(locations)/max(dist(locations)))
  diag(wei) <- 1
  return(wei)
}

##############################################################################################################################
#####    codes for bootstrap
##############################################################################################################################


bootstrap <- function(data,times){
  total=data
  totalsize<-nrow(total)
  samplesize<-totalsize-1000
  
  bootresult1<-list()
  cox1<-list()
  ######bootstrap
  for (n in 1:times) {
    ##########sample
    flag<-sample(1:totalsize,samplesize,replace = FALSE)
    samp<-total[flag,]
    
    i=1
    dat1<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=2
    dat2<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=3
    dat3<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=4
    dat4<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=5
    dat5<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=6
    dat6<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=7
    dat7<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=8
    dat8<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=9
    dat9<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
               samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
               samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=10
    dat10<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=11
    dat11<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=12
    dat12<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=13
    dat13<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=14
    dat14<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=15
    dat15<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    i=16
    dat16<-list(samp[samp$location==i & samp$time==1,],samp[samp$location==i & samp$time==2,],samp[samp$location==i & samp$time==3,],
                samp[samp$location==i & samp$time==4,],samp[samp$location==i & samp$time==5,],samp[samp$location==i & samp$time==6,],
                samp[samp$location==i & samp$time==7,],samp[samp$location==i & samp$time==8,],samp[samp$location==i & samp$time==9,])
    
    
    samp<-list(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11,dat12,dat13,dat14,dat15,dat16)
    
    
    
    ###########analysis
    m <- 15
    p <- 16
    q <- 9
    maxiter <- 1500
    epsilon <- 0.001
    
    realdat <- list()
    realdat1<- list()
    for (i in 1:p) {
      for (j in 1:q) {
        a=samp[[i]][[j]]
        colnames(a)[18] <- "y"
        error<-runif(nrow(a),0.0001,0.001)
        a[,18] <- a[,18]+error
        colnames(a)[19] <- "failed"
        realdat1[[j]] <- a[,-c(1,9)]
      }
      realdat[[i]] <- realdat1
    }
    
    beta_cox <- estimate.beta.cox(dat = realdat,  p = p, q = q, m = m)
    cox1[[n]]=beta_cox
    
    for (h in 1:m) {
      for (i in 1:p) {
        beta_cox[i,which(is.na(beta_cox[i,,h])),h]=0
      }
      
    }
    
    beta_est <- estimate.beta(dat = realdat, p = p, q =q,  beta0 = beta_cox, lambda1 = 40, 
                              lambda2 = 1.2, weight = weight, maxiter = maxiter, epsilon = epsilon)
    
    bootresult1[[n]]=beta_est
  }
  return(list(cox1,bootresult1))
}

























