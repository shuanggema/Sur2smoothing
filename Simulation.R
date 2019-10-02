#####################################################################################################
#####
#####    codes for the simulation of survival Analysis via spatial- and temporal-smoothing
##### 
#####################################################################################################

source(file = "Functions.R")
library(mvtnorm)
library(coxed)
library(survival)
library(foreach)
library(doParallel)

########################################set parameters###############################################
m = 10                #number of covariates
p = 9                 #number of locations
q = 10                #number of time intervals (or points)
n = 400               #sample size,sample size=200, 400, 1000, and 2000
maxiter = 200         #the maximum number of iterations
epsilon = 0.001       #set the difference between two consecutive estimates
nsim = 200            #200 replicated
k<-5                  #set k-fold cross validation

########################################generate weights corresponding to locations##################
 
points <- matrix(runif(p*2, 0,1), p ,2)
weight <- as.matrix(sqrt(2)-dist(points))
diag(weight) <- 1

########################################beta generating###############################################

##############Scenario 1######################
beta <- array(NA, dim = c(p,q,m))
a <- runif(m, 3, 4)

for (cov in seq(m)) {
  beta[,1,cov] <- (sqrt((sqrt(2))^2-(sqrt(points[,1]^2+points[,2]^2))^2))/a[cov]+0.2
  for (i in 1:p) {
    beta[i,1:q,cov] <- 0.2*sin(seq(0,2*pi, length.out = q))+ beta[i,1,cov]
  }
}

##############Scenario 2######################
beta <- array(NA, dim = c(p,q,m))
a = runif(m, 3, 4)

for (cov in seq(m)) {
  beta[,1,cov] = (sqrt((sqrt(2))^2-(sqrt(points[,1]^2+points[,2]^2))^2))/a[cov]
  for (i in 1:p) {
    beta[i,1:q,cov] = beta[i,1,cov]+log10(seq(1,q, length.out = q))*0.5
  }
}


##############Scenario 3######################
beta <- array(NA, dim = c(p,q,m))
a = runif(m, 3, 4)

for (cov in seq(m)) {
  beta[,1,cov] = (sqrt((sqrt(2))^2-(sqrt(points[,1]^2+points[,2]^2))^2))/a[cov]
  for (i in 1:p) {
    beta[i,2:q,cov] <- 0.015*(seq(2,q, length.out = q-1)-5)^2+beta[i,1,cov]
    beta[i,1,cov]=beta[i,1,cov]+0.015*(1-5)^2
  }
}


##############Scenario 4######################
beta <- array(NA, dim = c(p,q,m))
a = seq(1.3,4, length.out = m)
#a = c(0.2661695, 0.2159384, 0.3612777)

for (cov in seq(m)) {
  beta[,1,cov] = (sqrt((sqrt(2)/2)^2-(sqrt((points[,1]-0.5)^2+(points[,2]-0.5)^2))^2))/a[cov]+0.2
  for (i in 1:p) {
    beta[i,1:q,cov] = 0.2*sin(seq(0,2*pi, length.out = q))+ beta[i,1,cov]
  }
}


##############Scenario 5######################
beta <- array(NA, dim = c(p,q,m))
a = seq(1.5,4, length.out = m)

for (cov in seq(m)) {
  beta[,1,cov] = (sqrt((sqrt(2)/2)^2-(sqrt((points[,1]-0.5)^2+(points[,2]-0.5)^2))^2))/a[cov]
  for (i in 1:p) {
    beta[i,1:q,cov] = beta[i,1,cov]+log10(seq(1,q, length.out = q))*0.5
  }
}


##############Scenario 6######################
beta <- array(NA, dim = c(p,q,m))
a = seq(1.3,4, length.out = m)

for (cov in seq(m)) {
  beta[,1,cov] = (sqrt((sqrt(2)/2)^2-(sqrt((points[,1]-0.5)^2+(points[,2]-0.5)^2))^2))/a[cov]
  for (i in 1:p) {
    beta[i,2:q,cov] <- 0.015*(seq(2,q, length.out = q-1)-5)^2+beta[i,1,cov]
    beta[i,1,cov]=beta[i,1,cov]+0.015*(1-5)^2
  }
}


##########################################code of simulation#################################################
beta_true = beta

myCluster <- makeCluster(10)
registerDoParallel(myCluster)

results <- foreach(1:nsim) %dopar% {

############data generating##############
  simdat <- list()
  for (i in seq_len(p)) {
    dalocation <- list()
    for (j in seq_len(q)) {
      
      datset <- gdat1(n,m,beta_true= beta[i,j,])
      dalocation[[j]] <- datset[[1]]
    }
    simdat[[i]] <- dalocation
    
  }
  
  lamb1 <- seq(from = 0 , to =200  ,by = 2)
  lamb2 <- seq(from = 0 , to =3  ,by = 0.2)
  lam_cv_err<-matrix(NA, length(lamb1)*length(lamb2), 3)
  
  for (lambda1 in lamb1) {
    for (lambda2 in lamb2) {
      ############cross validation##############
      folds <- cut(seq(1,n),breaks=k,labels=FALSE)
      cv_error<-rep(NA,k)
      
      for (cv in 1:k) {
        
        traindat <- list()
        testdat <-list()
        
        testIndexes <- which(folds==cv,arr.ind=TRUE)
        for (i in seq_len(p)) {
          train <- list()
          test<- list()
          for (j in seq_len(q)) {
            
            test[[j]] <- simdat[[i]][[j]][testIndexes, ] 
            train[[j]] <- simdat[[i]][[j]][-testIndexes, ] 
          }
          traindat[[i]] <- train
          testdat[[i]] <- test
        }
        
        beta0 <- array(0,dim=c(p,q,m))
        beta_cox1 <- estimate.beta.cox(dat = traindat,  p = p, q = q, m = m)
        result = estimate.beta(dat = traindat, p = p, q =q,  beta0 = beta_cox1, lambda1 = lambda1, 
                               lambda2 = lambda2, weight = weight, maxiter = maxiter, epsilon = epsilon)
        
        error<-matrix(NA,length(seq_len(p)),length(seq_len(q)))
        for (i in seq_len(p)) {
          for (j in seq_len(q)) {
            testij <- testdat[[i]][[j]]
            testij<- testij[order(testij$y),]
            deltaij <- as.numeric(testij$failed)
            yij <- as.numeric(testij$y)
            Xij <- as.matrix(subset(testij, select = -c(y, failed)))
            error[i,j]<-cox.logLik(time = yij, status = deltaij, X = Xij, beta = result[i,j,])$logLik
          }
        }
        cv_error1<- sum(error, na.rm = TRUE)
        
        betap <- array(0,dim=c(p,q,m))
        for (a in 1:m) {
          betap[,,a]<-sum((t(diff(t(result[,,a]))))^2)
        }
        p1<-sum(betap[1,1,])/2*lambda1
        
        p2 <- array(NA,dim=c(p*(p-1)/2))
        for (c in 1:(p-1)) {
          for (d in (c+1):p) {
            p2[which(is.na(p2))[1]]<-sum(weight[c,d]*(result[c,,]-result[d,,])^2)
          }
        }
        
        p2<-lambda2*sum(p2)/2    
        cverror<-cv_error1-p1-p2
        
        cv_error[cv]<- sum(cverror, na.rm = TRUE)
      }
      
      lam_cv_err[which(is.na(lam_cv_err))[1],] <- c(lambda1, lambda2, abs(cv_error))
    }
  }
  
  
  ########select the optimal lambdas##########
  
  lambda1<-lam_cv_err[which(lam_cv_err[,3]==min(lam_cv_err[,3],na.rm = TRUE)),1]
  lambda2<-lam_cv_err[which(lam_cv_err[,3]==min(lam_cv_err[,3],na.rm = TRUE)),2]
 
  #######results of the separate method with the optimal lambdas######
  beta_cox = estimate.beta.cox(dat = simdat,  p = p, q = q, m = m)
  
  #######results of the proposed method with the optimal lambdas######
  beta_est = estimate.beta(dat = simdat, p = p, q =q,  beta0 = beta_cox, lambda1 = lambda1, 
                           lambda2 = lambda2, weight = weight, maxiter = maxiter, epsilon = epsilon)
  
  list(lam_cv_err, beta_est, beta_cox, lambda1, lambda2)
  
}

stopCluster(myCluster)

save(results,file = "results.RData")  

############################compute the mean and variance of MSE (mean squared errors) ##########################

########compute the mean of mse#######

####the proposed method
mse1<-rep(NA,nsim)
for (i in 1:nsim) {
  a<-result[[i]][[2]]
  mse1[i]<-(sqrt(sum((a-beta_true)^2)))^2
}
mse_est<-(sum(mse1))/nsim

####the separate method
mse_cox1<-rep(NA,nsim)
for (i in 1:nsim) {
  a<-result[[i]][[3]]
  mse_cox1[i]<-(sqrt(sum((a-beta_true)^2)))^2
}
mse_cox<-(sum(mse_cox1))/nsim


#######compute the variance of MSE#######

####the proposed method
var1<-array(0, dim = c(p,q,m))
for (j in 1:m) {
  for (h in 1:p) {
    for (k in 1:q) {
      b<-rep(NA,nsim)
      for (i in 1:nsim) {
        a=result[[i]][[2]]
        b[i]=a[h,k,j]
        var1[h,k,j]=var(b) 
      }
    }
  }
}
var<-sum(var1)

####the separate method
var1<-array(0, dim = c(p,q,m))

for (j in 1:m) {
  for (h in 1:p) {
    for (k in 1:q) {
      b<-rep(NA,nsim)
      for (i in 1:nsim) {
        a=result[[i]][[3]]
        b[i]=a[h,k,j]
        var1[h,k,j]=var(b) 
      }
    }
  }
}

var_cox<-sum(var1)


