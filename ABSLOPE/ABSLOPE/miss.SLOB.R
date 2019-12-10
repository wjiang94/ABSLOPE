library(SLOPE)
library(truncdist)
library(nlshrink)
library(MASS)
library(glmnet)
library(missMDA)

miss.SLOB<-function(X, y, lambda, a, b, beta.start = NA, maxit = 300, case = 'MCAR', seed = NA, print_iter = FALSE, tol_em=1e-6, impute = 'PCA', sigma.known=NA, sigma.init=NA, method_na = 'MH'){
  
  # missing pattern
  rindic = as.matrix(is.na(X)) 
  if(sum(rindic) > 0){ # missing data exist
    whichcolmissing = (1:ncol(rindic))[apply(rindic,2,sum)>0] 
    missingcols = length(whichcolmissing) 
  }
  if(sum(rindic) == 0){missingcols = 0} # no missingness
  
  # delete rows completely missing
  if(missingcols != 0){
    if(any(apply(is.na(X),1,sum) == p)){
      i_allNA = which(apply(is.na(X),1,sum) == p)
      X = X[-i_allNA,]
      y = y[-i_allNA]
    }
    if(any((is.na(y)) == TRUE)){
      i_YNA = which(is.na(y) == TRUE)
      X = X[-i_YNA,]
      y = y[-i_YNA]
    }
  }
  
  p = ncol(X)
  n = length(y)
  
  if(missingcols != 0){
    # impute by PCA 
    if(impute == 'PCA'){
      X.sim = imputePCA(X)$completeObs
    } else{ # impute by mean
      X.mean = X
      for(i in 1:ncol(X.mean)){
        X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)
      }
      X.sim <- X.mean
    }
  } else{X.sim = X} # no missingness
  
  ##------------------------------------
  # initialization
  eps<-1
  
  # beta
  if(is.na(beta.start)){ # initialization not given
    # LASSO cv
    objstart = cv.glmnet(X.sim,y)
    beta = coef(objstart,s = "lambda.1se")
    beta = beta[2:(p+1)]
  } else{beta = beta.start} # initialization given
  
  beta_old <- beta
  
  # sigma
  if(is.na(sigma.known)){ # real value of sigma is unknown
    if(is.na(sigma.init)){ # initialization not given
      #sigma = sd(y - X.sim %*% beta) 
      sigma = sqrt(sum((y - X.sim %*% beta)^2)/(n-1))
    } else{sigma = sigma.init} # initialization given
  } else{sigma = sigma.known} # real value of sigma is known
  

  #rank
  rk <- p-rank(abs(beta),ties.method="max")+1;
  
  lambda_sigma=lambda * sigma
  
  lambda_sigma_inv <- lambda /sigma
  
  # c
  c_k <- min((sum(abs(beta)>0)+1)/sum(abs(beta[beta != 0]))/lambda_sigma_inv[p],1)
  if(!is.finite(c_k)) c_k <- 1
  
  # theta
  theta<-(sum(beta!=0)+a)/(p+b+a); 
  
  # pi
  pstar_k<-1/(1+(1-theta)/theta/c_k*exp(-abs(beta)*lambda_sigma_inv[rk]*(1-c_k)));
  
  
  k<-0
  
  
  seqsigma = seqeps= rep(NA,(maxit+1))
  while(eps>tol_em & k<maxit | k<20){ 
    k<- k+1
    if(missingcols != 0){
      # mu and Sigma
      mu = apply(X.sim,2,mean)
      #Sigma = var(X.sim)*(n-1)/n
      Sigma = linshrink_cov(X.sim)
    }
    
    # weigths
    wts<-pstar_k*c_k+(1-pstar_k);
    revwts<-1/wts;
    revwts<-as.vector(revwts);
    Xtemp<-sweep(X.sim,2,revwts,'*');
    utemp<-slope_admm(Xtemp, y, rep(0, p), rep(0, p), lambda_seq = lambda_sigma, 1)
    utemp1<-utemp$z;
    
    # beta
    beta<-revwts*utemp1;
    
    # rank
    rk<-p-rank(abs(utemp1),ties.method="max")+1;
    
    # sigma
    RSS<- sum((y-X.sim%*%beta)^2)
    sum_lamwbeta <- sum(lambda[rk]*abs(utemp1))
    sigma <- (sqrt(sum_lamwbeta^2+4*n*RSS)+sum_lamwbeta)/(2*n)
    
    lambda_sigma<-lambda * sigma
    lambda_sigma_inv <- lambda / sigma
    
    
    
    #########
    # Xmis
    
    if(missingcols != 0){
      if(method_na == 'MH'){
        S.inv = solve(Sigma)
        for (i in (1:n)) {
          yi = y[i]
          jna <- which(is.na(X[i,]))
          njna <- length(jna)
          if (njna>0) {
            xi <- X.sim[i,]
            Oi <- Sigma[jna,jna]
            mi <- mu[jna]
            if (njna<p) {
              jobs <- setdiff(1:p,jna)
              mi <- mi + Sigma[jna,jobs] %*% solve(Sigma[jobs,jobs]) %*% (xi[jobs] - mu[jobs])
              Oi <- Oi - Sigma[jna,jobs] %*% solve(Sigma[jobs,jobs]) %*% Sigma[jobs,jna]
            }
            nmcmc = 100
            avg.xina = 0
            xina  <- xi[jna]
            for (m in (1:nmcmc)){
              xina.c <- mvrnorm(n = 1, mu=mi, Sigma =Oi)
              xi[jna] = xina.c
              alpha = dnorm(yi, mean=xi%*%beta , sd=sigma, log = FALSE)/dnorm(yi, mean=xi[jobs]%*%beta[jobs] + xina%*%beta[jna] , sd=sigma, log = FALSE)
              if (runif(1) < alpha){
                xina <- xina.c
              }
              avg.xina = (avg.xina * (m-1) + xina)/m
            }
            X.sim[i,jna] <- avg.xina
          }
        }
      }
      if(method_na == 'lineq'){
        S.inv = solve(Sigma)
        m = S.inv %*% mu
        tau = sqrt(diag(S.inv) + (beta/sigma)^2)
        for (i in (1:n)) {
          yi = y[i]
          jna <- which(is.na(X[i,]))
          njna <- length(jna)
          if (njna>0) {
            xo <- X[i,-jna]
            betai = beta[jna]
            mi = m[jna]
            ui = S.inv[jna,-jna] %*% xo
            r = (yi - xo %*% beta[-jna])[1,1]
            taui = tau[jna]
            
            # linear equation Ax = cc 
            cc = (r * betai / sigma^2 + mi - ui)/taui
            A = (betai %*% t(betai) / sigma^2 + S.inv[jna,jna]) / (taui %*% t(taui))
            diag(A) <- 1
            
            #showEqn(A, cc)
            mu_tilde = solve(A, cc,tol=1e-20)
            
            X.sim[i,jna] <- mu_tilde / taui
          }
        }
      }
    }
    #########
    
    
    pstar_k<-1/(1+(1-theta)/theta/c_k*exp(-abs(beta)*lambda_sigma_inv[rk]*(1-c_k)))
    
    
    rate_temp<-sum(abs(beta)*lambda_sigma_inv[rk]*pstar_k)
    
    c_k<-min((sum(pstar_k)+1)/rate_temp,1)
    
    if(!is.finite(c_k)) c_k<-1
    
    theta<-(sum(pstar_k)+a)/(p+b+a); 
    
    
    eps<-sum((beta_old-beta)^2)
    seqeps[k] = eps
    seqsigma[k] = sigma
    beta_old<-beta
    if ((print_iter == TRUE) & (k%%100==0)){
      cat(sprintf('iteration = %i ', k))
      cat(sprintf('Distance from last iter ='),eps,'\n')
      #cat(sprintf(' c ='),c,'\n')
      #cat(sprintf(' pi.binom ='),pi.binom)
      #cat(sprintf(' gamma ='),gamma, '\n')
    }
  }
  
  result<-list(beta=beta,pgamma=pstar_k, sigma=sigma, c=c_k, theta=theta,seqeps=seqeps,seqsigma=seqsigma)
  return(result)
}

# Baseline of miss.SLOPE
baseline_miss.SLOB= function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300,impute='mean',mec='MCAR',sigma = 1,beta.start =NA,sigma.known=NA,sigma.init=NA,print_iter=FALSE,method_na='lineq'){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
  #start_time <- Sys.time()
  list.SLOB = miss.SLOB(X, y, lambda, a=a, b=b, beta.start=beta.start,maxit = maxit,print_iter =print_iter,tol_em=tol_em,impute=impute,sigma.known=sigma.known,sigma.init=sigma.init,method_na=method_na)
  #time <- Sys.time() - start_time
  pr = power(beta,which(list.SLOB$pgamma>1/2)) 
  fdr = fdp(beta,which(list.SLOB$pgamma>1/2)) 
  bias_beta_new =norm(list.SLOB$beta  - beta, "2")/norm(beta,"2")
  bias_sigma = abs(list.SLOB$sigma - sigma)/sigma
  MSE_new = norm(X.comp %*% list.SLOB$beta  -  X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta_new, bias_sigma = bias_sigma,MSE=MSE_new))
  
}



