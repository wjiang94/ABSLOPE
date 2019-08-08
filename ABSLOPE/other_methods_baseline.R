library(SSLASSO)
baseline_SSLASSO = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,sigma=1){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  if(p.miss > 0){
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
  
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
  
  # SSLASSO
  list.SSLASSO=SSLASSO(X.sim, y, nlambda = 100, variance = "unknown")
  model.SSLASSO = list.SSLASSO$model
  pick.nlam = which(colSums(abs(list.SSLASSO$beta)>=1e-6) == length(model.SSLASSO))[1]
  beta.SSLASSO = list.SSLASSO$beta[,pick.nlam]
  
  
  pr = power(beta,model.SSLASSO) 
  fdr = fdp(beta,model.SSLASSO) 
  bias_beta = norm(beta.SSLASSO -  beta, "2")/norm(beta,"2")
  #bias_sigma = abs(list.SLOB$sigma - sigma)/sigma
  MSE = norm(X.comp %*% beta.SSLASSO  -  X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}



baseline_SLOB = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,a = 2/p, b = 1-2/p,impute='mean',mec='MCAR',sigma = 1){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  if(p.miss > 0){
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
  
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
 
  objstart<-cv.glmnet(X.sim,y);
  start<-coef(objstart);
  start<-start[2:(p+1)];
  
  list.SLOB = SLOB(start,y,X.sim,lambda,a=a,b=b)
  pr = power(beta,which(list.SLOB$pgamma>1/2)) 
  fdr = fdp(beta,which(list.SLOB$pgamma>1/2)) 
  bias_beta = norm(list.SLOB$beta -  beta, "2")/norm(beta,"2")
  bias_sigma = abs(list.SLOB$sigma - sigma)/sigma
  MSE = norm(X.comp %*% list.SLOB$beta  -  X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}


baseline_SLOPE = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1, impute='mean',mec='MCAR', sigma = 1, sigma.known=NA,estimateReg=TRUE){

  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta

  if(p.miss > 0){
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
  
  if(is.na(sigma.known)){
    list.SLOPE = SLOPE(X.sim,y,fdr = 0.10)
  } else{list.SLOPE = SLOPE(X.sim,y,fdr = 0.10,sigma=sigma.known)}
  
  if(estimateReg){
    beta.est = rep(0,p)
    modelsel = list.SLOPE$selected
    if(length(modelsel)>0){
      beta.est[modelsel]= summary(lm(y ~ X.sim[,modelsel]))$coef[-1,1]
    }
  } else{
    beta.est = list.SLOPE$beta
  }
  
  pr = power(beta,list.SLOPE$selected) 
  fdr = fdp(beta,list.SLOPE$selected) 
  bias_beta = norm(beta.est -  beta, "2")/norm(beta, "2")
  bias_sigma = abs(list.SLOPE$sigma - sigma)/sigma
  MSE = norm(X.comp %*% beta.est  -  X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}

baseline_LASSO.cv = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,impute='mean',mec='MCAR', sigma = 1){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  if(p.miss > 0){
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
  
  cvfit <- glmnet::cv.glmnet(X.sim, y)
  # lambda.1se s.t. error is within 1 standard error of minimum error
  # to choose the simplest model whose accuracy is comparable with the best model.
  beta_lasso = coef(cvfit, s = "lambda.1se")#?coef.cv.glmnet
  beta_lasso=beta_lasso[2:(p+1)]
  
  pr = power(beta,which(beta_lasso!=0)) 
  fdr = fdp(beta,which(beta_lasso!=0)) 
  bias_beta = norm(beta_lasso -  beta, "2")/norm(beta,"2")
  MSE = norm(X.comp %*% beta_lasso - X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  bias_sigma = NA
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}



baseline_ncLASSO = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,mec='MCAR',sigma=1,oracle=TRUE){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  XNA = data.list$X.obs
  X0 = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  ## NClasso (nonconvex, Loh and Wainwright 2012)
  MAXITS <- 5000
  stepsize <- 0.05
  
  if(oracle){
    R <- sqrt(sum(beta^2)) * sqrt(nspr) # oracle
    maxSparsity = nspr
  } else {
    R  =  signallevel*p
    maxSparsity = p}
  lambda.max <- 1.2*max(abs(t(X0) %*% y))/n
  lambda.seq <- seq(lambda.max, 0.001*lambda.max, length=200)
  NClasso.res <- NClasso(y, XNA,  MAXITS, stepsize, R, lambda.seq, maxSparsity = maxSparsity)$betas
  lambda.ix <- min(which(apply(NClasso.res != 0, 2, sum) >= maxSparsity)) #oracle
  #non-oracle
  NClasso.res <- NClasso.res[, lambda.ix]
  
  pr = power(beta,which(NClasso.res!=0)) 
  fdr = fdp(beta,which(NClasso.res!=0)) 
  bias_beta = norm(NClasso.res -  beta, "2")/norm(beta,"2")
  MSE = norm(X0 %*% NClasso.res  - X0 %*% beta, "2")/norm(X0 %*% beta,"2")
  bias_sigma = NA
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}


baseline_adaLASSO = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,impute='mean',mec='MCAR', sigma = 1){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  if(p.miss > 0){
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
  
  
  ###LASSO to get initial estimate
  A=cv.glmnet(X.sim,y,intercept=FALSE,standardize=FALSE, thresh=10^(-12))
  # A$lambda
  # A$lambda.min
  # A$glmnet.fit$beta
  # coef(A,s=A$lambda.min)

  SL=coef(A,s=A$lambda.min)[2:(p+1)]; #SL=coefficients(A)[2:(p+1),];
  ###Adaptive LASSO 
  e=0.0000001
  
  W=abs(SL)+e   #W=abs(SL[,2])+e
  X1=X.sim%*%diag(W)
  lambda = create_lambda_bhq(ncol(X.sim),fdr=0.10)
  #### Solution adaptive lasso, Here lambda[1] is the first (largest) lambda of SLOPE
  A=glmnet(X1,y,intercept=FALSE,standardize=FALSE,lambda=lambda[1]/n,thresh=10^(-12))
  beta_lasso=coefficients(A)[2:(p+1)] * W
  
  pr = power(beta,which(beta_lasso!=0)) 
  fdr = fdp(beta,which(beta_lasso!=0)) 
  bias_beta = norm(beta_lasso -  beta, "2")/norm(beta,"2")
  MSE = norm(X.comp %*% beta_lasso  - X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  bias_sigma = NA
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}





# baseline_ensemble = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,mec='MCAR'){
#   set.seed(nb.seed)
#   print(nb.seed)
#   data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec) 
#   X = data.list$X.obs
#   y = data.list$y
#   beta = data.list$beta
#   res.algo<-algo(nnodes=detectCores()-1,#parallelisation
#                  X,#matrice des variables explicatives
#                  y,#variable reponse (compl?te)
#                  #path.outfile="experiment_results",#chemin pour exporter resultats d'affichage
#                  Nb.sim=1000,#parametre B dans article
#                  methods="lasso",#methode de selection (autre choix possibles : lasso et knockoff)
#                  proport=nspr,#parametre nspr dans article
#                  printflag=FALSE,#un peu d'affichage
#                  method.na="norm"#imputation par modele gaussien
#   )
#   pr = power(beta,res.algo$res$selection[[1]]) 
#   fdr = fdp(beta,res.algo$res$selection[[1]]) 
#   bias = norm(colMeans(res.algo$res$beta[[1]]) -  beta, "2")
#   return(list(pr=pr,fdr=fdr,bias=bias))
# }
