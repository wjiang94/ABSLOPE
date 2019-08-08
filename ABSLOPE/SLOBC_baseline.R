library(Rcpp)
Rcpp::sourceCpp('ABSLOPE/SLOBE_cpp_missing.cpp')

baseline_SLOBC = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,a = 2/p, b = 1-2/p,impute='mean',mec='MCAR',sigma = 1){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  X.mean = X
  for(i in 1:ncol(X.mean)){
    X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)
  }
  X.sim <- X.mean
  
  cv.lasso<-cv.glmnet(X.sim,y) 
  start = as.numeric(coef(cv.lasso))[-1]
  
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
  
  objstart<-cv.glmnet(X.sim,y);
  start<-coef(objstart);
  start<-start[2:(p+1)];
  
  list.SLOB = SLOBE_ADMM_approx_missing(start,X,y,a_prior = 1, b_prior = 10,FDR=0.1,verbose = F  )
  pr = power(beta,which(list.SLOB$beta!=0)) 
  fdr = fdp(beta,which(list.SLOB$beta!=0)) 
  bias_beta = norm(list.SLOB$beta -  beta, "2")/norm(beta,"2")
  bias_sigma = NA
  MSE = norm(X.comp %*% list.SLOB$beta  -  X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta, bias_sigma = bias_sigma,MSE=MSE))
}