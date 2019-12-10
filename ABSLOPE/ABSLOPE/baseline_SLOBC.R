# Baseline of ABSLOPEC
library(MASS)
library(glmnet)
library(Rcpp)
Rcpp::sourceCpp('ABSLOPE/SLOBE_cpp_missing.cpp')

baseline_SLOBC = function(n=100, p=20, nspr=6, p.miss=0.1, mu = rep(0,p), Sigma = diag(p), signallevel=3, nb.seed=1,a = 2/p, b = 1-2/p,mec='MCAR',sigma = 1,max_iter=100,a_prior = 1, b_prior = 10,FDR=0.1){
  set.seed(nb.seed)
  print(nb.seed)
  data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel,mec=mec,sigma=sigma) 
  X = data.list$X.obs
  X.comp = data.list$X
  y = data.list$y
  beta = data.list$beta
  
  Xim = X
  impute_mean(Xim,n,p)
  cv.lasso<-cv.glmnet(Xim,y) 
  start = as.numeric(coef(cv.lasso))[-1]
  list.SLOBC = SLOBE_ADMM_approx_missing(start,X,y,a_prior = a_prior, b_prior = b_prior,FDR=FDR,verbose = FALSE,max_iter = max_iter)
  pr = power(beta,which(list.SLOBC$beta!=0)) 
  fdr = fdp(beta,which(list.SLOBC$beta!=0)) 
  bias_beta_new =norm(list.SLOBC$beta  - beta, "2")/norm(beta,"2")
  bias_sigma = abs(list.SLOBC$sigma - sigma)/sigma
  MSE_new = norm(X.comp %*% list.SLOBC$beta  - X.comp %*% beta, "2")/norm(X.comp %*% beta,"2")
  
  return(list(pr=pr,fdr=fdr,bias_beta=bias_beta_new, bias_sigma = bias_sigma,MSE=MSE_new))
}


