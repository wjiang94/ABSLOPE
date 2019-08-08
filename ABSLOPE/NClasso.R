NClasso <- function(y, XNA,  MAXITS, stepsize, R, lambdaseq, maxSparsity) {
  
  # Non convex Lasso for missing data (Loh, Wainwright, 2012)
  # 
  # y: response vector
  # XNA: incomplete matrix
  # MAXITS, stepsize: parameters for composite gradient descent
  # R: bound for l1-norm of beta
  # lambdaseq: sequence of tuning parameters
  # 
  # Compute the solution path of the Lagrangian version
  # for a sequence of lambda values
  
  
  ######################## functions of Po-Ling Loh ##########################
  ## (my only modification: doCompGrad: add initiat beta value as argument) ##
  project_onto_l1_ball <- function(v, b){
    u = sort(abs(v), decreasing = TRUE);
    sv = cumsum(u);
    rho = tail(which(u > (sv - b) / c(1:length(u))), n=1);
    theta = pmax(0, (sv[rho] - b) / rho);
    return(sign(v) * pmax(abs(v) - theta, 0));
  }
  
  softThresh <- function(eta, lambda){
    # thresholding operator used for composite gradient descent
    p = length(eta);
    beta = c(rep(0, p));
    for(i in 1:p){
      if(eta[i] > lambda){
        beta[i] = eta[i] - lambda;
      } else if(eta[i] < -lambda){
        beta[i] = eta[i] + lambda;
      } else
        beta[i] = 0;
    }
    return(beta);
  }
  
  doCompGrad <- function(Gamma, gamma, MAXITS, stepsize, R, lambda, betastart){
    # composite gradient descent for Lagrangian version of l_1-constrained QP
    # min {1/2 * beta^T Gamma beta - gamma^T beta + lambda * ||beta||_1} s.t. ||beta||_1 <= R
    p = length(gamma);
    #betaold = c(rnorm(p));
    betaold <- betastart
    betaold = betaold / sqrt(sum(betaold^2))
    #norm(as.matrix(betaold), "2");
    TOL = 1e-4;
    its = 1;
    change = 1;
    while((its <= MAXITS) && (change > TOL)){
      grad = Gamma %*% betaold - gamma;
      betanew = softThresh(betaold - stepsize * grad, lambda * stepsize);
      if(norm(as.matrix(betanew), "1") > R){
        betanew = project_onto_l1_ball(betanew, R);
      }
      diff = sqrt(sum((betanew - betaold)^2))
      #norm(as.matrix(betanew - betaold), "2");
      change = diff / min(stepsize, 0.1);
      betaold = betanew;
      its = its + 1;
    }
    if(its == MAXITS + 1){
      sprintf('max iterations');
    }
    return(list(betanew, its));
  }
  ###########################################################################
  
  ## compute surrogate matrix and vector:
  
  n <- nrow(XNA)
  p <- ncol(XNA)
  Nobs <- apply(!is.na(XNA), 2, sum) # number of observed values in each column
  Z <- XNA
  Z[is.na(XNA)] <- 0
  # compute gamma:
  gamma <- t(Z) %*% y / Nobs
  # compute Gamma:
  Z <- t(t(Z) * sqrt(n) / Nobs)
  Gamma <- t(Z) %*% Z
  # correction on diagonal:
  diag(Gamma) <- diag(Gamma) * Nobs / n
  
  ## compute sol for each lambda
  lambdaseq <- sort(lambdaseq, decreasing = TRUE)
  betas <- matrix(0, p, length(lambdaseq))
  betastart <- rnorm(p, 0, 1/sqrt(n))
  k <- 0
  cont <- TRUE
  while (cont) {
    k <- k+1
    lambda <- lambdaseq[k]
    betastart <- doCompGrad(Gamma, gamma, MAXITS, stepsize, R, lambda, betastart)[[1]]
    betas[, k] <- betastart
    cont <- (sum(betas[, k] != 0) < maxSparsity) & (k < length(lambdaseq))
    if(sum(betastart^2) < 1.e-4) betastart <- rnorm(p, 0, 1/sqrt(n)) # her code starts with nonzero inital value for beta
  }
  if (k < length(lambdaseq)) {
    lambdaseq <- lambdaseq[1:k] 
    betas <- betas[, 1:k]
  }
  
  out <- NULL
  out$betas <- betas
  out$lambdaseq <- lambdaseq
  out$gamma <- gamma
  out$Gamma <- Gamma
  out$Nobs <- Nobs
  return(out)
  
}



prox_l1 <- function(x, tau) {
  #computes the projection on the set M.*X = X_na
  res = max(0, 1 - tau / max(1e-15, abs(x))) * x
  return(res)
}


soft_thresh_svd <- function(X, gamma) {
  #performs soft thresholding in the value decomposition
  res=svd(X)
  U = res$u
  V = res$v
  S = diag(res$d)
  s = sapply(diag(S), prox_l1, gamma)
  S[1:length(s), 1:length(s)] = diag(s)
  X_res = U %*% S %*% t(V)
  max_sing=max(s)
  nuc_norm_res=sum(s)
  return(list(X_res=X_res,nuc_norm_res=nuc_norm_res,max_sing=max_sing))
}

min_FISTA_nuc <- function(Y,mask,niter,lambda){
  X=matrix(0,ncol=ncol(mask),nrow=nrow(mask))
  W=matrix(0,ncol=ncol(mask),nrow=nrow(mask))
  L=1
  for (i in 1:niter){
    theta=2/(i+1)
    X_=X
    XX = (1-theta)* X + theta*W
    X = soft_thresh_svd(XX-(1/L)*mask*(XX-Y),1/L*lambda)$X_res
    W = X_ + (X-X_)/theta
  }
  return(X)
}

min_FISTA_l1 <- function(y, X, niter, lambda) {
  beta = rep(0, dim(X)[2])
  v = rep(0, dim(X)[2])
  resSVD = svd(X)
  U = resSVD$u
  S = resSVD$d
  V = resSVD$v
  L = max(S) ^ 2
  for (i in 1:niter) {
    theta = 2 / (i + 1)
    beta0 = beta
    xx = (1 - theta) * beta + theta * v
    beta = prox_l1(xx - (1 / L) * t(X) %*% (X %*% xx - y), 1 / L * lambda)
    v = beta0 + (beta - beta0) / theta
  }
  return(beta)
}


supp_identifier <- function(beta,s){
  vec = sort(abs(beta),index.return=TRUE,decreasing=TRUE)
  betaNew=rep(0,n)
  betaNew[vec$ix[1:s]] = beta[vec$ix[1:s]]
  res=betaNew
  return(res)
}


test_alternate_min_palm <-
  function(X_na, y, M, beta0, s, Xnoisy) {
    ## 2-step procedure and also initialization of alternate min
    lambda_X_vec = seq(1e-4,1e-2,length=10)
    errX = c()
    for (i in 1:length(lambda_X_vec)){
      X_ = min_FISTA_nuc(X_na, M, 50, lambda_X_vec[i])
      errX = c(errX, norm(X_-X0,"F"))
    }
    X_ = min_FISTA_nuc(X_na, M, 50, lambda_X_vec[which.min(errX)])
    lambda_vec = seq(1,100,length=10)
    err=c()
    for (i in 1:length(lambda_vec)){
      beta_ = min_FISTA_l1(y, X_, 100, lambda_vec[i])
      err = c(err,sqrt(sum((beta_ - beta0) ^ 2)))
    }
    beta_ = min_FISTA_l1(y, X_, 100, lambda_vec[which.min(err)])
    
    # Brute force parameters tuning
    lambda_beta_vec = c(1e-3, 2e-2, 1e-1, 1e1, 1e1, 1e2, 1e3 )
    lambda_X_vec = c(1e-3, 2e-2, 1e-1, 1e1, 1e1, 1e2, 1e3 )
    lambda3_vec = c(1e-3, 2e-2, 1e-1, 1e1, 1e1, 1e2, 1e3 )
    
    comb3 <- function(...)
      abind(..., along = 3)
    comb4 <- function(...)
      abind(..., along = 4)
    
    result <-
      foreach (i = 1:length(lambda3_vec), .combine = "comb4")  %:%
      foreach (j = 1:length(lambda_X_vec), .combine = "comb3")  %:%
      foreach (k = 1:length(lambda_beta_vec), .combine = "cbind")  %dopar% {
        
        lambda3 = lambda3_vec[i]
        #Initialisation
        X = X_
        beta = beta_
        tau = 1 / (lambda3 + sum(beta ^ 2))
        sigma = 1 / norm(t(X) %*% X)
        lambda_X = lambda_X_vec[j]
        lambda_beta = lambda_beta_vec[k]
        flag_X = 0
        Tt =200
        for (t in 1:Tt) {
          Xaux = X - tau * (lambda3 * M * (X - Xnoisy) + (X %*% beta - y) %*% t(beta))
          if (sum(is.na(Xaux)) > 0) {
            flag_X = 1
            break
          }
          if (sum(is.infinite(Xaux)) > 0) {
            flag_X = 1
            break
          }
          res = soft_thresh_svd(Xaux, tau * lambda_X)
          X = res$X_res
          nucX = res$nuc_norm_res
          max_sing_X = res$max_sing
          beta = prox_l1(beta - sigma * t(X) %*% (X %*% beta - y), sigma *
                           lambda_beta)
          tau = 1 / (lambda3 + 2 * sum(beta ^ 2))
          sigma = 1 / max_sing_X ^ 2
          if (is.na(tau)) {
            flag_X = 1
          }
        }
        if (flag_X == 1) {
          c(Inf,-1)
        } else{
          rbind(sqrt(sum((supp_identifier(beta,s) - beta0) ^ 2)),sum(beta0 != 0 & supp_identifier(beta,s) != 0) / s)
        }
      }
    abc <-
      array(dim = c(
        length(lambda_beta_vec),
        length(lambda_X_vec),
        length(lambda3_vec)
      ))
    resultvec <- as.vector(result[1,1:length(lambda_beta_vec),1:length(lambda_X_vec),1:length(lambda3_vec)])
    resultvecsupp <- as.vector(result[2,1:length(lambda_beta_vec),1:length(lambda_X_vec),1:length(lambda3_vec)])
    
    ind <- which.min(resultvec[resultvecsupp==max(resultvecsupp)])
    resind <- arrayInd(which.min(resultvec), dim(abc))
    
    Tt = 200
    X = X_
    lambda3 = lambda3_vec[resind[3]]
    beta = beta_
    tau = 1 / (lambda3 + sum(beta ^ 2))
    sigma = 1 / norm(t(X) %*% X)
    lambda_X = lambda_X_vec[resind[2]]
    lambda_beta = lambda_beta_vec[resind[1]]
    flag_X = 0
    for (t in 1:Tt) {
      Xaux = X - tau * (lambda3 * M * (X - Xnoisy) + (X %*% beta - y) %*% t(beta))
      if (sum(is.na(Xaux)) > 0) {
        flag_X = 1
        break
      }
      if (sum(is.infinite(Xaux)) > 0) {
        flag_X = 1
        break
      }
      res = soft_thresh_svd(Xaux, tau * lambda_X)
      X = res$X_res
      nucX = res$nuc_norm_res
      max_sing_X = res$max_sing
      beta = prox_l1(beta - sigma * t(X) %*% (X %*% beta - y), sigma *
                       lambda_beta)
      tau = 1 / (lambda3 + 2 * sum(beta ^ 2))
      sigma = 1 / max_sing_X ^ 2
      if (is.na(tau)) {
        flag_X = 1
      }
    }
    
    if (flag_X == 1) {
      beta = c()
      X = c()
    }
    
    return(list(beta_al_min = beta, X_al_min = X, param = resind, res_grilleParam = result))
  }