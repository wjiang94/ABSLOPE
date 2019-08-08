power <- function(beta,selected){sum(selected %in% which(beta!=0)) / max(1, length(which(beta!=0)))}
fdp <- function(beta,selected){sum(beta[selected] == 0) / max(1, length(selected))}
TP <- function(beta,gamma){sum(beta != 0 & gamma == 1)}
FP <- function(beta,gamma){sum(beta == 0 & gamma != 0)}
FN <- function(beta,gamma){sum(beta != 0 & gamma == 0)}
TN <- function(beta,gamma){sum(beta == 0 & gamma == 0)}

# Generate dataset with missing values (case MCAR)
library(mice)
data.generation <- function(n=100, p=100, nspr=10, p.miss=0.1, mu = rep(0,p), Sigma=toeplitz(0^(0:(p-1))), signallevel=3, mec='MCAR',sigma=1,scale=TRUE){

  amplitude = signallevel*sqrt(2*log(p))  # signal amplitude (for noise level = 1)

  # Design matrix
  # normal distribution
  X <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
  if(scale){
    X = scale(X)/sqrt(n)
  }
  
  # Coefficient and response vectors
  nonzero = sample(p, nspr)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta +  sigma*rnorm(n)
  
  # Missing values
  X.obs <- X

  if(p.miss>0){
    if(mec=='MCAR'){
    patterns <- runif(n*p)< p.miss # missing completely at random
    X.obs[patterns] <- NA
    }
    if(mec=='MAR'){
      nb.pat = 5
      patterns = 1 - matrix(rbinom(p*nb.pat, 1, p.miss*1.5), ncol = p, nrow = nb.pat)
      list.amp <- ampute(X, prop = p.miss, mech = 'MAR',bycases = FALSE, patterns = patterns)
      X.obs <- as.matrix(list.amp$amp)
      # library(VIM)
      # matrixplot(as.data.frame(X.obs), interactive = T,axes =TRUE, xlab ='index of variables' ,ylab='index of observations')
    }
  }
  return(list(X.obs=X.obs, X=X, y=y, beta=beta,c = 1/ amplitude))
}
# data.list = data.generation()
# X = data.list$X.obs
# y = data.list$y

#-------------------------------------------------------------------
#BHQ lambda 
create_lambda_bhq <- function(p, fdr) {
  q = (1:p) * fdr / (2*p)
  qnorm(1 - q)
}

#Gaussian lambda
create_lambda_gaussian <- function(n, p, fdr) {
  w <- function(nspr) 1 / max(1, n - nspr - 1)
  lambda.bhq = create_lambda_bhq( p, fdr)
  lambda = rep(0,p)
  lambda[1] <- lambda.bhq[1]
  if (p >= 2) {
    sum_sq <- 0
    for (i in 2:p) {
      sum_sq <- sum_sq + lambda[i-1]^2
      lambda[i] <- lambda.bhq[i] * sqrt(1 + w(i-1) * sum_sq)
    }
  }
  return(lambda)
}
#lambda = create_lambda_bhq(ncol(X),fdr=0.1)
#lambda = create_lambda_gaussian(nrow(X),ncol(X),fdr=0.1)
#-------------------------------------------------------------------
library(SLOPE)

slope_admm <- function(A, b, z, u, lambda_seq, rho, 
                       max_iter = 500, tol_infeas = 1e-8,
                       verbose = FALSE){ 
  M <- solve(crossprod(A) + diag(rho, ncol(A)),tol=1e-22)
  MtAb <- M %*% crossprod(A,b)
  lambda_seq_rho <- lambda_seq/rho
  z_new <- NULL
  for(iter in 1:max_iter){ #just until we do not choose some reasonable convergence criterion
    
    x <- MtAb + crossprod(M, (rho*(z - u)))
    z_new <- SLOPE::prox_sorted_L1(x = as.vector(x + u), lambda = lambda_seq_rho,method = 'c')
    u <- u + x - z_new
    
    dual_feasibility <- sqrt(sum((rho*(z_new-z))^2))
    primal_feasibility <- sqrt(sum((z_new - x)^2))
    
    z <- z_new
    
    if(verbose)
      message(sprintf("Iter %i\nprimal: %f\ndual: %f\n", 
                      iter, primal_feasibility, dual_feasibility))
    
    if(dual_feasibility < tol_infeas & primal_feasibility < tol_infeas){
      break;
    }
  }
  return(list(x = x, z = z, u = u, 
              primal_feasibility = primal_feasibility, 
              dual_feasibility = dual_feasibility))
}



rowMedian <- function(x, na.rm = FALSE){
apply(x, 1, median, na.rm = na.rm) 
}
colMedian <- function(x, na.rm = FALSE){
  apply(x, 2, median, na.rm = na.rm) 
}

