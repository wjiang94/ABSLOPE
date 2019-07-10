#' Data generation
#'
#' Generate dataset with missing values
#' @param n Observation number
#' @param p Dimenaion
#' @param nspr Non-sparsity
#' @param p.miss Percentage of missingness
#' @param mu Mean of covariates
#' @param Sigma Covariance of covariates
#' @param signallevel Signal magnitude
#' @param mec Missing mechanism
#' @param sigma noise variance
#' @param scale If TRUE, perform scaling after generating dataset.
#' @import mice stats
#' @examples
#' data.list = data.generation()
#' X = data.list$X.obs
#' y = data.list$y
#' @export
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



#' Create BHq sequence
#'
#' Create BHq sequence
#' @param p Dimenaion
#' @param fdr Expected False Discovery Rate
#' @import stats
#' @examples
#' create_lambda_bhq(50, 0.1)
#' @export
create_lambda_bhq <- function(p, fdr) {
  q = (1:p) * fdr / (2*p)
  qnorm(1 - q)
}

#' Optimization solver
#'
#' Solve optimation with proximal gradient descent
#' @param A Weighted design matrix
#' @param b Response
#' @param z Initailization
#' @param u Initailization
#' @param lambda_seq Penalization coefficient
#' @param rho Step-size
#' @param max_iter Maximum iteration number
#' @param tol_infeas Tolerance
#' @param verbose If TRUE, print results in each iteration
#' @import SLOPE
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

