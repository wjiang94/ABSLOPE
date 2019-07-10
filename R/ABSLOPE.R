#' ABSLOPE
#'
#' This function uses algorithm SAEM to fit the ABSLOPE model with missing data.
#' @param X Design matrix with missingness \eqn{n \times p}{n * p}
#' @param y Response vector \eqn{n \times 1}{n * 1}
#' @param lambda Penalization coefficient
#' @param a Coefficient for sparsity prior
#' @param b Coefficient for sparsity prior
#' @param beta.start Initialization for \eqn{\beta}{\tau}
#' @param maxit Maximum number of iteration. The default is maxit = 300.
#' @param case Missing mechanism
#' @param seed An integer as a seed set for the radom generator.
#' @param tol_em The tolerance to stop SAEM. The default is tol_em = 1e-6.
#' @param impute Imputation method for the initilization. The default is impute = 'PCA'.
#' @param print_iter If TRUE, miss.saem will print the estimated parameters in each iteration of SAEM.
#' @param sigma.known If TRUE, \eqn{\sigma}{\sigma} is known.
#' @param sigma.init Value for known \eqn{\sigma}{\sigma}.
#' @param scale If TRUE, perform scaling in each iteration.
#' @param method_na Method for dealing with missing values. Default is method_na = 'lineq'.
#' @return A list with components
#' \item{X.sim}{Simulated dataset in the last iteration.}
#' \item{beta}{Estiamated \eqn{\beta}{\beta}.}
#' \item{seqbeta}{Sequence of \eqn{\beta}{\beta} estimated in each iteration.}
#' \item{beta.new}{Estiamated \eqn{\beta}{\beta} * model selection indicator \eqn{\gamma}{\gamma}.}
#' \item{gamma}{Simulated model selection indicator \eqn{\gamma}{\gamma} in the last iteration.}
#' \item{gamma.avg}{Average \eqn{\gamma}{\gamma} during the iterations.}
#' \item{seqgamma}{Sequence of \eqn{\gamma}{\gamma} in each iteration.}
#' \item{pi.binom}{Probability for selecting each variable.}
#' \item{theta}{Simulated sparsity \eqn{\theta}{\theta} in the last iteration.}
#' \item{seqtheta}{Sequence of \eqn{\theta}{\theta} in each iteration.}
#' \item{sigma}{Estimated \eqn{\sigma}{\sigma}.}
#' \item{seqsigma}{Sequence of \eqn{\sigma}{\sigma} estimated in each iteration.}
#' \item{mu}{Estimated \eqn{\mu}{\mu}.}
#' \item{Sigma}{Estimated \eqn{\Sigma}{\Sigma}.}
#' \item{c}{Simulated inverse of average signal magnitude c in the last iteration.}
#' \item{seqc}{Sequence of c in each iteration.}
#' @import truncdist nlshrink MASS glmnet missMDA mice stats
#' @examples
#' # Generate dataset
#' n=50
#' p=50
#' nspr=10 # non-sparsity
#' p.miss=0.1 # percentage of missingness
#' mu = rep(0,p)
#' Sigma=toeplitz(0^(0:(p-1)))
#' signallevel=3 # signal strength
#' sigma=1 # noise
#' amplitude = signallevel*sqrt(2*log(p))  # signal amplitude
#' # Design matrix
# normal distribution
#' X <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
#' X = scale(X)/sqrt(n) #scale
#'
#' # Coefficient and response vectors
#' nonzero = sample(p, nspr)
#' beta = amplitude * (1:p %in% nonzero)
#' y = X %*% beta +  sigma*rnorm(n)
#'
#' # Missing values
#' X.obs <- X
#' p.miss <- 0.10
#' patterns <- runif(n*p)<p.miss #missing completely at random
#' X.obs[patterns] <- NA
#'
#' # ABSLOPE
#' lambda = create_lambda_bhq(ncol(X),fdr=0.10)
#' list.SLOB = ABSLOPE(X, y, lambda, a=2/p, b=1-2/p)
#' @export

ABSLOPE <- function(X, y, lambda, a, b, beta.start = NA, maxit = 300, case = 'MCAR', seed = NA, print_iter = FALSE, tol_em=1e-6, impute = 'PCA', sigma.known=NA, sigma.init=NA, scale=FALSE, method_na = 'lineq'){

  if(!is.na(seed)){set.seed(seed)}
  # missing pattern
  rindic = as.matrix(is.na(X))
  if(sum(rindic) > 0){ # missing data exist
    whichcolmissing = (1:ncol(rindic))[apply(rindic,2,sum)>0]
    missingcols = length(whichcolmissing)
  }
  if(sum(rindic) == 0){missingcols = 0} # no missingness

  p = ncol(X)

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


  # define functions to perform scaling
  scale_mean.w <- function(V, weight) {
    res <- sum(V * weight, na.rm = TRUE)/sum(weight[!is.na(V)])
  }
  scale_std.w <- function(V, weight) {
    res <- sqrt(sum(V^2 * weight, na.rm = TRUE)/sum(weight[!is.na(V)]))
  }
  # scaling
  if(scale){
    row.w = rep(1, nrow(X.sim))/nrow(X.sim)
    mean.w <- apply(X.sim, 2, scale_mean.w, row.w)
    X.sim <- t(t(X.sim) - mean.w)
    std.w <- apply(X.sim, 2, scale_std.w, row.w)*sqrt(n)
    X.sim<- t(t(X.sim)/std.w)
  }



  ##------------------------------------
  # initialization

  # beta
  if(is.na(beta.start)){ # initialization not given
    # LASSO cv
    objstart = cv.glmnet(X.sim,y)
    coeff = coef(objstart,s = "lambda.1se")
    beta = coeff[2:(p+1)]
  } else{beta = beta.start} # initialization given

  # sigma
  if(is.na(sigma.known)){ # real value of sigma is unknown
    if(is.na(sigma.init)){ # initialization not given
      #sigma = sd(y - X.sim %*% beta)
      sigma = sqrt(sum((y - X.sim %*% beta )^2)/(n-1))
    } else{sigma = sigma.init} # initialization given
  } else{sigma = sigma.known} # real value of sigma is known

  # theta
  theta = (sum(beta!=0)+a) / (a+b+p)


  # rank
  rk = p-rank(abs(beta),ties.method="max")+1

  # lambda
  lambda_sigma = lambda * sigma
  lambda_sigma_inv = lambda / sigma

  # c
  c = min( (sum(abs(beta)>0)+1) / sum (abs(beta[beta!= 0])) / lambda_sigma_inv[p], 1)
  if(!is.finite(c)){c = 1}

  # gamma
  gamma = (abs(beta)>0)*1
  pi.binom = 1/(1+(1-theta)/theta/c*exp(-abs(beta)*lambda_sigma_inv[rk]*(1-c)))
  #gamma = rbinom(1,1,pi.binom)

  # if NA -> mu & Sigma
  if(missingcols != 0){
    mu = apply(X.sim,2,mean)
    Sigma = linshrink_cov(X.sim)
  }



  ##------------------------------------

  # Iteration of SAEM
  # store the estimation of each iteration
  seqbeta = seqgamma = matrix(NA,nrow=ncol(X),ncol=(maxit+1))
  seqsigma = seqc = seqtheta = matrix(NA,nrow=1,ncol=(maxit+1))
  seqbeta[,1]=beta
  seqsigma[,1]=sigma
  seqgamma[,1]=gamma
  seqc[,1] = c
  seqtheta[,1] = theta

  cstop=1
  t=0
  while ((cstop>tol_em)*(t<maxit)|(t<20)){
    t=t+1

    if ((print_iter == TRUE) & (t%%100==0)){
      cat(sprintf('iteration = %i ', t))
      cat(sprintf('beta ='),beta)
      cat(sprintf(' sigma ='),sigma,'\n')
      cat(sprintf('Distance from last iter ='),cstop,'\n')
      #cat(sprintf(' c ='),c,'\n')
      #cat(sprintf(' pi.binom ='),pi.binom)
      #cat(sprintf(' gamma ='),gamma, '\n')
    }

    # step size
    if(t <= 20){eta = 1}else{eta = 1/(t - 20)}

    # Simulation step
    # 1) gamma
    W <- gamma * c+(rep(1,p) - gamma)
    pi.binom =  1/(1+(1-theta)/theta/c*exp(-abs(beta)*lambda_sigma_inv[rk]*(1-c)))
    gamma = rbinom(p,1,pi.binom)

    # 2) c
    W <- gamma * c+(rep(1,p) - gamma)
    a.gamma = 1 + sum(gamma == 1)
    b.gamma = sum(abs(beta)*lambda_sigma_inv[rk]*(gamma==1))
    # value0 = pgamma(0, shape=a.gamma, rate=b.gamma)
    # value1 = pgamma(1, shape=a.gamma, rate=b.gamma)
    # if(value0 != value1){c = rtrunc(1, "gamma", 0, 1, shape=a.gamma, rate=b.gamma)
    # }else{u = runif(1, min = 0, max = 1); c = u^(1/a.gamma)}
    if(a.gamma>1){
      if(b.gamma>0) c = rtrunc(1, "gamma", 0, 1, shape=a.gamma, rate=b.gamma)
      else c = rbeta(1,shape1 = a.gamma, shape2 = 1)
    }
    else c = runif(1,0,1)

    # 3) if NA -> Xmis
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
            nmcmc = 20
            xina <- xi[jna]
            for (m in (1:nmcmc)){
              xina.c <- mvrnorm(n = 1, mu=mi, Sigma =Oi)
              xi[jna] = xina.c
              alpha = dnorm(yi, mean=xi%*%beta , sd=sigma, log = FALSE)/dnorm(yi, mean=xi[jobs]%*%beta[jobs] + xina%*%beta[jna] , sd=sigma, log = FALSE)
              if (runif(1) < alpha){
                xina <- xina.c
              }
            }
            X.sim[i,jna] <- xina
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
            mu_tilde = solve(A, cc)

            B_inv = solve(A)

            Z <- mvrnorm(n = 1, mu=mu_tilde, Sigma=B_inv)
            X.sim[i,jna] <- Z / taui
          }
        }
      }

    }

    if(scale){
      X.sim = t(t(X.sim) * std.w)
      X.sim <- t(t(X.sim) + mean.w)
      mean.w <- apply(X.sim, 2, scale_mean.w, row.w)
      X.sim <- t(t(X.sim) - mean.w)
      std.w <- apply(X.sim, 2, scale_std.w, row.w)*sqrt(n)
      X.sim <- t(t(X.sim)/std.w)
    }

    # Stochastic Approximation step & Maximisation step
    beta.old = beta; sigma.old = sigma; theta.old =theta
    if(missingcols != 0){mu.old = mu; Sigma.old = Sigma}

    # MLE of completed data likelihood
    # 1) beta
    W <- gamma * c+(rep(1,p) - gamma)
    revW <- 1/W
    Xtemp <- sweep(X.sim,2,revW,'*'); # or: X.sim %*% diag(revW)
    z = slope_admm(Xtemp, y, rep(0, p), rep(0, p), lambda_seq = lambda_sigma, 1)$z
    beta<-revW*z
    cstop = sum((beta-beta.old)^2)

    # 2) sigma
    rk <- p-rank(abs(z),ties.method="max")+1
    RSS<- sum((y-X.sim %*%beta)^2)
    sum_lamwbeta <- sum(lambda[rk]*abs(z))
    if(is.na(sigma.known)){
      sigma <- (sqrt(sum_lamwbeta^2+4*n*RSS)+sum_lamwbeta)/(2*n)
    }

    #lambda
    lambda_sigma = lambda * sigma
    lambda_sigma_inv = lambda / sigma

    # 3) theta
    theta<-(sum(gamma==1)+a)/(a+b+p)

    # 4) mu et Sigma
    if(missingcols != 0){
      mu = apply(X.sim,2,mean)
      #Sigma = var(X.sim)*(n-1)/n
      Sigma = linshrink_cov(X.sim) # linear shrinkage
    }

    if(eta!=1){
      beta = beta.old + eta*(beta - beta.old)
      sigma = sigma.old + eta*(sigma - sigma.old)
      theta = theta.old + eta*(theta - theta.old)
      if(missingcols != 0){
        mu = mu.old + eta*(mu - mu.old)
        Sigma = Sigma.old + eta*(Sigma - Sigma.old)
      }
    }
    seqbeta[,t+1] = beta
    seqsigma[,t+1] = sigma
    seqgamma[,t+1] = gamma
    seqc[,t+1] = c
    seqtheta[,t+1] = theta
  }
  gamma.avg = rowMeans(seqgamma[,-(1:20)], na.rm = TRUE)
  #if(!beta.scale) {beta.new = beta * (gamma.avg>1/2)}
  beta.new = beta * (gamma.avg>1/2)
  return(list(X.sim=X.sim, beta=beta, seqbeta=seqbeta, beta.new =beta.new, gamma= gamma, gamma.avg = gamma.avg, seqgamma= seqgamma, pi.binom =pi.binom, theta = theta, seqtheta=seqtheta, sigma= sigma, seqsigma=seqsigma, c=c, seqc=seqc,mu=mu, Sigma = Sigma))
}
