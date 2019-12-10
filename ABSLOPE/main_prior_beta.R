# plot prior SLOPE & ABLSOPE
source('ABSLOPE/ABSLOPE.R')
source('ABSLOPE/fct_others.R')
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())
library(gbm)

fun_SLOPE <- function(b,beta,j,sig,X,y) {
  #beta[j]=b
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
  rk <- p-rank(abs(beta),ties.method="max")+1
  lambda_ord = lambda[rk]
  exp(-lambda_ord[j]/sig*abs(b))
}

fun_LASSO <- function(b,beta,j,sig,X,y) {
  #beta[j]=b
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
  lambda_ord = max(lambda)/2
  exp(-lambda_ord/sig*abs(b))
}

fun_ABSLOPE <- function(b,beta,j,sig,X,y) {
  #beta[j]=b
  gamma = (abs(beta)>0)*1
  lambda = create_lambda_bhq(ncol(X),fdr=0.10)
  c =  (sum(abs(beta)>0)+1) / sum (abs(beta[beta!= 0]))
  W <- gamma * c+(rep(1,p) - gamma)
  z = W * beta
  rk <- p-rank(abs(z),ties.method="max")+1
  lambda_ord = lambda[rk]
  exp(-lambda_ord[j]/sig*W[j]*abs(b))
}

set.seed(100)
n = 100
p = 100
signallevel = 3
p.miss = 0
mu = rep(0,p)
Sigma = diag(p)
nspr=10
data.list = data.generation(n, p, nspr, p.miss,mu,Sigma, signallevel)
X = data.list$X.obs
y = data.list$y
beta = data.list$beta
sigma = 1
j=12
b  <- seq(-10,10,by=0.1)
f1 <- fun_SLOPE(b=b,beta=beta,j=j,sig=sigma,X=X,y=y)
f1 <- f1/sum(f1)
f2 <- fun_ABSLOPE(b=b,beta=beta,j=j,sig=sigma,X=X,y=y)
f2 <- f2/sum(f2)
f3 <- fun_LASSO(b=b,beta=beta,j=j,sig=sigma,X=X,y=y)
f3 <- f3/sum(f3)

df <- data.frame(b,f1,f2)
png("experiment_results/prior_compare.png", width = 6, height = 5, units = 'in',res = 300)
ggplot(df, aes(b)) +                    
  geom_line(aes(y=f1, color="SLOPE")) +  
  geom_line(aes(y=f2, color="ABSLOPE")) +
  scale_color_manual(name = "Prior", 
                     values = c("SLOPE" = "blue", "ABSLOPE" = "red"))+
  xlab('beta') + ylab('prior density')
dev.off()


j=11
b  <- seq(-10,10,by=0.1)
f1 <- fun_SLOPE(b=b,beta=beta,j=j,sig=sigma,X=X,y=y)
f1 <- f1/sum(f1)
f2 <- fun_ABSLOPE(b=b,beta=beta,j=j,sig=sigma,X=X,y=y)
f2 <- f2/sum(f2)
f3 <- fun_LASSO(b=b,beta=beta,j=j,sig=sigma,X=X,y=y)
f3 <- f3/sum(f3)

df <- data.frame(b,f1,f2)
png("experiment_results/prior_compare0.png", width = 6, height = 5, units = 'in',res = 300)
ggplot(df, aes(b)) +                    
  geom_line(aes(y=f1, color="ABSLOPE & SLOPE")) +  
  #geom_line(aes(y=f2, color="ABSLOPE")) +
  scale_color_manual(name = "Prior", 
                     values = c("ABSLOPE & SLOPE" = "red"))+
  xlab('beta') + ylab('prior density')
dev.off()