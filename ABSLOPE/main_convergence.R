source('ABSLOPE/fct_others.R')
source('ABSLOPE/ABSLOPE.R')
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())
library(gbm)

# Several curves of convergence 
#----------with NA-------------------------

# Convergence of ABSLOPE

# beta
set.seed(100)
n=100
p=100
nspr=0.1*p 
mu = rep(0,p)
Sigma = diag(p)
signallevel=3
p.miss = 0.1
amplitude = signallevel*sqrt(2*log(p))  # signal amplitude (for noise level = 1)
nonzero = sample(p, nspr)
beta = amplitude * (1:p %in% nonzero)
lambda = create_lambda_bhq(ncol(X),fdr=0.05)
maxit = 500

# Design matrix 1
set.seed(100)
# normal distribution
X1 <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
X1 = scale(X1)/sqrt(n)
# response vectors
y1 = X1 %*% beta +  sigma*rnorm(n)
# Missing values
X.obs <- X1
patterns <- runif(n*p)< p.miss # missing completely at random
X.obs[patterns] <- NA
list.ABSLOPE1 = ABSLOPE(X.obs,y1,lambda, a = 2/p, b = 1-2/p, maxit = maxit,print_iter = FALSE,method_na = 'lineq')

# Design matrix 2
set.seed(250)
# normal distribution
X2 <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
X2 = scale(X2)/sqrt(n)
# response vectors
y2 = X2 %*% beta +  sigma*rnorm(n)
# Missing values
X.obs <- X2
patterns <- runif(n*p)< p.miss # missing completely at random
X.obs[patterns] <- NA
list.ABSLOPE2 = ABSLOPE(X.obs,y2,lambda, a = 2/p, b = 1-2/p, maxit = maxit,print_iter = FALSE,method_na = 'lineq')

# Design matrix 3
set.seed(400)
# normal distribution
X3 <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma) + matrix(rep(mu,n), nrow=n, byrow = TRUE)
X3 = scale(X3)/sqrt(n)
# response vectors
y3 = X3 %*% beta +  sigma*rnorm(n)
# Missing values
X.obs <- X3
patterns <- runif(n*p)< p.miss # missing completely at random
X.obs[patterns] <- NA
list.ABSLOPE3 = ABSLOPE(X.obs,y3,lambda, a = 2/p, b = 1-2/p, maxit = maxit,print_iter = FALSE,method_na = 'lineq')


df1 = t(list.ABSLOPE1$seqbeta) %>% as.data.frame() %>% dplyr::select(V45,V46,V47) %>% mutate(iteration=0:maxit)  %>%  melt(id.vars = 'iteration', variable.name = 'beta') 
df1$beta <-factor(df1$beta,labels=c("V45"=paste0('\u03b2',45),
                                    "V46"=paste0('\u03b2',46),
                                    "V47"=paste0('\u03b2',47)))
df1$Label = as.factor("Est. of 1st simulation")
df2 = t(list.ABSLOPE2$seqbeta) %>% as.data.frame() %>% dplyr::select(V45,V46,V47) %>% mutate(iteration=0:maxit)  %>%  melt(id.vars = 'iteration', variable.name = 'beta') 
df2$beta <-factor(df2$beta,labels=c("V45"=paste0('\u03b2',45),
                                    "V46"=paste0('\u03b2',46),
                                    "V47"=paste0('\u03b2',47)))
df2$Label = as.factor("Est. of 2nd simulation")
df3 = t(list.ABSLOPE3$seqbeta) %>% as.data.frame() %>% dplyr::select(V45,V46,V47) %>% mutate(iteration=0:maxit)  %>%  melt(id.vars = 'iteration', variable.name = 'beta') 
df3$beta <-factor(df3$beta,labels=c("V45"=paste0('\u03b2',45),
                                    "V46"=paste0('\u03b2',46),
                                    "V47"=paste0('\u03b2',47)))
df3$Label = as.factor("Est. of 3rd simulation")
truebeta <- data.frame(iteration=rep(0:maxit,3), beta = c(rep(paste0('\u03b2',45),maxit+1), rep(paste0('\u03b2',46),maxit+1), rep(paste0('\u03b2',47),maxit+1)), value = c(rep(beta[45],maxit+1), rep(beta[46],maxit+1),rep(beta[47],maxit+1)))
truebeta$Label = as.factor("True value")
df <- rbind(truebeta,df1, df2, df3)
png("experiment_results/convergence_beta.png", width = 8, height = 4, units = 'in',res = 300)
df %>%  ggplot() + aes(iteration,value,linetype=Label,color=Label) + geom_line() + facet_grid(. ~ beta) +  
  ylab(expression(paste("Estimate of ",beta))) + 
  #theme(strip.text = element_text(size=12), axis.title=element_text(size=14))+
  scale_linetype_manual(values=c(2,1,1,1)) +
  scale_color_manual(values=c(1,2,3,4))
dev.off()


