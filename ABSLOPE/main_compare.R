source('ABSLOPE/fct_others.R')
source('ABSLOPE/ABSLOPE.R')

source('ABSLOPE/miss.SLOB.R')
source('ABSLOPE/NClasso.R')
source('ABSLOPE/other_methods_baseline.R')
library(doParallel)

cl = makeCluster(40)
registerDoParallel(cl)

# -----------------------


n = 500
p = 500

simu.list = 1:40

p.miss.list = c(0.1,0.2)
corr.list = c(0,0.5,0.9)
signallevel.list = c(1,2,3,4)
nspr.list =c(10,20,30,40)


for(p.miss in p.miss.list){
  for(corr in corr.list){
    for(signallevel in signallevel.list){
      for(nspr in nspr.list){
        cat(sprintf('p.miss ='),p.miss,sprintf(' corr ='),corr,
            sprintf(' signallevel ='),signallevel,sprintf(' nspr ='),nspr,'\n')
        print('ABSLOPE')
        res.ABSLOPE <- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice'))  %dopar% unlist(baseline_ABSLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=500 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq',mec='MCAR',scale=TRUE))
        assign(paste0("ABSLOPE.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.ABSLOPE)
        save(list=ls(pattern=paste0("ABSLOPE.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr)), file = paste0("RData/mcar500/ABSLOPE.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr,'.Rdata'))
        
        print('miss.SLOB')
        res.miss.SLOB <- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice'))  %dopar% unlist(baseline_miss.SLOB(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq',mec='MCAR'))
        assign(paste0("miss.SLOB.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.miss.SLOB)
        save(list=ls(pattern=paste0("miss.SLOB.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr)), file = paste0("RData/mcar500/miss.SLOB.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr,'.Rdata'))
        
        print('SLOPE')
        res.SLOPE <- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice'))  %dopar% unlist(baseline_SLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,  sigma = 1, sigma.known=1,mec='MCAR'))
        assign(paste0("oth.SLOPE.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.SLOPE)
        
        print('LASSO.cv')
        res.LASSO.cv <- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice'))  %dopar% unlist(baseline_LASSO.cv(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,mec='MCAR'))
        assign(paste0("oth.LASSO.cv.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.LASSO.cv)
        
        print('ncLASSO')
        res.ncLASSO<- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice'))  %dopar% unlist(baseline_ncLASSO(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,mec='MCAR',oracle=TRUE))
        assign(paste0("ncLASSO.mcar.oracle.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.ncLASSO)
        save(list=ls(pattern=paste0("ncLASSO.mcar.oracle.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr)), file = paste0("RData/mcar500/ncLASSO.mcar.oracle.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr,'.Rdata'))
        
        print('adaLASSO')
        res.adaLASSO<- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice'))  %dopar% unlist(baseline_adaLASSO(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,mec='MCAR'))
        assign(paste0("oth.adaLASSO.cv.mcar.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.adaLASSO)
        save(list=ls(pattern="^oth"), file = paste0("RData/mcar500/oth.res.mcar.Rdata"))
        
      }
    }
  }
}


stopCluster(cl)


#--------
file_SLOB=as.list(dir(path = 'RData/mcar500/', pattern="miss.SLOB.*"))
file_SLOB=lapply(file_SLOB, function(x) paste0('RData/mcar500/', x))
lapply(file_SLOB,load,.GlobalEnv)
file_ABSLOPE=as.list(dir(path = 'RData/mcar500/', pattern="ABSLOPE.*"))
file_ABSLOPE=lapply(file_ABSLOPE, function(x) paste0('RData/mcar500/', x))
lapply(file_ABSLOPE,load,.GlobalEnv)
file_ncLASSO=as.list(dir(path = 'RData/mcar500/', pattern="ncLASSO.*"))
file_ncLASSO=lapply(file_ncLASSO, function(x) paste0('RData/mcar500/', x))
lapply(file_ncLASSO,load,.GlobalEnv)
file_oth.res=as.list(dir(path = 'RData/mcar500/', pattern="oth.*"))
file_oth.res=lapply(file_oth.res, function(x) paste0('RData/mcar500/', x))
lapply(file_oth.res,load,.GlobalEnv)

p.miss.list = c(0.1)
corr.list = c(0,0.5,0.9)
signallevel.list = c(1,2,3,4)
nspr.list =c(10,20,30,40)

crit.list <- c('pr', 'fdr', 'bias_beta', 'bias_sigma', 'MSE')
simu.list = 1:40
algo.names.read <- c('ABSLOPE.mcar','miss.SLOB.mcar','ncLASSO.oracle.mcar','oth.SLOPE.mcar','oth.LASSO.cv.mcar','oth.adaLASSO.cv.mcar')
algo.names <- c('ABSLOPE','SLOBE','ncLASSO','MeanImp+SLOPE','MeanImp+LASSO','MeanImp+adaLASSO')

# aggregate them in a single array
# results (7 dim): 1 method ; 2 p.miss ; 3 corr ; 4 signallevel ; 5 nspr ; 6 crit ; 7 nsimu
results <- array(NA, dim = c(length(algo.names), length(p.miss.list), 
                             length(corr.list), length(signallevel.list),
                             length(nspr.list), 5, length(simu.list)),
                 dimnames = list(method = algo.names, p.miss = p.miss.list,
                                 corr = corr.list, signallevel = signallevel.list,
                                 nspr = nspr.list, crit = crit.list ,nsimu = simu.list))
for(c1 in 1:length(algo.names)){
  for(c2 in 1:length(p.miss.list)){
    for(c3 in 1:length(corr.list)){
      for(c4 in 1:length(signallevel.list)){
        for(c5 in 1:length(nspr.list)){
          method = algo.names.read[c1]
          p.miss = p.miss.list[c2]
          corr = corr.list[c3]
          signallevel = signallevel.list[c4]
          nspr = nspr.list[c5]
          res = get(paste0(method,".", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr))
          results[c1,c2,c3,c4,c5, , ] <- res[,1:40]
        }
      }
    }
  }
}

save(results, file="RData/mcar500/results2.Rda")
load(file="RData/mcar500/results2.Rda")

Mresults <- melt(results, value.name="value")
# head(Mresults)

plot_addr= "experiment_results/mcar500bis2/"
ifelse(!dir.exists(file.path(plot_addr)),
       dir.create(file.path(plot_addr)),
       FALSE)
# 10% missingness, no correlation
pm = 0.1
cr = 0 
plot_tile(Mresults, pm, cr, plot_addr)

# 10% missingness, 0.5 correlation
pm = 0.1
cr = 0.5 
plot_tile(Mresults, pm, cr, plot_addr)

# 10% missingness, 0.9 correlation
pm = 0.1
cr = 0.9
plot_tile(Mresults, pm, cr, plot_addr)

