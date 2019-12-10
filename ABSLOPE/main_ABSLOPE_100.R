source('ABSLOPE/fct_others.R')
source('ABSLOPE/ABSLOPE.R')
source('ABSLOPE/other_methods_baseline.R')
library(doParallel)

cl = makeCluster(12)
registerDoParallel(cl)

# -----------------------


n = 100
p = 100

simu.list = 1:200

p.miss.list = c(0.1,0.2,0.3)
corr.list = c(0,0.5,0.9)
signallevel.list = c(1,2,3,4)
nspr.list =c(5,10,15,20)


for(p.miss in p.miss.list){
  for(corr in corr.list){
    for(signallevel in signallevel.list){
      for(nspr in nspr.list){
        cat(sprintf('p.miss ='),p.miss,sprintf(' corr ='),corr,
            sprintf(' signallevel ='),signallevel,sprintf(' nspr ='),nspr,'\n')
        print('ABSLOPE')
        res.ABSLOPE <- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE'))  %dopar% unlist(baseline_ABSLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq'))
        assign(paste0("ABSLOPE.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.ABSLOPE)
        save(list=ls(pattern=paste0("ABSLOPE.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr)), file = paste0("RData/MCAR100/ABSLOPE.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr,'.Rdata'))
      }
    }
  }
}

stopCluster(cl)


file_ABSLOPE=as.list(dir(path = 'RData/ABSLOPE100/', pattern="ABSLOPE.*"))
file_ABSLOPE=lapply(file_ABSLOPE, function(x) paste0('RData/ABSLOPE100/', x))
lapply(file_ABSLOPE,load,.GlobalEnv)

p.miss.list = c(0.1,0.2,0.3)
corr.list = c(0,0.5,0.9)
signallevel.list = c(1,2,3,4)
nspr.list =c(5,10,15,20)
crit.list <- c('pr', 'fdr', 'bias_beta', 'bias_sigma', 'MSE')
simu.list = 1:200
algo.names.read <- c('ABSLOPE')
algo.names <- c('ABSLOPE')

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
          results[c1,c2,c3,c4,c5, , ] <- res
        }
      }
    }
  }
}

save(results, file="RData/ABSLOPE100/results.Rda")

source('ABSLOPE/plot_tile.R')
Mresults <- melt(results, value.name="value")
plot_addr= "experiment_results/ABSLOPE100/"

source('ABSLOPE/plot_line.R')
#----
# varying signal strength
# 10% missingness, no correlation
pm = 0.1
cr = 0
plot_line_signallevel(Mresults, pm, cr, plot_addr)
# 10% missingness, 0.5 correlation
pm = 0.1
cr = 0.5 
plot_line_signallevel(Mresults, pm, cr, plot_addr)
# 10% missingness, 0.9 correlation
pm = 0.1
cr = 0.9
plot_line_signallevel(Mresults, pm, cr, plot_addr)


#---
# vary p.miss
# strong signal strength, no correlation
sl = 3
cr = 0
plot_line_percentna(Mresults, sl, cr, plot_addr)

# strong signal strength, 0.5 correlation
sl = 3
cr = 0.5
plot_line_percentna(Mresults, sl, cr, plot_addr)

# strong signal strength, 0. 9 correlation
sl = 3
cr = 0.9
plot_line_percentna(Mresults, sl, cr, plot_addr)


#------
#vary corr
pm = 0.1
sl = 3
plot_line_corr(Mresults, pm, sl, plot_addr)