source('ABSLOPE/fct_others.R')
source('ABSLOPE/baseline_SLOBC.R')
library(doParallel)
# cl = makeCluster(12)
# #clusterCall(cl, function() { library(Rcpp);Rcpp::sourceCpp('ABSLOPE/SLOBE_cpp_missing.cpp') })
# registerDoParallel(cl)

# -----------------------


n = 500
p = 500

simu.list = 1:200

p.miss.list = c(0.1,0.2,0.3)
corr.list = c(0,0.5,0.9)
signallevel.list = c(1,2,3,4)
nspr.list =seq(10,60,10)


for(p.miss in p.miss.list){
  for(corr in corr.list){
    for(signallevel in signallevel.list){
      for(nspr in nspr.list){
        destfile = paste0("RData/SLOBC500/SLOBC.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr,'.Rdata')
        if(!file.exists(destfile)){
          cat(sprintf('p.miss ='),p.miss,sprintf(' corr ='),corr,
              sprintf(' signallevel ='),signallevel,sprintf(' nspr ='),nspr,'\n')
          print('ABSLOPE')
          res.ABSLOPE <- foreach(simu = simu.list, .combine=cbind, .packages=c('glmnet','truncdist','nlshrink','MASS','missMDA','SLOPE','mice','Rcpp')
          )  %do% unlist(baseline_SLOBC(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,sigma = 1,mec='MCAR'))
          assign(paste0("SLOBC.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr),res.ABSLOPE)
          save(list=ls(pattern=paste0("SLOBC.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr)), file = paste0("RData/SLOBC500/SLOBC.", p.miss*100,'.',corr*10,'.',signallevel,'.',nspr,'.Rdata'))
        }
      }
    }
  }
}

#stopCluster(cl)

source('ABSLOPE/plot_tile.R')
source('ABSLOPE/plot_line.R')

file_SLOBC=as.list(dir(path = 'RData/SLOBC500/', pattern="SLOBC.*"))
file_SLOBC=lapply(file_SLOBC, function(x) paste0('RData/SLOBC500/', x))
lapply(file_SLOBC,load,.GlobalEnv)

p.miss.list = c(0.1,0.2,0.3)
corr.list = c(0.5)
signallevel.list = c(3)
nspr.list =seq(10,60,10)
crit.list <- c('pr', 'fdr', 'bias_beta', 'bias_sigma', 'MSE')
simu.list = 1:200
algo.names.read <- c('SLOBC')
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

Mresults <- melt(results, value.name="value")
plot_addr= "experiment_results/SLOBC500/"
# vary %NA
# strong signal strength, 0.5 correlation
sl = 3
cr = 0.5
plot_line_percentna(Mresults, sl, cr, plot_addr)

#----

p.miss.list = c(0.1)
corr.list = c(0,0.5,0.9)
signallevel.list = c(3)
nspr.list =seq(10,60,10)
crit.list <- c('pr', 'fdr', 'bias_beta', 'bias_sigma', 'MSE')
simu.list = 1:200
algo.names.read <- c('SLOBC')
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

Mresults <- melt(results, value.name="value")
plot_addr= "experiment_results/SLOBC500/"
#vary corr
pm = 0.1
sl = 3
plot_line_corr(Mresults, pm, sl, plot_addr)

#----

p.miss.list = c(0.1)
corr.list = c(0)
signallevel.list = c(1,2,3,4)

nspr.list =seq(10,60,10)
crit.list <- c('pr', 'fdr', 'bias_beta', 'bias_sigma', 'MSE')
simu.list = 1:200
algo.names.read <- c('SLOBC')
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

Mresults <- melt(results, value.name="value")
plot_addr= "experiment_results/SLOBC500/"

# varying signal strength
# 10% missingness, no correlation
pm = 0.1
cr = 0
plot_line_signallevel(Mresults, pm, cr, plot_addr)