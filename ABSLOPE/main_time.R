source('ABSLOPE/fct_others.R')
source('ABSLOPE/ABSLOPE.R')
source('ABSLOPE/miss.SLOB.R')
source('ABSLOPE/NClasso.R')
source('ABSLOPE/other_methods_baseline.R')
source('ABSLOPE/baseline_SLOBC.R')



library(microbenchmark)

n = 100
p = 100
p.miss = 0.1
corr = 0
signallevel =3
nspr=5
sigma = 1
simu = 29

res.time = microbenchmark(#baseline_ABSLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE),
               #baseline_miss.SLOB(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE),
               baseline_SLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,  sigma = 1, sigma.known=1),
               baseline_ncLASSO(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu),
               baseline_LASSO.cv(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu),
               baseline_adaLASSO(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu),
               times = 200)
save(res.time, file='res_time.RData')

res.time.ABSLOPE =  microbenchmark(baseline_ABSLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq'),
                                   baseline_miss.SLOB(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq'),
                                   times = 200)
save(res.time.ABSLOPE, file='RData/time_compare/res_time_ABSLOPE.RData')

res.time.SLOBC =  microbenchmark(baseline_SLOBC(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu ,sigma = 1),
                                   times = 200)
save(res.time.SLOBC, file='RData/time_compare/res_time_SLOBC.RData')


#-------------------------------------------------------
n = 500
p = 500
p.miss = 0.1
corr = 0
signallevel =3
nspr=10
sigma = 1
simu = 29
res.time500 = microbenchmark(baseline_ABSLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq'),
                             baseline_miss.SLOB(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu, a = 2/p, b = 1-2/p,tol_em=1e-3,maxit=300 ,sigma = 1 ,print_iter=FALSE,method_na = 'lineq'),
                             baseline_SLOBC(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu ,sigma = 1),
                             baseline_SLOPE(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu,  sigma = 1, sigma.known=1),
  baseline_ncLASSO(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu),
  baseline_LASSO.cv(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu),
  baseline_adaLASSO(n=n, p=p, nspr=nspr, p.miss=p.miss, mu=rep(0,p), Sigma=toeplitz(corr^(0:(p-1))), signallevel=signallevel, nb.seed=simu),
  times = 200)
save(res.time500, file='RData/time_compare/res_time_500.RData')


