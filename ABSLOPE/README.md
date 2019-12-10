# ABSLOPE

Codes for "Adaptive Bayesian SLOPE --  High-dimensional Model Selection with Missing Values" (2019, Jiang W., Bogdan M., Josse J., Miasojedow B., Rockova V., TraumaBase Group) 
 

In the main folder:
* main_prior_beta.R : Code to reproduce Figure 1 (page 10) --  the SLOPE prior and the ABLSOPE prior on a single coefficient.
* main_convergence.R : Code to reproduce Figure 3 (page 21) -- convergence plots for the coefficients with the proposed method on three simulated datasets.
* main_ABSLOPE_100.R : Code to reproduce Figure 4 (page 23) and Figure 5 (page 24) -- simulation study with n=p=100 and varied parametrisation. 
* main_SLOBC_500.R : Code to reproduce Figure 6 (page 25) and Figure 7 (page 26) -- simulation study with n=p=500 and varied parametrisation. 
* main_compare.R : Code to reproduce Figure 8 (page 29) -- method comparison.
* main_time.R : Code to reproduce Table 1 (page 30) -- comparison of average execution time.

Corresponding results of RData and the visualisation are included in the folders "RData" and "experiment_results" respectively.

The folder "ABSLOPE" includes the necessary functions of the algorithm:
* ABSLOPE.R : the main function containing the implementation of the proposed methodology, as introduced in Algorithm 1 (page 44).
* miss.SLOB.R : the function implementing the simplified algorithm, as introduced in Algorithm 2 (page 45).
* SLOBE_cpp_missing.cpp : C++ functions implementing the simplified algorithm, to accelerate the computing time.
 baseline_SLOBC.R : corporating C++ functions of the proposed methodology in R.
* NClasso.R : functions implementing the non-convex LASSO with missing values (Loh and Wainwright, 2012).
* other_methods_baseline.R : functions implementing other compared methods based on imputation.
* fct_others.R : other necessary functions including data generation, sequence of penalised parameters creation and optimisation solver.
* plot_tile.R and plot_line.R : functions for visualising the results.

