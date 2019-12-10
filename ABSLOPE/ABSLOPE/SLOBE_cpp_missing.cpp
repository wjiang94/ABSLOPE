#include<RcppArmadillo.h>
#include<math.h>
#include<stdlib.h>
#include<numeric>
#include<algorithm>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

int evaluateProx(double *y, double *lambda, double *x, size_t n, int *order);

/* ----------------------------------------------------------------------- */
int evaluateProx(double *y, double *lambda, double *x, size_t n, int *order)
/* ----------------------------------------------------------------------- */
{  double   d;
   double  *s     = NULL;
   double  *w     = NULL;
   size_t  *idx_i = NULL;
   size_t *idx_j = NULL;
   size_t  i,j,k;
   int      result = 0;

   /* Allocate memory */
   s     = (double *)malloc(sizeof(double) * n);
   w     = (double *)malloc(sizeof(double) * n);
   idx_i = (size_t *)malloc(sizeof(size_t) * n);
   idx_j = (size_t *)malloc(sizeof(size_t) * n);

   if ((s != NULL) && (w != NULL) && (idx_i != NULL) && (idx_j != NULL))
   {
      k = 0;
      for (i = 0; i < n; i++)
      {
         idx_i[k] = i;
         idx_j[k] = i;
         s[k]     = y[i] - lambda[i];
         w[k]     = s[k];

         while ((k > 0) && (w[k-1] <= w[k]))
         {  k --;
            idx_j[k] = i;
            s[k]    += s[k+1];
            w[k]     = s[k] / (i - idx_i[k] + 1);
         }

         k++;
      }

      if (order == NULL)
      {  for (j = 0; j < k; j++)
         {  d = w[j]; if (d < 0) d = 0;
            for (i = idx_i[j]; i <= idx_j[j]; i++)
            {  x[i] = d;
            }
         }
      }
      else
      {  for (j = 0; j < k; j++)
         {  d = w[j]; if (d < 0) d = 0;
            for (i = idx_i[j]; i <= idx_j[j]; i++)
            {  x[order[i]] = d;
            }
         }
      }
   }
   else
   {  result = -1;
   }

   /* Deallocate memory */
   if (s     != NULL) free(s);
   if (w     != NULL) free(w);
   if (idx_i != NULL) free(idx_i);
   if (idx_j != NULL) free(idx_j);

   return result;
}

// Comparator class for argsort function
class CompareByNumericVectorValues {
	private:
		const NumericVector* _values;
  
	public:
		CompareByNumericVectorValues(const NumericVector* values) { _values = values;}
		bool operator() (const int& a, const int& b) {return ((*_values)[a] > (*_values)[b]);}
};


// Writes down in IntegerVector ord sequnce of indexes,
// such that w[ord[i]] >= w[ord[j]] whenever i <= j 
// for the given NumericVector w.
void argsort(const NumericVector& w, IntegerVector ord) {
	std::iota(ord.begin(), ord.end(), 0);
	CompareByNumericVectorValues comp = CompareByNumericVectorValues(&w);
	std::sort(ord.begin(), ord.end(), comp);
}

// Computes proximal step of SLOPE(lambda) norm from point y
// NumericVector y is not assumed to be sorted, sorting is performed within function
NumericVector prox_sorted_L1_C(NumericVector y, NumericVector lambda) {
 	size_t n = y.size();
 	NumericVector x(n);
 	IntegerVector order(n);
 	argsort(abs(y),order);
 	IntegerVector sign_y= sign(y);
 	y = abs(y);
 	y.sort(true);
 	evaluateProx(y.begin(), lambda.begin(), x.begin(), n, NULL);
 	NumericVector res(n);
	for(int k=0;k<n;k++){
		res[order[k]]= sign_y[order[k]]*x[k];
	}
	return res;
}


void create_lambda(NumericVector& lam, int p, double FDR) {
	NumericVector h(p);
	for (double i = 0.0; i < h.size(); ++i) {
		h[i] = 1 - (FDR* (i+1)/(2*p));
	}
	lam = qnorm(h);
}

double EX_trunc_gamma(double a ,double b ){
  double c = exp(Rf_pgamma(1.0, a+1, 1.0/b, 1, 1) - Rf_pgamma(1.0, a, 1.0/b, 1, 1));
  c /= b;
  c *= a;
  return c;
}


arma::vec slope_admm(const mat& X, const vec& Y, NumericVector& lambda,
				const int& p, const double& rho, int max_iter=500, double tol_inf = 1e-08) {
	
		// Precompute M = (X^TX + rho I)^{-1} 
		// and MXtY = M * X^T * Y for proximal steps of quadratic part
		mat M = X.t() * X;
		for (int i=0; i<p; ++i) {
			M.at(i,i) += rho;
		}
		M = M.i();
		vec MXtY = M * (X.t() * Y);
		NumericVector lam_seq_rho = lambda/rho;
		
		// Prepare variables before starting ADMM loop
		int i=0;
		vec x = zeros(p);
		vec z = zeros(p);
		vec u = zeros(p);
		NumericVector z_new = NumericVector(p);
		vec z_new_arma = zeros(p);
		NumericVector x_plus_u(p);
		double dual_feas, primal_feas;
		
		// ADMM loop
		while (i < max_iter) {
			x = MXtY + M*(rho*(z - u));
			x_plus_u = as<NumericVector>(wrap(x+u));
			z_new = prox_sorted_L1_C(x_plus_u, lam_seq_rho);
			z_new_arma = as<arma::vec>(z_new);
			u += (x - z_new_arma);
			
			dual_feas = arma::norm(rho*(z_new_arma - z));
			primal_feas = arma::norm(z_new_arma - x);
			
			z = z_new_arma;
			if (primal_feas < tol_inf && dual_feas < tol_inf){
				i = max_iter;
			}

			++i;
		}
		
		return z;
}

void div_X_by_w(mat& X_div_w, const mat& X, const vec& w_vec, const int& n, const int& p) {
	for (int i=0; i<n; ++i) {
		for (int j=0; j<p; ++j) {
			X_div_w.at(i,j) = X.at(i,j) / w_vec(j);
		}
	}
}


// [[Rcpp::export]]
void impute_mean(mat& X, const int& n, const int &p) {
	double colmean;
	int non_na_rows;
	for (int c_num=0; c_num<p; c_num++) {
		colmean = 0.0;
		non_na_rows = 0;
		for (int r_num=0; r_num<n; r_num++) {
			if (arma::is_finite(X.at(r_num, c_num))) {
				colmean += X.at(r_num, c_num);
				non_na_rows += 1;
			}
		}
		colmean /= non_na_rows;
		for (int r_num=0; r_num<n; r_num++) {
			if (!arma::is_finite(X.at(r_num, c_num))) {
				X.at(r_num, c_num) = colmean;
			}
		}
	}
}

void linshrink_cov(mat &X, mat &S, const int& n, const int& p) {
	
	rowvec means = sum(X)/n;
	S = X.t() * X;
	for (int i=0; i<p; ++i) {
		for (int j=0; j<p; ++j) {
			S.at(i,j) -= n*means[i]*means[j];
			S.at(i,j) /= (n-1);
		}
	}
	double m = trace(S)/p;
	double d2 = norm(S, "fro");
	double b_bar2 = pow(d2,2);
	d2 = (b_bar2 - p*pow(m,2))/p;
	b_bar2 *= n;
	double prod;
	double sum_prod;
	for (int i=0; i<p; ++i) {
		for (int j=0; j<p; ++j) {
			sum_prod = 0.0;
			for (int k=0; k<n; ++k){
				prod = (X.at(k,i) - means[i]) * (X.at(k,j) - means[j]);
				b_bar2 += prod * prod;
				sum_prod += prod;
			}
			b_bar2 -= 2 * S.at(i,j) * sum_prod;
		}
	}
	b_bar2 /= (p * pow(n-1,2));
	double b2 = b_bar2 < d2 ? b_bar2 : d2;
	double a2 = d2 - b2;
	S *= (a2/d2);
	m *= (b2/d2);
	for (int i=0; i<p; ++i) {
		S.at(i,i) += m;
	}
}

void impute_row_advance(const vec& beta, mat& X, const vec& Y, const mat& S, const double& sigma_sq, 
				const LogicalMatrix& XisFin, const int& n, const int& p, const int& row, 
				const std::vector<int> nanCols, const vec& m, const vec& tau_sq) {
	
	int l = nanCols.size();	
	int s,t;
	mat A = zeros(l,l);
	vec u = zeros(l);
	double r = Y[row];
	int u_ind = 0;
	for (int i=0; i<p; ++i) {
		if (!XisFin.at(row, i)) {
			for (int j=0; j<p; ++j) {
				if(XisFin.at(row,j)) {
					//Rcout << "(row, j): (" << row << "," << j << ")\n";
					u[u_ind] += X.at(row,j) * S.at(j,i); 
				}
			}
			u_ind += 1;
		}
		else {
			r -= beta[i] * X.at(row,i);
		}
	}
	//Rcout << "r: " << r << "\n";
	//Rcout << "u: " << u << "\n";
	for (int i=0; i<l; ++i) {
		for (int j=0; j<l; ++j) {
			if (i == j) {
				A.at(i,j) = 1.0;
			}
			else {
				s = nanCols[i];
				t = nanCols[j];
				A.at(i,j) = (beta[s]*beta[t]/sigma_sq + S.at(s,t))/tau_sq[s];
			}
		}
	}

	vec b = zeros(l);
	for (int i=0; i<l; ++i) {
		t = nanCols[i];
		b[i] =  ((r * beta[t])/sigma_sq + m[t] - u[i])/tau_sq[t];
	}	
		
	vec sol = solve(A,b);
	//Rcout << "A: " << A << "\n";
	//Rcout << "b: " << b << "\n";
	//Rcout << "sol: " << sol << "\n";
	for (int i=0; i<l; ++i) {
		t = nanCols[i];
		//Rcout << "imputed column: " << t << "\n";
		X.at(row, t) = sol[i];
	}
}


// Imputation procedure for SLOBE with missing values
void impute_advance(const vec &beta, mat& X, const vec &Y, const mat& S, const double& sigma_sq,
				const int& n, const int& p, const vec& mu, const LogicalMatrix& XisFin, 
				const std::vector<int> anyNanXrows, const std::vector<std::vector<int> > nanIndInRows) {

	vec tau_sq = zeros(p);
	tau_sq = square(beta)/sigma_sq;
	for (int i=0; i<p; ++i) {
		tau_sq[i] = tau_sq[i] + S.at(i,i);
	}
	
	//Rcout << "tau_sq: " << tau_sq << "\n";
	vec m = zeros(p);
	for (int i=0; i<p; ++i) {
		for (int j=0; j<p; ++j) {
			m[i] += mu[j] * S.at(i,j);
		}
	}
	//Rcout << "m: " << m << "\n";

	for (int i=0; i<anyNanXrows.size(); ++i) {
		//Rcout << "rownum: " << i << "\n";
		impute_row_advance(beta, X, Y, S, sigma_sq, XisFin, n, p, anyNanXrows[i], nanIndInRows[i], m, tau_sq);
	}			
		
}




// [[Rcpp::export]]
List SLOBE_ADMM_approx_missing(NumericVector start, mat Xmis, vec Y, double a_prior, double b_prior, double sigma = 1.0, 
				double FDR = 0.05, double tol = 1e-04, bool known_sigma = false, int max_iter=100, bool verbose = true) {

	// Initialize variables
	int p = start.length();
	int n = Y.size();
	NumericVector beta = clone(start);
	NumericVector beta_new(p);
	vec beta_arma = as<vec>(beta);
	NumericVector w(p, 1.0);
	vec w_vec = ones(p);
	NumericVector wbeta(p);
	NumericVector gamma(p);
	NumericVector gamma_h(p);
	NumericVector b_sum_h(p);
	NumericVector lambda_sigma(p);
	IntegerVector order(p);
	mat X_div_w = zeros(n,p);
  	double error = 0.0;
	double swlambda = 0.0;
	double RSS = 0.0;
	
        // First imputation
        bool anyNanInRow;
  
        LogicalMatrix XisFin(n,p);
        std::vector<int> anyNanXrows;
        std::vector<std::vector<int> > nanIndicesInRow; 
        for (int i=0; i<n; ++i) {
            anyNanInRow = false;
            std::vector<int> nanInd;
            for (int j=0; j<p; ++j) { 
                XisFin.at(i,j) = arma::is_finite(Xmis.at(i,j));
                if (!XisFin.at(i,j)) {
                nanInd.push_back(j);
                anyNanInRow = true;
               }
	      }
          if (anyNanInRow) {
             anyNanXrows.push_back(i);
             nanIndicesInRow.push_back(nanInd);
          }
       }
  
       mat X = mat(Xmis);
       impute_mean(X, n, p);
       X = normalise(X);
       mat Sigma = zeros(p,p);
       linshrink_cov(X, Sigma, n, p);
  
       mat S = inv_sympd(Sigma);
       vec mu = zeros(p);
       for (int i=0; i<p; ++i) {
           for (int j=0; j<n; ++j) {
               mu.at(i) += X.at(j,i); 
            }
        mu[i] /= n;
       }
  
  

	// Compute vector lambda based on BH procedure
	NumericVector lambda(p);
	create_lambda(lambda, p, FDR);

	// Initialize c, theta
	double sstart = sum(start != 0);
	double c = 0.0;
	if (sstart > 0) {
		double h = (sstart+1)/(abs(sstart * lambda[p-1] * sigma));
		c = (h < 0.9) ? h : 0.9;
	}
	else
		c = 0.9;	
  
	double theta = (sstart + a_prior)/(a_prior + b_prior + p);

	// Start main loop
	bool converged = false;
	int iter = 0;
	while (iter < max_iter) {
	  if(verbose){
		  Rcout << "Iteracja: " << iter <<"\n" ;
	  }
		wbeta = w * abs(beta);
		argsort(wbeta, order);

		// For the version with unknown sigma, compute it first
		if (!known_sigma) {
			RSS = sum(pow((X * beta_arma - Y),2));
			swlambda = sum(wbeta.sort(true) * lambda);
			sigma = (swlambda + sqrt(pow(swlambda, 2.0) + 4*n*RSS))/(2*n);
		}

		// compute new gamma
		gamma_h = abs(beta[order]) * (c-1) * lambda / sigma ;
		gamma_h = (theta * c)/(theta * c + (1-theta) * exp(gamma_h));

		// update c
		double sum_gamma = sum(gamma);
		b_sum_h = gamma[order];
		b_sum_h = b_sum_h * abs(beta[order]);
		b_sum_h = b_sum_h * lambda;
		double b_sum = sum(b_sum_h)/sigma;
		if (sum_gamma > 0) {
			if (b_sum > 0) 
				c = EX_trunc_gamma(sum_gamma, b_sum);	
			else
				c = sum_gamma/(sum_gamma + 1);
			}
		else
			c = 0.5;
	
		// update theta, gamma, w based on previous calculations
		theta = (sum_gamma + a_prior)/(p + a_prior + b_prior);
		//std::copy(gamma_h.begin(), gamma_h.end(), gamma.begin());
		for(int i=0; i<p;++i){
		  gamma[order[i]]=gamma_h[i];
		}
		w = 1.0 - (1.0 - c) * gamma;

		
		// Compute rewieghted SLOPE estimator using computed weights and sigma
		lambda_sigma = lambda * sigma;
		w_vec = as<vec>(w);
		div_X_by_w(X_div_w, X, w_vec, n, p);
		beta_arma = slope_admm(X_div_w, Y, lambda_sigma, p, 1.0);
		for (int i=0; i<p; ++i) {
			beta_arma[i] /= w_vec[i];
		}
		beta_new = as<NumericVector>(wrap(beta_arma));
	        linshrink_cov(X, Sigma, n, p);
    
                S = inv_sympd(Sigma);
                mu = zeros(p);
                for (int i=0; i<p; ++i) {
                    for (int j=0; j<n; ++j) {
                      mu.at(i) += X.at(j,i); 
                    }
                    mu[i] /= n;
                }
    
    
                impute_advance(beta_new, X, Y, S, sigma, n, p, mu, XisFin, anyNanXrows, nanIndicesInRow);
	        X= normalise(X); 
		// Check stop condition
		error= sum(abs(beta-beta_new));
		if (error < tol) {
			iter = max_iter;
			converged = true;
		}
		if(verbose){
		
	 	  Rcout<< "Error =  "<< error <<" sigma = "<< sigma <<" theta = "<< theta<<" c = "<< c<<"\n";
		}
		std::copy(beta_new.begin(), beta_new.end(), beta.begin()) ;
		++iter;
	}


	return List::create(Named("beta")=beta, Named("sigma")=sigma, Named("theta")=theta, Named("c")=c,
                           Named("w")=w, Named("converged")=converged, Named("X")=X,Named("Sigma")=Sigma,Named("mu")=mu);
}
