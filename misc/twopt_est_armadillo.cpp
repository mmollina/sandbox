#include <RcppArmadillo.h>
#include "outcross_estimators_arma.h"
using namespace Rcpp;
using namespace arma;
using namespace std;
#define TOL 1e-05

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP est_rf_arma(NumericVector x, NumericVector segreg_type, int n_ind) {
  int n_mar=((int)x.size()/n_ind);
  NumericVector r(8);
  NumericMatrix r1(n_mar,n_mar);  
  NumericMatrix r2(n_mar,n_mar);
  NumericMatrix r3(n_mar,n_mar);
  NumericMatrix r4(n_mar,n_mar);
  arma::mat n;
  for(int i=0; i < n_mar-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < n_mar; j++)
        {
	  std::fill(r.begin(), r.end(), 0);
	  n=count_genotypes(x, i, j, n_ind);
	  if(segreg_type(i) < segreg_type(j))
	    r=est_rf_pair(n, n_mar, n_ind, n(0,0), segreg_type(i), segreg_type(j));
	  else if((segreg_type(i) > segreg_type(j)))
	    r=est_rf_pair(trans(n), n_mar, n_ind, n(0,0), segreg_type(j), segreg_type(i));
	  r1(j,i)=r[0];
	  r1(i,j)=r[1];
	  r2(j,i)=r[2];
	  r2(i,j)=r[3];
	  r3(j,i)=r[4];
	  r3(i,j)=r[5];
	  r4(j,i)=r[6];
	  r4(i,j)=r[7];
	}
    } 
  List z = List::create(wrap(r1), wrap(r2), wrap(r3), wrap(r4));
  return(z);
}
