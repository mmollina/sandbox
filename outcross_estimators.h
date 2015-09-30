#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.000001

Rcpp::NumericVector rf_A_A(Rcpp::NumericMatrix n,
			     int n_ind,
			   int mis);
Rcpp::NumericVector rf_A_B3(Rcpp::NumericMatrix n,
			  int n_ind,
			      int mis);
