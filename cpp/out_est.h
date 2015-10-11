#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;
using namespace std;

Rcpp::NumericVector rf_A_A(Rcpp::NumericMatrix n,
			     int n_ind,
			   int mis);
Rcpp::NumericVector rf_A_B1(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_B2(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_B3(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_C(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis);
Rcpp::NumericVector rf_A_D1(Rcpp::NumericMatrix n,
			   int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_D2(Rcpp::NumericMatrix n,
			   int n_ind,
			    int mis);
Rcpp::NumericVector rf_B1_B1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_C(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_B2_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B3_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B3_C(Rcpp::NumericMatrix n,
			     int n_ind,
			    int mis);
Rcpp::NumericVector rf_B3_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B3_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_C_C(Rcpp::NumericMatrix n,
			     int n_ind,
			   int mis);
Rcpp::NumericVector rf_C_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			    int mis);
Rcpp::NumericVector rf_C_D2(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_D1_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_D2_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
