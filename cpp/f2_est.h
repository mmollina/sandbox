#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Rcpp::NumericVector est_rf_C_C(std::vector<int> k_sub,
			       std::vector<int> k1_sub,
			       int n_ind);
Rcpp::NumericVector est_rf_C_D_43(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_C_D_51(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_D_D_43(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_D_D_51(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_D_D_43_51(std::vector<int> k_sub,
				     std::vector<int> k1_sub,
				     int n_ind);
Rcpp::NumericVector est_rf_A_A(std::vector<int> k_sub,
			       std::vector<int> k1_sub,
			       int n_ind);
