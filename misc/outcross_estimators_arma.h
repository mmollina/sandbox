#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
#define TOL 1e-05

arma::mat Tr(double r, int ph);

arma::mat Nr(int ph);

arma::mat count_genotypes(Rcpp::NumericVector x,
			    int i,
			    int j,
			    int n_ind);

arma::mat Mt(int type);

Rcpp::NumericVector est_rf_pair(arma::mat n,
				  int n_mar,
                                  int n_ind,
				  int mis,
				  int sg1,
				  int sg2);
