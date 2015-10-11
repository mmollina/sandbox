#include <Rcpp.h>
#include "f2_est.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP est_rf_f2(NumericVector x, NumericVector type, int n_ind) {
  int n_mar=((int)x.size()/n_ind), k1, k2;
  NumericMatrix r(n_mar, n_mar);
  Rcpp::NumericVector rtemp(2);   
  for(int i=0; i < n_mar-1; i++)
    {
      R_CheckUserInterrupt(); // check for ^C 
      for(int j=(i+1); j  < n_mar; j++)
        {
	  std::vector<int> k_sub(&x[i*n_ind],&x[i*n_ind+n_ind]);
	  std::vector<int> k1_sub(&x[j*n_ind],&x[j*n_ind+n_ind]);
	  k1=type(i); k2=type(j);
	  // Rcpp::Rcout << k1 << "--" << k2 << "\n";
	  switch(k1){
	  case 1: //C
	    switch(k2){
	    case 1: 
	      rtemp = est_rf_C_C(k_sub, k1_sub, n_ind);  //Markers: Codominant - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_C_D_43(k_sub, k1_sub, n_ind);  //Markers: Codominant - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_C_D_51(k_sub, k1_sub, n_ind);  //Markers Codominant - Dominant( (1, (23)=5) 
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 2: //D_43
	    switch(k2){
	    case 1: 
	      rtemp = est_rf_C_D_43(k1_sub, k_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Codominant
	      break;  
	    case 2:
	      rtemp = est_rf_D_D_43(k_sub, k1_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_D_D_43_51(k_sub, k1_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Dominant( (1, (23)=5) 
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 3: //D51
	    switch(k2){
	    case 1: 
	      rtemp = est_rf_C_D_51(k1_sub, k_sub, n_ind); //Markers  Dominant( (1, (23)=5) - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_D_D_43_51(k1_sub, k_sub, n_ind); //Markers:  Dominant( (1, (23)=5) - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_D_D_51(k_sub, k1_sub, n_ind); //Markers:  Dominant( (1, (23)=5) - Dominant( (1, (23)=5)
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 4:
	    switch(k2){
	    case 1: case 2: case 3: case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind); //Any type of markers - slower function
	      break;
	    }
	    break;
	  }
	  r(j,i)=rtemp(0);
	  r(i,j)=rtemp(1);
	  rtemp(0)=rtemp(1)=0;
	}
    }
  return wrap(r);
}
