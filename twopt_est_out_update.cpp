#include <Rcpp.h>
#include "outcross_estimators.h"
using namespace Rcpp;
using namespace std;
#define TOL 0.000001

// [[Rcpp::export]]
SEXP est_rf_out_new(NumericVector x, NumericVector segreg_type, NumericVector n_ind_R) {
  int ct=0, nrec=0, n_ind=n_ind_R(0);
  int n1, n2, n3, n4, n_mar=((int)x.size()/n_ind), k1, k2;
  double q, rold, rnew, l1, l2, l3, l4, l0, l01, l02;
  NumericMatrix n(5,5);
  NumericVector r(8);
  NumericMatrix r1(n_mar,n_mar);  
  NumericMatrix r2(n_mar,n_mar);
  NumericMatrix r3(n_mar,n_mar);
  NumericMatrix r4(n_mar,n_mar);
  for(int i=0; i < n_mar-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < n_mar; j++)
        {
	  n=count_genotypes( x, i, j, n_ind);
	  std::fill(r.begin(), r.end(), 0);
	  k1=segreg_type(i); k2=segreg_type(j);
	  switch(k1){
	  case 1:
	    switch(k2){
	    case 1: 
	      r=rf_A_A(n, n_ind, n(0,0)); 	      /*Markers A - A */
	      break;
	    case 2:
	      r=rf_A_B1(n, n_ind, n(0,0));	      /*Markers A - B1 */
	      break;
	    case 3:
	      r=rf_A_B2(n, n_ind, n(0,0));	      /*Markers A - B2 */
	      break;
	    case 4:
	      r=rf_A_B3(n, n_ind, n(0,0));	      /*Markers A - B3 */
	      break;
	    case 5:
	      r=rf_A_C(n, n_ind, n(0,0));	      /*Markers A - C */
	      break;
	    case 6:
	      r=rf_A_D1(n, n_ind, n(0,0));	      /*Markers A - D1 */
	      break;
	    case 7:
	      r=rf_A_D2(n, n_ind, n(0,0));	      /*Markers A - D2 */
	      break;
	    }
	    break;
	  case 2:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_B1(n,n_ind, n(0,0));	      /*Markers B1 - A*/
	      break;
	    case 2: 
	      r=rf_B1_B1(n,n_ind, n(0,0));	      /*Markers B1 - B1*/
	      break;	    
	    case 3: 
	      r=rf_B1_B2(n,n_ind, n(0,0));	      /*Markers B1 - B2*/
	      break;	    
	    case 4: 
	      r=rf_B1_B3(n,n_ind, n(0,0));	      /*Markers B1 - B3*/
	      break;	    
	    case 5: 
	      r=rf_B1_C(n,n_ind, n(0,0));	      /*Markers B1 - c*/
	      break;	    
	    }
	    break;
	  case 3:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_B2(n,n_ind, n(0,0));	      /*Markers B2 - A*/
	      break;	    
	    case 2: 
	      n=transpose_counts(n);
	      r=rf_B1_B2(n,n_ind, n(0,0));	      /*Markers B2 - B1*/
	      break;	    

	    }
	    break;
	  case 4:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_B3(n,n_ind, n(0,0));	      /*Markers B3 - A*/
	      break;	    
	    case 2: 
	      n=transpose_counts(n);
	      r=rf_B1_B3(n,n_ind, n(0,0));	      /*Markers B3 - B1*/
	      break;	    	      
	    }
	    break;
	  case 5:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_C(n,n_ind, n(0,0));	      /*Markers C - A*/
	      break;	    
	    case 2: 
	      n=transpose_counts(n);
	      r=rf_B1_C(n,n_ind, n(0,0));	      /*Markers C - B1*/
	      break;	    
	    }
	    break;
	  case 6:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_D1(n,n_ind, n(0,0));	      /*Markers D1 - A*/
	      break;	    
	    }
	    break;
	  case 7:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_D2(n,n_ind, n(0,0));	      /*Markers D2 - A*/
	      break;	    
	    }
	    break;
	  }
	  r1(j,i)=r[0];
	  r2(j,i)=r[1];
	  r3(j,i)=r[2];
	  r4(j,i)=r[3];
	  r1(i,j)=r[4];
	  r2(i,j)=r[5];
	  r3(i,j)=r[6];
	  r4(i,j)=r[7];
	}
    }
  List z = List::create(wrap(r1), wrap(r2), wrap(r3), wrap(r4));
  return(z);
}
