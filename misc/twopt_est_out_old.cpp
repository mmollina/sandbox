#include <Rcpp.h>
#include "utils.h"
#include "outcross_estimators.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP est_rf_out_old(NumericVector x, NumericVector segreg_type, int n_ind) {
  int n_mar=((int)x.size()/n_ind), k1, k2;
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
	  std::fill(r.begin(), r.end(), 0);
	  n=count_genotypes_outcross( x, i, j, n_ind);
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
	    case 6: 
	      r=rf_B1_D1(n,n_ind, n(0,0));	      /*Markers B1 - D1*/
	      break;	    
	    case 7: 
	      r=rf_B1_D2(n,n_ind, n(0,0));	      /*Markers B1 - D2*/
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
	    case 3: 
	      r=rf_B2_B2(n,n_ind, n(0,0));	      /*Markers B2 - B2*/
	      break;	    
	    case 4: 
	      r=rf_B2_B3(n,n_ind, n(0,0));	      /*Markers B2 - B3*/
	      break;	    
	    case 5: 
	      r=rf_B2_C(n,n_ind, n(0,0));	      /*Markers B2 - C*/
	      break;
	    case 6: 
	      r=rf_B2_D1(n,n_ind, n(0,0));	      /*Markers B2 - D1*/
	      break;
	    case 7: 
	      r=rf_B2_D2(n,n_ind, n(0,0));	      /*Markers B2 - D2*/
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
	    case 3: 
	      n=transpose_counts(n);
	      r=rf_B2_B3(n,n_ind, n(0,0));	      /*Markers B3 - B2*/
	      break;	    	      
	    case 4: 
	      r=rf_B3_B3(n,n_ind, n(0,0));	      /*Markers B3 - B3*/
	      break;	    	      
	    case 5: 
	      r=rf_B3_C(n,n_ind, n(0,0));	      /*Markers B3 - C*/
	      break;
	    case 6: 
	      r=rf_B3_D1(n,n_ind, n(0,0));	      /*Markers B3 - D1*/
	      break;
	    case 7: 
	      r=rf_B3_D2(n,n_ind, n(0,0));	      /*Markers B3 - D2*/
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
	    case 3: 
	      n=transpose_counts(n);
	      r=rf_B2_C(n,n_ind, n(0,0));	      /*Markers C - B2*/
	      break;	    
	    case 4: 
	      n=transpose_counts(n);
	      r=rf_B3_C(n,n_ind, n(0,0));	      /*Markers C - B3*/
	      break;	    
	    case 5: 
	      r=rf_C_C(n,n_ind, n(0,0));	      /*Markers C - C*/
	      break;	    
	    case 6: 
	      r=rf_C_D1(n,n_ind, n(0,0));	      /*Markers C - D1*/
	      break;
	    case 7: 
	      r=rf_C_D2(n,n_ind, n(0,0));	      /*Markers C - D2*/
	      break;
	    }
	    break;
	  case 6:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_D1(n,n_ind, n(0,0));	      /*Markers D1 - A*/
	      break;	    
	    case 2: 
	      n=transpose_counts(n);
	      r=rf_B1_D1(n,n_ind, n(0,0));	      /*Markers D1 - B1*/
	      break;	    
	    case 3: 
	      n=transpose_counts(n);
	      r=rf_B2_D1(n,n_ind, n(0,0));	      /*Markers D1 - B1*/
	      break;	    
	    case 4: 
	      n=transpose_counts(n);
	      r=rf_B3_D1(n,n_ind, n(0,0));	      /*Markers D1 - B3*/
	      break;	    
	    case 5: 
	      n=transpose_counts(n);
	      r=rf_C_D1(n,n_ind, n(0,0));	      /*Markers D1 - C*/
	      break;	    
	    case 6: 
	      r=rf_D1_D1(n,n_ind, n(0,0));	      /*Markers D1 - D1*/
	      break;	   
	    case 7: 
	      r= rep( NumericVector::get_na(), 8 );   /*Markers D1 - D2 - Impossible to compute*/
	      break;	    
	    }
	    break;
	  case 7:
	    switch(k2){
	    case 1: 
	      n=transpose_counts(n);
	      r=rf_A_D2(n,n_ind, n(0,0));	      /*Markers D2 - A*/
	      break;	    
	    case 2: 
	      n=transpose_counts(n);
	      r=rf_B1_D2(n,n_ind, n(0,0));	      /*Markers D2 -  B1*/
	      break;
	    case 3: 
	      n=transpose_counts(n);
	      r=rf_B1_D2(n,n_ind, n(0,0));	      /*Markers D2 -  B2*/
	      break;	    
	    case 4: 
	      n=transpose_counts(n);
	      r=rf_B3_D2(n,n_ind, n(0,0));	      /*Markers D2 -  B3*/
	      break;	    
	    case 5: 
	      n=transpose_counts(n);
	      r=rf_C_D2(n,n_ind, n(0,0));	      /*Markers D2 -  C*/
	      break;	
	    case 6: 
	      r= rep( NumericVector::get_na(), 8 );   /*Markers D2 -  D1 - Impossible to compute*/
	      break;	
	    case 7: 
	      n=transpose_counts(n);
	      r=rf_D2_D2(n,n_ind, n(0,0));	      /*Markers D2 -  D2*/
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
