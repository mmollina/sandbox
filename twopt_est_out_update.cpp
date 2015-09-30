#include <Rcpp.h>
#include "outcross_estimators.h"
using namespace Rcpp;
using namespace std;
#define TOL 0.000001

// [[Rcpp::export]]
SEXP est_rf_out(NumericVector x, NumericVector segreg_type, NumericVector n_ind_R) {
  int ct=0, mis=0, nrec=0, n_ind=n_ind_R(0);
  int n1, n2, n3, n4, temp, n_mar=((int)x.size()/n_ind);
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
	  std::fill(n.begin(), n.end(), 0);
	  mis=0;
          for(int k=0; k < n_ind; k++)
            {
              if(x(i*n_ind+k)==1) 
                {
		  if (x(j*n_ind+k)==1) n(1,1)++;
		  else if (x(j*n_ind+k)==2) n(1,2)++;
		  else if (x(j*n_ind+k)==3) n(1,3)++;
		  else if (x(j*n_ind+k)==4) n(1,4)++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==2) 
                {
		  if (x(j*n_ind+k)==1) n(2,1)++;
		  else if (x(j*n_ind+k)==2) n(2,2)++;
		  else if (x(j*n_ind+k)==3) n(2,3)++;
		  else if (x(j*n_ind+k)==4) n(2,4)++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==3) 
                {
		  if (x(j*n_ind+k)==1) n(3,1)++;
		  else if (x(j*n_ind+k)==2) n(3,2)++;
		  else if (x(j*n_ind+k)==3) n(3,3)++;
		  else if (x(j*n_ind+k)==4) n(3,4)++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==4) 
                {
		  if (x(j*n_ind+k)==1) n(4,1)++;
		  else if (x(j*n_ind+k)==2) n(4,2)++;
		  else if (x(j*n_ind+k)==3) n(4,3)++;
		  else if (x(j*n_ind+k)==4) n(4,4)++;
		  else mis++;
                }
	      else mis++;
	    }
	  
	  if(segreg_type(i) == 1 && segreg_type(j) == 1)
	    { 
	      r=rf_A_A(n,n_ind, mis);
	      r1(j,i)=r[0];
	      r2(j,i)=r[1];
	      r3(j,i)=r[2];
	      r4(j,i)=r[3];
	      r1(i,j)=r[4];
	      r2(i,j)=r[5];
	      r3(i,j)=r[6];
	      r4(i,j)=r[7];
	    }

	  if((segreg_type(i) == 4 && segreg_type(j) == 1) ||(segreg_type(i) == 1 && segreg_type(j) == 4) )
	    { 
	      if(segreg_type(i) == 1 && segreg_type(j) == 4)
		{
		  temp=n(1,2); n(1,2)=n(2,1); n(2,1)=temp;
		  temp=n(1,3); n(1,3)=n(3,1); n(3,1)=temp;
		  temp=n(1,4); n(1,4)=n(4,1); n(4,1)=temp;
		  temp=n(2,3); n(2,3)=n(3,2); n(3,2)=temp;
		  temp=n(2,4); n(2,4)=n(4,2); n(4,2)=temp;
		  temp=n(3,4); n(3,4)=n(4,3); n(4,3)=temp;
		}
	      r=rf_A_B3(n,n_ind, mis);
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
    }
  List z = List::create(wrap(r1), wrap(r2), wrap(r3), wrap(r4));
  return(z);





}



 





