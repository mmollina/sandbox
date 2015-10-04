#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.00001
#define LN_5 log(1.0/2.0)

// [[Rcpp::export]]
SEXP est_rf_bc(NumericVector x, int n_ind) {
  int ct=0, mis=0, nr=0;
  double rtemp, l, l0;
  NumericMatrix r(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  for(int i=0; i < (int)(x.size()/n_ind)-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < (int)x.size()/n_ind; j++)
        {
	  nr=mis=0;
          for(int k=0; k < n_ind; k++)
            {
              if(x(i*n_ind+k)==1 && x(j*n_ind+k)==2 || x(i*n_ind+k)==2 && x(j*n_ind+k)==1) nr++;
	      else if(x(i*n_ind+k)==0 || x(j*n_ind+k)==0 ) mis++;
	    }
	  rtemp=(double)nr/(double)(n_ind-mis);
	  l=(n_ind-(mis+nr))*log(1-rtemp)+nr*log(rtemp);
	  l0=LN_5*(n_ind-mis);
	  r(j,i)=rtemp;
	  r(i,j)=(l-l0)/log(10.0);
	}
    }
  return wrap(r);
}
