#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.00001

// [[Rcpp::export]]
SEXP est_rf_f2_old(NumericVector x, NumericVector n) {
  int ct=0, mis=0, nrec=0, n_ind=n(0), n0=0, n11=0, n12=0, n21=0, n22=0, n31=0, n32=0, n41=0, n42=0, n5=0, n6=0, n7=0, n8=0;
  double q, rold, rnew, r0, r1, r2, l, l0;
  NumericMatrix r(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  for(int i=0; i < (int)(x.size()/n_ind)-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < (int)x.size()/n_ind; j++)
        {
          rold=0, rnew=0.01;
	  n0=n11=n12=n21=n22=n31=n32=n41=n42=n5=n6=n7=n8=0;
	  mis=0;
          for(int k=0; k < n_ind; k++)
            {
              if(x(i*n_ind+k)==1) 
                {
		  if (x(j*n_ind+k)==1) n0++;
		  else if (x(j*n_ind+k)==2) n11++;
		  else if (x(j*n_ind+k)==3) n8++;
		  else if (x(j*n_ind+k)==4) n41++;
		  else if (x(j*n_ind+k)==5) n21++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==2) 
                {
		  if (x(j*n_ind+k)==1) n12++;
		  else if (x(j*n_ind+k)==2) n7++;
		  else if (x(j*n_ind+k)==3) n12++;
		  else if (x(j*n_ind+k)==4) n31++;
		  else if (x(j*n_ind+k)==5) n31++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==3) 
                {
		  if (x(j*n_ind+k)==1) n8++;
		  else if (x(j*n_ind+k)==2) n11++;
		  else if (x(j*n_ind+k)==3) n0++;
		  else if (x(j*n_ind+k)==4) n21++;
		  else if (x(j*n_ind+k)==5) n41++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==4) 
                {
		  if (x(j*n_ind+k)==1) n42++;
		  else if (x(j*n_ind+k)==2) n32++;
		  else if (x(j*n_ind+k)==3) n22++;
		  else if (x(j*n_ind+k)==4) n5++;
		  else if (x(j*n_ind+k)==5) n6++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==5) 
                {
		  if (x(j*n_ind+k)==1) n22++;
		  else if (x(j*n_ind+k)==2) n32++;
		  else if (x(j*n_ind+k)==3) n42++;
		  else if (x(j*n_ind+k)==4) n6++;
		  else if (x(j*n_ind+k)==5) n5++;
		  else mis++;
                }
	      else mis++;
	    }
	  
	  /*EM algorithm*/
	  while(abs(rold-rnew) > TOL)
	    {
	      rold=rnew;
	      r0=(1-rold)*(1-rold);
	      r1=(1-rold)*rold;
	      r2=rold*rold;
	      rnew=((n11 + n12) + 
		    (n21 + n22)*2*r1/(r2+2*r1) +
		    (n31 + n32)*r1/(r1+r0+r2) +
		    (n41 + n42)*2*r1/(r0+2*r1) +
		    n5*4*r1/(3*r0 + 4*r1 + 2*r2) +
		    n6*4*r1/(2*r0 + 4*r1 + 3*r2) +
		    2*(n8 +
		       (n21 + n22)*r2/(r2+2*r1) +
		       n7*r2/(r0+r2) +
		       (n31 + n32)*r2/(r1+r0+r2) +
		       n5*2*r2/(3*r0 + 4*r1 + 2*r2) +
		       n6*3*r2/(2*r0 + 4*r1 + 3*r2)
		       ))/(2*(n_ind-mis));
	    }
	  r0=(1.0-rnew)*(1.0-rnew);
	  r1=rnew*(1.0-rnew);
	  r2=rnew*rnew;
	  /*Likelihood*/
	  l=n0 * (2.0*log(1.0-rnew)) + 
	    n11  *  (M_LN2 + log(rnew) + log(1.0-rnew)) + 
	    n12  *  (log(r1)) + 
	    n21  *  (log(1.0-r0)) + 
	    n22  *  (log((1.0-r0) / 3.0)) + 
	    n31  *  (log(1.0-r1)) + 
	    n32  *  (log((1.0-r1) / 3.0) + M_LN2) + 
	    n41  *  (log(1.0-r2)) + 
	    n42  *  (log((1.0-r2) / 3.0)) + 
	    n5 * (log((r0 + 2.0) / 3.0)) + 
	    n6 * (log((r2+2.0) / 3.0)) + 
	    n7 * (log(r2+r0)) + 
	    n8 * (2.0*log(rnew));
	  l0= -2 * M_LN2 *(n8 + n42 + n22 + n12 + n0) - M_LN2 *(n7 + n32 + n11 ) - 0.28768207245178 * (n6 + n5 + n41 + n31 + n21);
	  r(j,i)=rnew;
	  r(i,j)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
	}
    }
  return wrap(r);
}
