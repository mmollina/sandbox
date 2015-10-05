#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.00001
#define LN_75 -0.28768207245178 

// [[Rcpp::export]]
SEXP est_rf_f2(NumericVector x, int n_ind) {
  //NumericVector n(14);
  int n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13;
  double rold, rnew, r0, r1, r2, l, l0;
  NumericMatrix r(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  for(int i=0; i < (int)(x.size()/n_ind)-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < (int)x.size()/n_ind; j++)
        {
          rold=0, rnew=0.01;
	  n0=n1=n2=n3=n4=n5=n6=n7=n8=n9=n10=n11=n12=n13=0;
	  //std::fill(n.begin(), n.end(), 0);
	  for(int k=0; k < n_ind; k++)
	    {
	      if(x(i*n_ind+k)==1) 
		{
		  if (x(j*n_ind+k)==1) n1++;
		  else if (x(j*n_ind+k)==2) n2++;
		  else if (x(j*n_ind+k)==3) n3++;
		  else if (x(j*n_ind+k)==4) n4++;
		  else if (x(j*n_ind+k)==5) n5++;
		  else n0++;
		}
	      else if(x(i*n_ind+k)==2) 
		{
		  if (x(j*n_ind+k)==1) n6++;
		  else if (x(j*n_ind+k)==2) n7++;
		  else if (x(j*n_ind+k)==3) n6++;
		  else if (x(j*n_ind+k)==4) n8++;
		  else if (x(j*n_ind+k)==5) n8++;
		  else n0++;
		}
	      else if(x(i*n_ind+k)==3) 
		{
		  if (x(j*n_ind+k)==1) n3++;
		  else if (x(j*n_ind+k)==2) n2++;
		  else if (x(j*n_ind+k)==3) n1++;
		  else if (x(j*n_ind+k)==4) n5++;
		  else if (x(j*n_ind+k)==5) n4++;
		  else n0++;
		}
	      else if(x(i*n_ind+k)==4) 
		{
		  if (x(j*n_ind+k)==1) n9++;
		  else if (x(j*n_ind+k)==2) n10++;
		  else if (x(j*n_ind+k)==3) n11++;
		  else if (x(j*n_ind+k)==4) n12++;
		  else if (x(j*n_ind+k)==5) n13++;
		  else n0++;
		}
	      else if(x(i*n_ind+k)==5) 
		{
		  if (x(j*n_ind+k)==1) n11++;
		  else if (x(j*n_ind+k)==2) n10++;
		  else if (x(j*n_ind+k)==3) n9++;
		  else if (x(j*n_ind+k)==4) n13++;
		  else if (x(j*n_ind+k)==5) n12++;
		  else n0++;
		}
	      else n0++;
	    }
	  //n=count_genotypes_f2( x, i, j, n_ind);
	  //EM algorithm
	  while(abs(rold-rnew) > TOL)
	    {
	      rold=rnew;
	      r0=(1-rold)*(1-rold);
	      r1=(1-rold)*rold;
	      r2=rold*rold;
	      rnew=((n2 + n6)+ 
		    (n5 + n11)*2*r1/(r2+2*r1) +
		    (n8 + n10)*r1/(r1+r0+r2) +
		    (n4 + n9)*2*r1/(r0+2*r1) + 
		    n12*4*r1/(3*r0 + 4*r1 + 2*r2) +
		    n13*4*r1/(2*r0 + 4*r1 + 3*r2) +
		    2*(n3 +
		       (n5 + n11)*r2/(r2+2*r1) +
		       n7*r2/(r0+r2) +
		       (n8 + n10)*r2/(r1+r0+r2) + 
		       n12*2*r2/(3*r0 + 4*r1 + 2*r2) +
		       n13*3*r2/(2*r0 + 4*r1 + 3*r2)
		       ))/(2*(n_ind-n0));
	    }
	  r0=(1.0-rnew)*(1.0-rnew);
	  r1=rnew*(1.0-rnew);
	  r2=rnew*rnew;
	  //Likelihood
	  l=n1 * (2.0*log(1.0-rnew)) + 
	    n2  *  (M_LN2 + log(rnew) + log(1.0-rnew)) + 
	    n6  *  (log(r1)) + 
	    n5  *  (log(1.0-r0)) + 
	    n11  *  (log((1.0-r0) / 3.0)) + 
	    n8  *  (log(1.0-r1)) + 
	    n10  *  (log((1.0-r1) / 3.0) + M_LN2) + 
	    n4  *  (log(1.0-r2)) + 
	    n9  *  (log((1.0-r2) / 3.0)) + 
	    n12 * (log((r0 + 2.0) / 3.0)) + 
	    n13 * (log((r2+2.0) / 3.0)) + 
	    n7 * (log(r2+r0)) + 
	    n3 * (2.0*log(rnew));
	  //Likelihood unde H0: r=0.5
	  l0= -2 * M_LN2 *(n3 + n9 + n11 + n6 + n1) - M_LN2 *(n7 + n10 + n2 ) +
	    LN_75*(n13 + n12 + n4 + n8 + n5);
	  r(j,i)=rnew;
	  r(i,j)=(l-l0)/log(10.0); //transforming to base 10 logarithm
	}
    }
  return wrap(r);
}









