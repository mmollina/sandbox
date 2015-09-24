#include <Rcpp.h>
using namespace Rcpp;
using namespace std;



/* FUNCTION: nChoosek
 -----------------------------------------------------
 The famous binomial coefficient
 */

int nChoosek(int n, int k)
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  int result = n;
  for( int i = 2; i <= k ; ++i )
  {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}


// [[Rcpp::export]]
SEXP est_rf(NumericVector x) {
  int ct=0, n_ind=250;
  double q, rold, rnew;
  //NumericVector r(nChoosek((int)x.size()/250, 2));
  NumericMatrix r(((int)x.size()/250), ((int)x.size()/250));
  NumericVector n(4);
  for(int i=0; i < (int)(x.size()/250)-1; i++)
    {
//      for(int j=(i+1); j  < (int)x.size()/250; j++)
  for(int j=(i+1); j  < (int)x.size()/250; j++)
        {
          rold=0, rnew=0.03;
          n[1]=n[2]=n[3]=0;
          for(int k=0; k < n_ind; k++)
            {
              if(x(i*n_ind+k)==1) 
                {
		  if (x(j*n_ind+k)==1) n[0]++;
		  else if (x(j*n_ind+k)==2) n[1]++;
		  else if(x(j*n_ind+k)==3) n[2]++;
		  else if(x(j*n_ind+k)==4) n[3]++;
		  else if(x(j*n_ind+k)==5) n[4]++;
                }
              else if(x(i*n_ind+k)==2) 
                {
		  if (x(j*n_ind+k)==1) n[5]++;
		  else if (x(j*n_ind+k)==2) n[6]++;
		  else if(x(j*n_ind+k)==3) n[7]++;
		  else if(x(j*n_ind+k)==4) n[8]++;
		  else if(x(j*n_ind+k)==5) n[9]++;
                }
              else if(x(i*n_ind+k)==3) 
                {
		  if (x(j*n_ind+k)==1) n[10]++;
		  else if (x(j*n_ind+k)==2) n[11]++;
		  else if(x(j*n_ind+k)==3) n[12]++;
		  else if(x(j*n_ind+k)==4) n[13]++;
		  else if(x(j*n_ind+k)==5) n[14]++;
                }
              else if(x(i*n_ind+k)==4) 
                {
		  if (x(j*n_ind+k)==1) n[15]++;
		  else if (x(j*n_ind+k)==2) n[16]++;
		  else if(x(j*n_ind+k)==3) n[17]++;
		  else if(x(j*n_ind+k)==4) n[18]++;
		  else if(x(j*n_ind+k)==5) n[19]++;
                }
              else if(x(i*n_ind+k)==5) 
                {
		  if (x(j*n_ind+k)==1) n[20]++;
		  else if (x(j*n_ind+k)==2) n[21]++;
		  else if(x(j*n_ind+k)==3) n[22]++;
		  else if(x(j*n_ind+k)==4) n[23]++;
		  else if(x(j*n_ind+k)==5) n[24]++;
                }
              }
          while(abs(rold-rnew) > 0.000000001)
            { 
            rold=rnew;
            q=pow(rold,2.0)/(2.0*(pow(rold,2.0)/2+pow((1.0-rold),2.0)/2.0));
            rnew=(n[1] + 2.0*(n[2] + q*n[3]))/(2.0*(double)n_ind);
            }
            r(j,i)=r(i,j)=rnew;
          ct++;
        }
    }
  return wrap(r);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
