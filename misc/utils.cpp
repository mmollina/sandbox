#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Rcpp::NumericMatrix transpose_counts(Rcpp::NumericMatrix n)
{
  int temp;
  temp=n(1,2); n(1,2)=n(2,1); n(2,1)=temp;
  temp=n(1,3); n(1,3)=n(3,1); n(3,1)=temp;
  temp=n(1,4); n(1,4)=n(4,1); n(4,1)=temp;
  temp=n(2,3); n(2,3)=n(3,2); n(3,2)=temp;
  temp=n(2,4); n(2,4)=n(4,2); n(4,2)=temp;
  temp=n(3,4); n(3,4)=n(4,3); n(4,3)=temp;
  return(n);
}

Rcpp::NumericVector count_genotypes_f2(Rcpp::NumericVector x,
				       int i,
				       int j,
				       int n_ind)
{
  NumericVector n(14);
  std::fill(n.begin(), n.end(), 0);
   for(int k=0; k < n_ind; k++)
    {
      if(x(i*n_ind+k)==1) 
	{
	  if (x(j*n_ind+k)==1) n(1)++;
	  else if (x(j*n_ind+k)==2) n(2)++;
	  else if (x(j*n_ind+k)==3) n(3)++;
	  else if (x(j*n_ind+k)==4) n(4)++;
	  else if (x(j*n_ind+k)==5) n(5)++;
	  else n(0)++;
	}
      else if(x(i*n_ind+k)==2) 
	{
	  if (x(j*n_ind+k)==1) n(6)++;
	  else if (x(j*n_ind+k)==2) n(7)++;
	  else if (x(j*n_ind+k)==3) n(6)++;
	  else if (x(j*n_ind+k)==4) n(8)++;
	  else if (x(j*n_ind+k)==5) n(8)++;
	  else n(0)++;
	}
      else if(x(i*n_ind+k)==3) 
	{
	  if (x(j*n_ind+k)==1) n(3)++;
	  else if (x(j*n_ind+k)==2) n(2)++;
	  else if (x(j*n_ind+k)==3) n(1)++;
	  else if (x(j*n_ind+k)==4) n(5)++;
	  else if (x(j*n_ind+k)==5) n(4)++;
	  else n(0)++;
	}
      else if(x(i*n_ind+k)==4) 
	{
	  if (x(j*n_ind+k)==1) n(9)++;
	  else if (x(j*n_ind+k)==2) n(10)++;
	  else if (x(j*n_ind+k)==3) n(11)++;
	  else if (x(j*n_ind+k)==4) n(12)++;
	  else if (x(j*n_ind+k)==5) n(13)++;
	  else n(0)++;
	}
      else if(x(i*n_ind+k)==5) 
	{
	  if (x(j*n_ind+k)==1) n(11)++;
	  else if (x(j*n_ind+k)==2) n(10)++;
	  else if (x(j*n_ind+k)==3) n(9)++;
	  else if (x(j*n_ind+k)==4) n(13)++;
	  else if (x(j*n_ind+k)==5) n(12)++;
	  else n(0)++;
	}
      else n(0)++;
    }
  return(n);
}

Rcpp::NumericMatrix count_genotypes_outcross(Rcpp::NumericVector x,
					     int i,
					     int j,
					     int n_ind)
{
  NumericMatrix n(5,5);
  std::fill(n.begin(), n.end(), 0);
  for(int k=0; k < n_ind; k++)
    {
      if(x(i*n_ind+k)==1) 
	{
	  if (x(j*n_ind+k)==1) n(1,1)++;
	  else if (x(j*n_ind+k)==2) n(1,2)++;
	  else if (x(j*n_ind+k)==3) n(1,3)++;
	  else if (x(j*n_ind+k)==4) n(1,4)++;
	  else n(0,0)++;
	}
      else if(x(i*n_ind+k)==2) 
	{
	  if (x(j*n_ind+k)==1) n(2,1)++;
	  else if (x(j*n_ind+k)==2) n(2,2)++;
	  else if (x(j*n_ind+k)==3) n(2,3)++;
	  else if (x(j*n_ind+k)==4) n(2,4)++;
	  else n(0,0)++;
	}
      else if(x(i*n_ind+k)==3) 
	{
	  if (x(j*n_ind+k)==1) n(3,1)++;
	  else if (x(j*n_ind+k)==2) n(3,2)++;
	  else if (x(j*n_ind+k)==3) n(3,3)++;
	  else if (x(j*n_ind+k)==4) n(3,4)++;
	  else n(0,0)++;
	}
      else if(x(i*n_ind+k)==4) 
	{
	  if (x(j*n_ind+k)==1) n(4,1)++;
	  else if (x(j*n_ind+k)==2) n(4,2)++;
	  else if (x(j*n_ind+k)==3) n(4,3)++;
	  else if (x(j*n_ind+k)==4) n(4,4)++;
	  else n(0,0)++;
	}
      else n(0,0)++;
    }
  return(n);
}
