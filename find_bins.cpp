#include <Rcpp.h>
using namespace Rcpp;

int check_occurrence(std::vector<int>& v, int x)
{
    for(int i = 0; i < v.size(); i++)
        if(v[i]==x)
            return(i);
    return(-1);
}

// [[Rcpp::export]]
IntegerVector getbins(IntegerMatrix geno){ 
  int n_mar = geno.ncol();
  int n_ind = geno.nrow();
  std::vector<int> b_vec(1);
  b_vec[0]=1;
  std::vector<int> b(n_mar);
  std::fill(b.begin(), b.end(), 1);
  int flag, l;
  for(int i = 0; i < n_mar; i++)
    {
      for(int j = b_vec.size() - 1; j >= 0 ; j--)
	{
	  flag=0;
	  l=check_occurrence(b, b_vec[j]);
	  for(int k = 0; k < n_ind; k++)
	    {
	      if(geno(k,i)!=geno(k,l) && geno(k,i)!=0 && geno(k,l)!=0)
		{
		  flag=1;
		  break;
		}
	    }
	  if(flag==0)
	    {
	      b[i]=b_vec[j];
	      break;
	    }
	  }
	if(flag==1)
	  {
	    b[i] = *std::max_element(b_vec.begin(), b_vec.end())+1;
	    b_vec.push_back(b[i]);
	  }
	Rcpp::Rcout.precision(2);
	Rcpp::Rcout << std::fixed << 100.0*i/n_mar << " ";
	if(i % 20 == 0)
	  Rcpp::Rcout << "\n";
      }
  return(wrap(b));
}


