#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.000001

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

Rcpp::NumericMatrix count_genotypes(Rcpp::NumericVector x,
				    int i,
				    int j,
				    int n_ind)
{
  int mis=0;
  NumericMatrix n(5,5);
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
  n(0,0)=mis;
  return(n);
}


Rcpp::NumericVector rf_A_A(Rcpp::NumericMatrix n,
			  int n_ind,
			  int mis)
{
  NumericVector r(8);
  int n1, n2, n3, n4;
  double l0, l;
  l0=-2.0*M_LN2*(n_ind-mis);  
  n1=n(4,1)+n(3,2)+n(2,3)+n(1,4);   
  n2=n(4,4)+n(3,3)+n(2,2)+n(1,1);
  n3=n(4,3)+n(3,4)+n(2,1)+n(1,2);
  n4=n(4,2)+n(3,1)+n(2,4)+n(1,3);
  
  r[0]=(2.0*(n1)+n3+n4)/(2.0*(n_ind-mis));
  l=(n3+n4)*log((1-r[0])*r[0])+2.0*(n1)*log(r[0])+2.0*(n2)*log(1-r[0]);
  r[4]=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  r[1]=(2.0*(n4)+n2+n1)/(2.0*(n_ind-mis));
  l=(n2+n1)*log((1-r[1])*r[1])+2.0*(n4)*log(r[1])+2.0*(n3)*log(1-r[1]);
  r[5]=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  r[2]=(2.0*(n3)+n2+n1)/(2.0*(n_ind-mis));
  l=(n2+n1)*log((1-r[2])*r[2])+2.0*(n3)*log(r[2])+2.0*(n4)*log(1-r[2]);
  r[6]=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  r[3]=(2.0*(n2)+n3+n4)/(2.0*(n_ind-mis));
  l=(n3+n4)*log((1-r[3])*r[3])+2.0*(n2)*log(r[3])+2.0*(n1)*log(1-r[3]);
  r[7]=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/

  return(r);
}

Rcpp::NumericVector rf_A_B3(Rcpp::NumericMatrix n,
			  int n_ind,
			  int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;

  /*Likelihoods under h0: r=0.5*/
  l0= -2.0*M_LN2*(n_ind-mis);

  /*EM algorithm*/
  rold=0, rnew=0.01;

  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(((n(3,2)+n(2,2))*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 	    
	    2.0*(n(1,3)+n(4,1)) + n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,2)+n(2,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))*log(rnew-(rnew*rnew))+
    2.0*(n(1,3)+n(4,1))*log(rnew)+
    2.0*(n(4,3)+n(1,1))*log(1.0-rnew);
  r(4)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;
	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(((n(4,2)+n(1,2))*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    2.0*(n(2,3)+n(3,1)) + n(4,3)+n(1,3)+n(3,2)+n(2,2)+n(4,1)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,2)+n(1,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(4,3)+n(1,3)+n(3,2)+n(2,2)+n(4,1)+n(1,1))*log(rnew-(rnew*rnew))+
    2.0*(n(2,3)+n(3,1))*log(rnew)+
    2.0*(n(3,3)+n(2,1))*log(1.0-rnew);
  r(5)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(((n(4,2)+n(1,2))*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    2.0*(n(3,3)+n(2,1)) + n(4,3)+n(1,3)+n(3,2)+n(2,2)+n(4,1)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(2)=rnew;
  l=(n(4,2)+n(1,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(4,3)+n(1,3)+n(3,2)+n(2,2)+n(4,1)+n(1,1))*log(rnew-(rnew*rnew))+
    2.0*(n(3,3)+n(2,1))*log(rnew)+
    2.0*(n(2,3)+n(3,1))*log(1.0-rnew);
  r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(((n(3,2)+n(2,2))*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    2.0*(n(4,3)+n(1,1)) + n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))/(2.0*(n_ind-mis));
    }
  r(3)=rnew;
  l=(n(3,2)+n(2,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))*log(rnew-(rnew*rnew))+
    2.0*(n(4,3)+n(1,1))*log(rnew)+
    2.0*(n(1,3)+n(4,1))*log(1.0-rnew);
  r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/    
  return(r);
}





