#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 1e-05
#define LN3 1.098612288668109
#define LN4 1.38629436111989
#define LN_75 -0.28768207245178

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
  r(0)=(2.0*(n1)+n3+n4)/(2.0*(n_ind-mis));
  l=(n3+n4)*log((1-r(0))*r(0))+2.0*(n1)*log(r(0))+2.0*(n2)*log(1-r(0));
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/ 
  r(1)=(2.0*(n4)+n2+n1)/(2.0*(n_ind-mis));
  l=(n2+n1)*log((1-r(1))*r(1))+2.0*(n4)*log(r(1))+2.0*(n3)*log(1-r(1));
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));

  r(3)=abs(1.0-r(0));
  
  return(r);
}
Rcpp::NumericVector rf_A_B1(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = - (M_LN2 * (n(2,1)+n(1,1) + n(4,1)+n(3,1)+2*n(2,2)+2*n(1,3)) +
	  2 * M_LN2*(n(4,2)+n(3,3)+n(2,3)+n(1,2)+n(4,3)+n(3,2)));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;

  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) +
	    2*(n(2,2)+n(1,3)) + 
	    n(4,2)+n(4,1)+n(3,3)+n(3,1)+n(2,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(4,2)+n(3,3)+n(2,3)+n(1,2))*log(rnew-(rnew*rnew))+(n(4,1)+n(3,1)+2*n(2,2)+2*n(1,3))*log(rnew)+(2*n(4,3)+2*n(3,2)+n(2,1)+n(1,1))*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) + 
	    2*(n(2,3)+n(1,2))+
	    n(4,3)+n(4,1)+n(3,2)+n(3,1)+n(2,2)+n(1,3))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,3)+n(3,2)+n(2,2)+n(1,3))*log(rnew-(rnew*rnew))+(n(4,1)+n(3,1)+2*n(2,3)+2*n(1,2))*log(rnew)+(2*n(4,2)+2*n(3,3)+n(2,1)+n(1,1))*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_A_B2(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = - (M_LN2 * (n(2,1)+n(1,1) + n(4,1)+n(3,1)+2*n(2,2)+2*n(1,3)) +
	  2 * M_LN2*(n(4,2)+n(3,3)+n(2,3)+n(1,2)+n(4,3)+n(3,2)));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;

  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) +
	    2*(n(3,2)+n(1,3)) + 
	    n(4,2)+n(4,1)+n(3,3)+n(2,1)+n(2,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(4,2)+n(3,3)+n(2,3)+n(1,2))*log(rnew-(rnew*rnew))+(n(4,1)+2*n(3,2)+n(2,1)+2*n(1,3))*log(rnew)+(2*n(4,3)+2*n(2,2)+n(3,1)+n(1,1))*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) + 
	    2*(n(2,3)+n(4,2))+
	    n(4,3)+n(3,2)+n(3,1)+n(2,2)+n(1,3)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,3)+n(3,2)+n(2,2)+n(1,3))*log(rnew-(rnew*rnew))+(2*n(4,2)+n(3,1)+2*n(2,3)+n(1,1))*log(rnew)+(n(4,1)+2*n(3,3)+n(2,1)+2*n(1,2))*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
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
	    2.0*(n(1,3)+n(4,1))+n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,2)+n(2,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))*log(rnew-(rnew*rnew))+
    2.0*(n(1,3)+n(4,1))*log(rnew)+
    2.0*(n(4,3)+n(1,1))*log(1.0-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
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
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/

  r(2)=abs(1.0-r(1));

  r(3)=abs(1.0-r(0));  

  return(r);
}

Rcpp::NumericVector rf_A_C(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis)
{
  NumericVector r(8);
  double l, l0, r0, r1, r2, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0=-2*M_LN2*n(4,2)-(LN4-LN3)*n(4,1)-2*M_LN2*(n(3,2)+n(2,2))-(LN4-LN3)*(n(3,1)+n(2,1))-2*M_LN2*n(1,2)-(LN4-LN3)*n(1,1);
  //EM algorithm
  rold=0, rnew=0.01;
  
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((2.0*n(1,1)*r1)/(2.0*r1+r0) +
	    (n(4,1)*(2.0*r2+2.0*r1))/(r2+2.0*r1) +
	    ((n(3,1)+n(2,1))*(2.0*r2+r1))/(r2+r1+r0) + 
	    n(3,2)+n(2,2)+2.0*n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,1)+n(2,1))*log(rnew*rnew-rnew+1)+n(4,1)*log(2.0*rnew-rnew*rnew)+(n(3,2)+n(2,2))*log(rnew-rnew*rnew)+n(1,1)*log(1-rnew*rnew)+2.0*n(1,2)*log(rnew)+2.0*n(4,2)*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((2.0*n(2,1)*r1)/(2.0*r1+r0) +
	    (n(3,1)*(2.0*r2+2.0*r1))/(r2+2.0*r1) +
	    ((n(4,1)+n(1,1))*(2.0*r2+r1))/(r2+r1+r0) + 
	    n(4,2)+2.0*n(2,2)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,1)+n(1,1))*log(rnew*rnew-rnew+1)+n(3,1)*log(2.0*rnew-rnew*rnew)+(n(4,2)+n(1,2))*log(rnew-rnew*rnew)+n(2,1)*log(1-rnew*rnew)+2.0*n(2,2)*log(rnew)+2.0*n(3,2)*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm

  r(2)=abs(1.0-r(1));

  r(3)=abs(1.0-r(0));  

  return(r);
} 
Rcpp::NumericVector rf_A_D1(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis)
{
  NumericVector r(8);
  double l, l0,  rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0, rnew=0.01;
  
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,2)+n(3,2)+n(2,1)+n(1,1))+(rold+1)*(n(4,1)+n(3,1)+n(2,2)+n(1,2)))/(2.0*(n_ind-mis));
    }
  r(1)=r(0)=rnew;
  l=(n(4,1)+n(3,1)+n(2,2)+n(1,2))*log(rnew)+(n(4,2)+n(3,2)+n(2,1)+n(1,1))*log(1-rnew);
  r(4)=r(5)=r(6)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(3)=r(2)=abs(1.0-r(1));
  return(r);
}
Rcpp::NumericVector rf_A_D2(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis)
{
  NumericVector r(8);
  double l, l0,  rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0, rnew=0.01;
  
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,2)+n(3,1)+n(2,2)+n(1,1))+(rold+1)*(n(4,1)+n(3,2)+n(2,1)+n(1,2)))/(2.0*(n_ind-mis));
    }
  r(2)=r(0)=rnew;
  l=(n(4,1)+n(3,2)+n(2,1)+n(1,2))*log(rnew)+(n(4,2)+n(3,1)+n(2,2)+n(1,1))*log(1-rnew);
  r(4)=r(5)=r(6)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(3)=r(1)=abs(1.0-r(0));
  return(r);
}
Rcpp::NumericVector rf_B1_B1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2.0*n(3,3)+2.0*n(3,2)+n(3,1)+2.0*n(2,3)+2.0*n(2,2)+n(2,1)+2.0*n(1,3)+2.0*n(1,2)+n(1,1));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+n(1,1))*rold+n(3,2)+n(3,1)+n(2,3)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=((n(3,2)+n(2,3))*log(rnew-rnew*rnew)+(n(3,1)+n(2,1))*log(rnew)+(n(1,3)+n(1,2))*log(rnew/2.0)+(2.0*n(3,3)+2.0*n(2,2)+n(1,1))*log(1-rnew));
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+n(1,1))*rold+n(3,3)+n(3,1)+n(2,2)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=((n(3,3)+n(2,2))*log(rnew-rnew*rnew)+(n(3,1)+n(2,1))*log(rnew)+(n(1,3)+n(1,2))*log(rnew/2.0)+(2.0*n(3,2)+2.0*n(2,3)+n(1,1))*log(1-rnew));
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2*n(3,3)+n(2,1)+n(3,1)+2*n(2,2)+n(1,1))-2*M_LN2*(n(3,2)+n(2,3)+n(1,3)+n(1,2));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+2*n(1,1))*rold+n(3,2)+n(3,1)+n(2,3)+2.0*n(2,2)+n(1,3))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,2)+n(2,3))*log(rnew-rnew*rnew)+(n(3,1)+2*n(2,2))*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,2)*log(-(rnew-1.0)/2.0)+(2.0*n(3,3)+n(2,1))*log(1-rnew)-M_LN2*n(1,1);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+2*n(1,1))*rold+n(3,3)+2*n(3,2)+n(2,2)+n(2,1)+n(1,3))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(3,3)+n(2,2))*log(rnew-rnew*rnew)+(n(2,1)+2*n(3,2))*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,2)*log(-(rnew-1.0)/2.0)+(2.0*n(2,3)+n(3,1))*log(1-rnew)-M_LN2*n(1,1);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2.0*n(3,3)+n(3,2)+2.0*n(3,1)+2.0*n(2,3)+n(2,2)+2.0*n(2,1)+2.0*n(1,3)+n(1,2)+2.0*n(1,1));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2.0*n(1,3)+4*n(1,2)+2.0*n(1,1))*(rold*rold*rold) +
	    (2.0*n(3,2)+4*n(3,1)+2.0*n(2,3)+2.0*n(2,2)+2.0*n(2,1)-4*n(1,2)-2.0*n(1,1))*(rold*rold) +
	    (-2.0*n(3,2)-4*n(3,1)-2.0*n(2,3)-2.0*n(2,1)-n(1,3)+2.0*n(1,2)+n(1,1))*rold + 
	    n(3,2)+2.0*n(3,1)+n(2,3)+n(2,1)+n(1,3))/((2.0*rold*rold-2.0*rold+1)*2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=n(2,2)*log(2.0*rnew*rnew-2.0*rnew+1)+(n(2,3)+n(2,1))*log(rnew-rnew*rnew)+n(3,2)*log(2.0*rnew-2.0*rnew*rnew)+2.0*n(3,1)*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,1)*log(-(rnew-1)/2.0)+2.0*n(3,3)*log(1-rnew)-M_LN2*n(1,2);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     


  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2.0*n(1,3)+4*n(1,2)+2.0*n(1,1))*(rold*rold*rold) +
	    (2.0*n(3,3)+2.0*n(3,2)+2.0*n(3,1)+2.0*n(2,2)+4*n(2,1)-4*n(1,2)-2.0*n(1,1))*(rold*rold) +
	    (-2*n(3,3)-2*n(3,1)-2*n(2,2)-4*n(2,1)-n(1,3)+2*n(1,2)+n(1,1))*rold + 
	    +n(3,3)+n(3,1)+n(2,2)+2*n(2,1)+n(1,3))/((2.0*rold*rold-2.0*rold+1)*2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=n(3,2)*log(2.0*rnew*rnew-2.0*rnew+1)+(n(3,3)+n(3,1))*log(rnew-rnew*rnew)+n(2,2)*log(2.0*rnew-2.0*rnew*rnew)+2.0*n(2,1)*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,1)*log(-(rnew-1)/2.0)+2.0*n(2,3)*log(1-rnew)-M_LN2*n(1,2);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = LN_75*(n(3,1) + n(2,1) + n(1,1)) - 2.0*M_LN2*(n(3,2)+n(2,2)+n(1,2));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (n(2,2)+n(2,1)-2*n(1,2)-4*n(1,1))*rold*rold*rold+
	    (-2*n(3,1)-3*n(2,2)-n(2,1)+4*n(1,1))*rold*rold+
	    (2*n(3,1)+3*n(2,2)-2*n(2,1)+n(1,2)-3*n(1,1))*rold-2*n(3,1)-2*n(2,2)-2*n(1,2))/
	(((rold-2)*(rold*rold-rold+1))*2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=n(2,1)*log(rnew*rnew-rnew+1)+
    n(3,1)*log(2*rnew-rnew*rnew)+
    n(2,2)*log(rnew-rnew*rnew)+
    n(1,2)*log(rnew/2)+
    n(1,1)*log(-(rnew-2)/2)+
    2*n(3,2)*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	       
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (n(3,2)+n(3,1)-2*n(1,2)-4*n(1,1))*rold*rold*rold+
	    (-3*n(3,2)-n(3,1)-2*n(2,1)+4*n(1,1))*rold*rold+
	    (3*n(3,2)-2*n(3,1)+2*n(2,1)+n(1,2)-3*n(1,1))*rold-2*n(3,2)-2*n(2,1)-2*n(1,2))/
	(((rold-2)*(rold*rold-rold+1))*2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=n(3,1)*log(rnew*rnew-rnew+1)+
    n(2,1)*log(2*rnew-rnew*rnew)+
    n(3,2)*log(rnew-rnew*rnew)+
    n(1,2)*log(rnew/2)+
    n(1,1)*log(-(rnew-2)/2)+
    2*n(2,2)*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}

/*
rf_B1_D1
rf_B1_D2
rf_B2_B2
rf_B2_B3
rf_B2_C
rf_B2_D1
rf_B2_D2
rf_B3_B3
rf_B3_C
rf_B3_D1
rf_B3_D2
rf_C_C
rf_C_D1
rf_C_D2
rf_D1_D1
rf_D1_D2
rf_D2_D2
*/



