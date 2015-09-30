#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.000001

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
  double l, l01, l02, rnew, rold;
  l01=-M_LN2*(2*n(3,4)+2*n(1,1))-2.0*M_LN2*(n(3,3)+n(3,2)+n(2,4)+n(2,1)+n(1,3)+n(1,2))-M_LN2*(2*n(3,1)+2*n(1,4))-2.0*M_LN2*(n(2,3)+n(2,2));
  l02=-2.0*M_LN2*(n(3,4)+n(3,1)+n(2,3)+n(2,2)+n(1,4)+n(1,1))-M_LN2*(2*n(3,3)+2*n(1,2))-M_LN2*(2*n(3,2)+2*n(1,3))-2.0*M_LN2*(n(2,4)+n(2,1));
  /*EM algorithm*/
  rold=0, rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,3)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    (n(2,2)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
	    2.0*(n(3,1)+n(1,4)) + n(3,3)+n(3,2)+n(2,4)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(2,3)+n(2,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,3)+n(3,2)+n(2,4)+n(2,1)+n(1,3)+n(1,2))*log(rnew-(rnew*rnew))+
    2.0*(n(3,1)+n(1,4))*log(rnew)+
    2.0*(n(3,4)+n(1,1))*log(1.0-rnew);
  r(4)=(l-l01)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,4)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    (n(2,1)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
	    2.0*(n(3,2)+n(1,3)) + n(3,4)+n(3,1)+n(2,3)+n(2,2)+n(1,4)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(2,4)+n(2,1))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,4)+n(3,1)+n(2,3)+n(2,2)+n(1,4)+n(1,1))*log(rnew-(rnew*rnew))+
    2.0*(n(3,2)+n(1,3))*log(rnew)+
    2.0*(n(3,3)+n(1,2))*log(1.0-rnew);
  r(5)=(l-l02)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,4)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    (n(2,1)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
	    2.0*(n(3,3)+n(1,2)) + n(3,4)+n(3,1)+n(2,3)+n(2,2)+n(1,4)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(2)=rnew;
  l=(n(2,4)+n(2,1))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,4)+n(3,1)+n(2,3)+n(2,2)+n(1,4)+n(1,1))*log(rnew-(rnew*rnew))+
    2.0*(n(3,3)+n(1,2))*log(rnew)+
    2.0*(n(3,2)+n(1,3))*log(1.0-rnew);
  r(6)=(l-l02)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,3)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    (n(2,2)*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
	    2.0*(n(3,4)+n(1,1)) + n(3,3)+n(3,2)+n(2,4)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(3)=rnew;
  l=(n(2,3)+n(2,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,3)+n(3,2)+n(2,4)+n(2,1)+n(1,3)+n(1,2))*log(rnew-(rnew*rnew))+
    2.0*(n(3,4)+n(1,1))*log(rnew)+
    2.0*(n(3,1)+n(1,4))*log(1.0-rnew);
  r(7)=(l-l01)/log(10.0); /*transforming to base 10 logarithm*/    
  return(r);
}





