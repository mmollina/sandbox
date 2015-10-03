#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
#define TOL 1e-05
#define LN3 1.098612288668109
#define LN4 1.38629436111989
#define LN_75 -0.28768207245178

arma::mat Tr(double r, int ph)
{
  arma::mat tr(4,4);
  switch(ph){
  case 1:  // Couple-couple
    tr(0,0)=tr(1,1)=tr(2,2)=tr(3,3)=(1-r)*(1-r);
    tr(0,1)=tr(0,2)=tr(1,0)=tr(1,3)=tr(2,0)=tr(2,3)=tr(3,1)=tr(3,2)=(1-r)*r;
    tr(0,3)=tr(1,2)=tr(2,1)=tr(3,0)=r*r;   
    break;
  case 2: //Couple-Repulsion
    tr(0,1)=tr(1,0)=tr(2,3)=tr(3,2)=(1-r)*(1-r);
    tr(0,3)=tr(1,2)=tr(2,1)=tr(3,0)=tr(0,0)=tr(1,1)=tr(2,2)=tr(3,3)=(1-r)*r;
    tr(0,2)=tr(1,3)=tr(2,0)=tr(3,1)=r*r;
    break;
  case 3: //Repulsion-Couple
    tr(0,1)=tr(1,0)=tr(2,3)=tr(3,2)=r*r;
    tr(0,3)=tr(1,2)=tr(2,1)=tr(3,0)=tr(0,0)=tr(1,1)=tr(2,2)=tr(3,3)=(1-r)*r;
    tr(0,2)=tr(1,3)=tr(2,0)=tr(3,1)=(1-r)*(1-r);    
    break;
  case 4: //Repulsion-Repulsion
    tr(0,0)=tr(1,1)=tr(2,2)=tr(3,3)=r*r;   
    tr(0,1)=tr(0,2)=tr(1,0)=tr(1,3)=tr(2,0)=tr(2,3)=tr(3,1)=tr(3,2)=(1-r)*r;
    tr(0,3)=tr(1,2)=tr(2,1)=tr(3,0)=(1-r)*(1-r);
    break;
  }
  return(tr);
}

arma::mat Nr(int ph)
{
  arma::mat D(4,4,fill::ones);
  switch(ph){
  case 1:  // Couple-couple
    D(0,0)=D(1,1)=D(2,2)=D(3,3)=0;
    D(0,3)=D(1,2)=D(2,1)=D(3,0)=2;   
    break;
  case 2: //Couple-Repulsion
    D(0,1)=D(1,0)=D(2,3)=D(3,2)=0;
    D(0,2)=D(1,3)=D(2,0)=D(3,1)=2;
    break;
  case 3:  //Repulsion-Couple
    D(0,1)=D(1,0)=D(2,3)=D(3,2)=2;
    D(0,2)=D(1,3)=D(2,0)=D(3,1)=0;
    break;
  case 4: //Repulsion-Repulsion
    D(0,0)=D(1,1)=D(2,2)=D(3,3)=2;   
    D(0,3)=D(1,2)=D(2,1)=D(3,0)=0;
    break;
  }
  return(D);
}

arma::mat count_genotypes(Rcpp::NumericVector x,
			  int i,
			  int j,
			  int n_ind)
{
  int mis=0;
  arma::mat n(5,5);
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

arma::mat Mt(int type)
{
  if(type==2)
    {
      arma::mat I(3,4,fill::zeros);
      I(0,0)=I(0,1)=I(1,2)=I(2,3)=1;
      return(I);
    }
  else if(type==5)
    {
      arma::mat I(2,4,fill::zeros);
      I(0,0)=I(0,1)=I(0,2)=I(1,3)=1;
      return(I);
    }
}

Rcpp::NumericVector est_rf_pair(arma::mat n,
				int n_mar,
				int n_ind,
				int mis,
				int sg1,
				int sg2)
{
  double rold=0, rnew=0.01, loglike, loglike_ho=.1;
  NumericVector r(8);
  arma::mat I4(1,3,fill::ones);
  arma::mat I3(1,3,fill::ones);
  arma::mat I2(1,2,fill::ones);
  arma::mat U(4,4,fill::ones);
  arma::mat T(4,4);
  arma::mat D(4,4);
  arma::mat Ik=Mt(sg1);
  arma::mat Ik1=Mt(sg2);
  arma::mat M(Ik.n_rows, Ik1.n_rows);
  arma::mat H(Ik.n_rows, Ik1.n_rows);
  arma::mat A=((U*trans(Ik1))/(trans(Ik)*Ik*U*trans(Ik1)));
  for(int i=1; i <=2; i++)
    {
      loglike_ho=arma::as_scalar(I3*(log(Ik*(A%((Tr(0.5,i))*trans(Ik1))))%n(span(1,3), span(1,2)))*trans(I2));
      rold=0, rnew=0.01;	     
      while(abs(rnew-rold) > TOL) 
	{
	  rold=rnew;
	  T=Tr(rnew, i);
	  D=Nr(i);
	  H=Ik*(A%((T)*trans(Ik1)));
	  M=Ik*(A%((T%D)*trans(Ik1)));
	  loglike=arma::as_scalar(I3*(log(H)%n(span(1,3), span(1,2)))*trans(I2));
	  rnew=arma::as_scalar((I3*((M % n(span(1,3),span(1,2)))/H) * trans(I2)) / (2.0*(n_ind-mis)));
	} 
      r(2*(i-1))=rnew;
      r(2*(i-1)+1)=(loglike-loglike_ho)/log(10.0);  
    }
  r(4)=abs(1.0-r(2));
  r(6)=abs(1.0-r(0));
  r(5)=r(3);
  r(7)=r(1);
  return(r);
}
