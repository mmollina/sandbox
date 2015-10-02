#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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

arma::mat Tr(double r)
{
  arma::mat tr(4,4);
  tr(0,0)=tr(1,1)=tr(2,2)=tr(3,3)=(1-r)*(1-r);
  tr(0,1)=tr(0,2)=tr(1,0)=tr(1,3)=tr(2,0)=tr(2,3)=tr(3,1)=tr(3,2)=(1-r)*r;
  tr(0,3)=tr(1,2)=tr(2,1)=tr(3,0)=r*r;
  return(tr);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat est_rf_arma(NumericVector x)
{
  int n_ind=250, n_mar=((int)x.size()/n_ind);
  arma::mat T(4,4);
  arma::mat n(5,5);
  arma::mat I3(3,1,fill::ones);
  arma::mat I2(2,1,fill::ones);
  arma::mat U(4,4,fill::ones);
  arma::mat D(4,4,fill::ones);
  D(0,0)=D(1,1)=D(2,2)=D(3,3)=0;
  D(0,3)=D(1,2)=D(2,1)=D(3,0)=2;
  arma::mat IB1(3,4,fill::zeros);
  IB1(0,0)=IB1(0,1)=IB1(1,2)=IB1(2,3)=1;
  arma::mat IC(2,4,fill::zeros);
  IC(0,0)=IC(0,1)=IC(0,2)=IC(1,3)=1;
  arma::mat M(IB1.n_rows, IC.n_rows);
  arma::mat H(IB1.n_rows, IC.n_rows);
  for(int i=0; i < n_mar-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < n_mar; j++)
        {
	  n=count_genotypes(x,i,j,n_ind);
	  T=Tr(.01);
	  H=IB1*(((U*trans(IC))/(trans(IB1)*IB1*U*trans(IC)))%((T)*trans(IC)));
	  M=IB1*(((U*trans(IC))/(trans(IB1)*IB1*U*trans(IC)))%((T%D)*trans(IC)));
	  //Rcpp::Rcout << "\nlikelihood: " << I3%(log(H)*n)%I2 << "\n"; //ERRO: tenho que pegar uma submatrix 
	}
    }
  //List z = List::create(H, n, U, IB1, C);
  return(M) ;
}
