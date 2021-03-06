#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 1e-6
#define LN_75 -0.28768207245178 

Rcpp::NumericVector est_rf_C_C(std::vector<int> k_sub,
			       std::vector<int> k1_sub,
			       int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==1) 
	{
	  if (k1_sub[k]==1) n1++;
	  else if (k1_sub[k]==2) n2++;
	  else if (k1_sub[k]==3) n3++;
	  else n0++;

	}
      else if(k_sub[k]==2) 
	{
	  if (k1_sub[k]==1) n2++;
	  else if (k1_sub[k]==2) n4++;
	  else if (k1_sub[k]==3) n2++;
	  else n0++;
	}
      else if(k_sub[k]==3) 
	{
	  if (k1_sub[k]==1) n3++;
	  else if (k1_sub[k]==2) n2++;
	  else if (k1_sub[k]==3) n1++;
	  else n0++;
	}
      else n0++;
    }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(n2+2*(n3+n4*rold*rold/((1-rold)*(1-rold)+rold*rold)))/(2*(n_ind-n0));
    }
  //Likelihood
  l=n4*log(rnew*rnew+(1-rnew)*(1-rnew))+
    2*n3*log(rnew)+n2*log(rnew*(1-rnew))+2*n1*log(1-rnew);
  //Likelihood unde H0: r=0.5
  l0=-M_LN2*(n4+2.0*n3+2.0*n2+2*n1);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}

Rcpp::NumericVector est_rf_C_D_43(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n3=0, n4=0, n5=0, n6=0, n8=0;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==1) 
	{
	  if (k1_sub[k]==3) n3++;
	  else if (k1_sub[k]==4) n4++;
	  else n0++;
	}
      else if(k_sub[k]==2) 
	{
	  if (k1_sub[k]==3) n6++;
	  else if (k1_sub[k]==4) n8++;
	  else n0++;
	}
      else if(k_sub[k]==3) 
	{
	  if (k1_sub[k]==3) n1++;
	  else if (k1_sub[k]==4) n5++;
	  else n0++;
	}
      else n0++;
    }
  //n=count_genotypes_f2( x, i, j, n_ind);
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(n6+ 
	    n5*2*(1-rold)*rold/(rold*rold+2*(1-rold)*rold) +
	    n8*(1-rold)*rold/((1-rold)*rold+(1-rold)*(1-rold)+rold*rold) +
	    n4*2*(1-rold)*rold/((1-rold)*(1-rold)+2*(1-rold)*rold) + 
	    2*(n3 + n5*rold*rold/(rold*rold+2*(1-rold)*rold) +
	       n8*rold*rold/((1-rold)*rold+(1-rold)*(1-rold)+rold*rold)))/(2*(n_ind-n0));
    }
  //Likelihood
  l=n1 * (2.0*log(1.0-rnew)) + 
    n6  *  (log(rnew*(1.0-rnew))) + 
    n5  *  (log(1.0-(1.0-rnew)*(1.0-rnew))) + 
    n8  *  (log(1.0-rnew*(1.0-rnew))) + 
    n4  *  (log(1.0-rnew*rnew)) + 
    n3 * (2.0*log(rnew));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n3 + n6 + n1) + LN_75*(n4 + n8 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_C_D_51(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n3=0, n4=0, n5=0, n6=0, n8=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==1) 
	{
	  if (k1_sub[k]==1) n1++;
	  else if (k1_sub[k]==5) n5++;
	  else n0++;
	}
      else if(k_sub[k]==2) 
	{
	  if (k1_sub[k]==1) n6++;
	  else if (k1_sub[k]==5) n8++;
	  else n0++;
	}
      else if(k_sub[k]==3) 
	{
	  if (k1_sub[k]==1) n3++;
	  else if (k1_sub[k]==5) n4++;
	  else n0++;
	}
    
      else n0++;
    }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    { 
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=(n6 + n5*2.0*r1/(r2+2*r1) +
	    n8*r1/(r1+r0+r2) +
	    n4*2*r1/(r0+2*r1) + 
	    2*(n3 + n5*r2/(r2+2*r1) +
	       n8*r2/(r1+r0+r2)))/(2.0*(n_ind-n0));
    }
  r0=(1.0-rnew)*(1.0-rnew);
  r1=rnew*(1.0-rnew);
  r2=rnew*rnew;
  //Likelihood
  l=n1 * (2.0*log(1.0-rnew)) + n6  *  (log(r1)) + 
    n5  *  (log(1.0-r0)) + n8  *  (log(1.0-r1)) + 
    n4  *  (log(1.0-r2)) + n3 * (2.0*log(rnew));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n3 + n6 + n1) + LN_75*(n4 + n8 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_D_D_43(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0, n6=0, n7=0, n8=0, n9=0, n10=0, n11=0, n12=0, n13=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==3) 
	{
	  if (k1_sub[k]==3) n1++;
	  else if (k1_sub[k]==4) n5++;
	  else n0++;
	}
      else if(k_sub[k]==4) 
	{
	  if (k1_sub[k]==3) n11++;
	  else if (k1_sub[k]==4) n12++;
	  else n0++;
	}
      else n0++;
    }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    { 
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((n5 + n11)*2.0*r1/(r2+2*r1) +
	    n12*4*r1/(3*r0 + 4*r1 + 2*r2) +
	    2*((n5 + n11)*r2/(r2+2*r1) +
	       n12*2*r2/(3*r0 + 4*r1 + 2*r2)))/(2.0*(n_ind-n0));
    }
  r0=(1.0-rnew)*(1.0-rnew);
  r1=rnew*(1.0-rnew);
  r2=rnew*rnew;
  //Likelihood
  l=n1 * (2.0*log(1.0-rnew)) + 
    n5  *  (log(1.0-r0)) + 
    n11  *  (log((1.0-r0) / 3.0)) + 
    n12 * (log((r0 + 2.0) / 3.0)); 

  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n11 + n1) + LN_75*(n12 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_D_D_51(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n5=0, n11=0, n12=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==1) 
	{
	  if (k1_sub[k]==1) n1++;
	  else if (k1_sub[k]==5) n5++;
	  else n0++;
	}
      else if(k_sub[k]==5) 
	{
	  if (k1_sub[k]==1) n11++;
	  else if (k1_sub[k]==5) n12++;
	  else n0++;
	}
      else n0++;
    }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    { 
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((n5 + n11)*2.0*r1/(r2+2*r1) +
	    n12*4*r1/(3*r0 + 4*r1 + 2*r2) +	 
	    2*((n5 + n11)*r2/(r2+2*r1) +
	       n12*2*r2/(3*r0 + 4*r1 + 2*r2)))/(2.0*(n_ind-n0));
    }

  //Likelihood
  l=n1 * (2.0*log(1.0-rnew)) + 
    n5  *  (log(1.0-(1.0-rnew)*(1.0-rnew))) + 
    n11  *  (log((1.0-(1.0-rnew)*(1.0-rnew)) / 3.0)) + 
    n12 * (log(((1.0-rnew)*(1.0-rnew) + 2.0) / 3.0));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n11 +  n1) + LN_75*(n12 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_D_D_43_51(std::vector<int> k_sub,
				     std::vector<int> k1_sub,
				     int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0, n6=0, n7=0, n8=0, n9=0, n10=0, n11=0, n12=0, n13=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==3) 
	{
	  if (k1_sub[k]==1) n3++;
	  else if (k1_sub[k]==5) n4++;
	  else n0++;
	}
      else if(k_sub[k]==4) 
	{
	  if (k1_sub[k]==1) n9++;
	  else if (k1_sub[k]==5) n13++;
	  else n0++;
	}
      else n0++;
    }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    { 
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((n4 + n9)*2*r1/(r0+2*r1) + 
	    n13*4*r1/(2*r0 + 4*r1 + 3*r2) +
	    2*(n3 + n13*3*r2/(2*r0 + 4*r1 + 3*r2)))/(2.0*(n_ind-n0));
    }
  r0=(1.0-rnew)*(1.0-rnew);
  r1=rnew*(1.0-rnew);
  r2=rnew*rnew;
  //Likelihood
  l=
    n4  *  (log(1.0-r2)) + 
    n9  *  (log((1.0-r2) / 3.0)) + 
    n13 * (log((r2+2.0) / 3.0)) + 
    n3 * (2.0*log(rnew));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n3 + n9) + LN_75*(n13 + n4);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_A_A(std::vector<int> k_sub,
			       std::vector<int> k1_sub,
			       int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0, n6=0, n7=0, n8=0, n9=0, n10=0, n11=0, n12=0, n13=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
    {
      if(k_sub[k]==1) 
	{
	  if (k1_sub[k]==1) n1++;
	  else if (k1_sub[k]==2) n2++;
	  else if (k1_sub[k]==3) n3++;
	  else if (k1_sub[k]==4) n4++;
	  else if (k1_sub[k]==5) n5++;
	  else n0++;
	}
      else if(k_sub[k]==2) 
	{
	  if (k1_sub[k]==1) n6++;
	  else if (k1_sub[k]==2) n7++;
	  else if (k1_sub[k]==3) n6++;
	  else if (k1_sub[k]==4) n8++;
	  else if (k1_sub[k]==5) n8++;
	  else n0++;
	}
      else if(k_sub[k]==3) 
	{
	  if (k1_sub[k]==1) n3++;
	  else if (k1_sub[k]==2) n2++;
	  else if (k1_sub[k]==3) n1++;
	  else if (k1_sub[k]==4) n5++;
	  else if (k1_sub[k]==5) n4++;
	  else n0++;
	}
      else if(k_sub[k]==4) 
	{
	  if (k1_sub[k]==1) n9++;
	  else if (k1_sub[k]==2) n10++;
	  else if (k1_sub[k]==3) n11++;
	  else if (k1_sub[k]==4) n12++;
	  else if (k1_sub[k]==5) n13++;
	  else n0++;
	}
      else if(k_sub[k]==5) 
	{
	  if (k1_sub[k]==1) n11++;
	  else if (k1_sub[k]==2) n10++;
	  else if (k1_sub[k]==3) n9++;
	  else if (k1_sub[k]==4) n13++;
	  else if (k1_sub[k]==5) n12++;
	  else n0++;
	}
      else n0++;
    }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
    { 
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((n2 + n6)+ 
	    (n5 + n11)*2.0*r1/(r2+2*r1) +
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
	       ))/(2.0*(n_ind-n0));
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
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}







/*

  if(i==10 && j==15)
  {
  for(int bla=0; bla < 5; bla++)
  { 
  for(int ble=0; ble < 5; ble++)
  { 
  Rcpp::Rcout << n(bla,ble) << " ";
  }
  Rcpp::Rcout << "\n";
  }
  }

  if(i==0 && j==1)
  {
  for(int ble=0; ble < n_ind; ble++)
  { 
  Rcpp::Rcout << k_sub[ble] << "--" << k1_sub[ble] << "\n";
  }

  }



*/
