#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.000001

// [[Rcpp::export]]
SEXP est_rf_out(NumericVector x, NumericVector segreg_type, NumericVector n) {
  int ct=0, mis=0, nrec=0, n_ind=n(0);
  int n11=0, n12=0, n13=0, n14=0, n21=0, n22=0, n23=0, n24=0, n31=0, n32=0, n33=0, n34=0, n41=0, n42=0, n43=0, n44=0;
  int n1, n2, n3, n4, temp;
  double q, rold, rnew, l1, l2, l3, l4, l0, l01, l02;
  NumericMatrix r1(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  NumericMatrix r2(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  NumericMatrix r3(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  NumericMatrix r4(((int)x.size()/n_ind), ((int)x.size()/n_ind));  
  for(int i=0; i < (int)(x.size()/n_ind)-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < (int)x.size()/n_ind; j++)
        {
	  n11=n12=n13=n14=n21=n22=n23=n24=n31=n32=n33=n34=n41=n42=n43=n44=0;
	  mis=0;
          for(int k=0; k < n_ind; k++)
            {
              if(x(i*n_ind+k)==1) 
                {
		  if (x(j*n_ind+k)==1) n11++;
		  else if (x(j*n_ind+k)==2) n12++;
		  else if (x(j*n_ind+k)==3) n13++;
		  else if (x(j*n_ind+k)==4) n14++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==2) 
                {
		  if (x(j*n_ind+k)==1) n21++;
		  else if (x(j*n_ind+k)==2) n22++;
		  else if (x(j*n_ind+k)==3) n23++;
		  else if (x(j*n_ind+k)==4) n24++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==3) 
                {
		  if (x(j*n_ind+k)==1) n31++;
		  else if (x(j*n_ind+k)==2) n32++;
		  else if (x(j*n_ind+k)==3) n33++;
		  else if (x(j*n_ind+k)==4) n34++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==4) 
                {
		  if (x(j*n_ind+k)==1) n41++;
		  else if (x(j*n_ind+k)==2) n42++;
		  else if (x(j*n_ind+k)==3) n43++;
		  else if (x(j*n_ind+k)==4) n44++;
		  else mis++;
                }
	      else mis++;
	    }
	  
	  if(segreg_type(i) == 1 && segreg_type(j) == 1)
	    { 

	      /*Four possible likage phases*/ 
	      l0=-2.0*M_LN2*(n_ind-mis);  
	      n1=n41+n32+n23+n14;   
	      n2=n44+n33+n22+n11;
	      n3=n43+n34+n21+n12;
	      n4=n42+n31+n24+n13;

	      r1(j,i)=(2.0*(n1)+n3+n4)/(2.0*(n_ind-mis));
	      l1=(n3+n4)*log((1-r1(j,i))*r1(j,i))+2.0*(n1)*log(r1(j,i))+2.0*(n2)*log(1-r1(j,i));
	      r1(i,j)=(l1-l0)/log(10.0); /*transforming to base 10 logarithm*/
	      
	      r2(j,i)=(2.0*(n4)+n2+n1)/(2.0*(n_ind-mis));
	      l2=(n2+n1)*log((1-r2(j,i))*r2(j,i))+2.0*(n4)*log(r2(j,i))+2.0*(n3)*log(1-r2(j,i));
	      r2(i,j)=(l2-l0)/log(10.0); /*transforming to base 10 logarithm*/
	      
	      r3(j,i)=(2.0*(n3)+n2+n1)/(2.0*(n_ind-mis));
	      l3=(n2+n1)*log((1-r3(j,i))*r3(j,i))+2.0*(n3)*log(r3(j,i))+2.0*(n4)*log(1-r3(j,i));
	      r3(i,j)=(l3-l0)/log(10.0); /*transforming to base 10 logarithm*/
	      
	      r4(j,i)=(2.0*(n2)+n3+n4)/(2.0*(n_ind-mis));
	      l4=(n3+n4)*log((1-r4(j,i))*r4(j,i))+2.0*(n2)*log(r4(j,i))+2.0*(n1)*log(1-r4(j,i));
	      r4(i,j)=(l4-l0)/log(10.0); /*transforming to base 10 logarithm*/
	    }


	  if((segreg_type(i) == 4 && segreg_type(j) == 1) ||(segreg_type(i) == 1 && segreg_type(j) == 4) )
	    { 
	      if(segreg_type(i) == 1 && segreg_type(j) == 4)
		{
		  temp=n12; n12=n21; n21=temp;
		  temp=n13; n13=n31; n31=temp;
		  temp=n14; n14=n41; n41=temp;
		  temp=n23; n23=n32; n32=temp;
		  temp=n24; n24=n42; n42=temp;
		  temp=n34; n34=n43; n43=temp;
		}
	      l01=-M_LN2*(2*n34+2*n11)-2.0*M_LN2*(n33+n32+n24+n21+n13+n12)-M_LN2*(2*n31+2*n14)-2.0*M_LN2*(n23+n22);
	      l02=-2.0*M_LN2*(n34+n31+n23+n22+n14+n11)-M_LN2*(2*n33+2*n12)-M_LN2*(2*n32+2*n13)-2.0*M_LN2*(n24+n21);
	      //Rcpp::Rcout << "-----------------------------------------" <<  "\n";
	      //Rcpp::Rcout << "i: " << i << "    j: " << j << "\n";
	      //Rcpp::Rcout << "\t" << n11 << "\t" << n12 << "\t" <<  n13 << "\t" << n14 <<  "\n";
	      //Rcpp::Rcout << "\t" << n21 << "\t" << n22 << "\t" << n23  << "\t" << n24 <<  "\n";
	      //Rcpp::Rcout << "\t" << n31 << "\t" << n32 << "\t" << n33 << "\t" << n34 <<  "\n";
	      //Rcpp::Rcout << "-----------------------------------------" <<  "\n";
	      /*EM algorithm*/
	      rold=0, rnew=0.01;
	      while(abs(rold-rnew) > TOL)
		{
		  rold=rnew;
		  rnew=((n23*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
			(n22*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
			2.0*(n31+n14) + n33+n32+n24+n21+n13+n12)/(2.0*(n_ind-mis));
		}
	      r1(j,i)=rnew;
	      l1=(n23+n22)*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
		(n33+n32+n24+n21+n13+n12)*log(rnew-(rnew*rnew))+
		2.0*(n31+n14)*log(rnew)+
		2.0*(n34+n11)*log(1.0-rnew);
	      r1(i,j)=(l1-l01)/log(10.0); /*transforming to base 10 logarithm*/
	      rold=0, rnew=0.01;	     
	      while(abs(rold-rnew) > TOL)
		{
		  rold=rnew;
		  rnew=((n24*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
			(n21*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
			2.0*(n32+n13) + n34+n31+n23+n22+n14+n11)/(2.0*(n_ind-mis));
		}
	      r2(j,i)=rnew;
	      l2=(n24+n21)*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
		(n34+n31+n23+n22+n14+n11)*log(rnew-(rnew*rnew))+
		2.0*(n32+n13)*log(rnew)+
		2.0*(n33+n12)*log(1.0-rnew);
	      r2(i,j)=(l2-l02)/log(10.0); /*transforming to base 10 logarithm*/
	      rold=0, rnew=0.01;	     
	      while(abs(rold-rnew) > TOL)
		{
		  rold=rnew;
		  rnew=((n24*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
			(n21*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
			2.0*(n33+n12) + n34+n31+n23+n22+n14+n11)/(2.0*(n_ind-mis));
		}
	      r3(j,i)=rnew;
	      l3=(n24+n21)*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
		(n34+n31+n23+n22+n14+n11)*log(rnew-(rnew*rnew))+
		2.0*(n33+n12)*log(rnew)+
		2.0*(n32+n13)*log(1.0-rnew);
	      r3(i,j)=(l3-l02)/log(10.0); /*transforming to base 10 logarithm*/
	      rold=0, rnew=0.01;	     
	      while(abs(rold-rnew) > TOL)
		{
		  rold=rnew;
		  rnew=((n23*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
			(n22*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) +
			2.0*(n34+n11) + n33+n32+n24+n21+n13+n12)/(2.0*(n_ind-mis));
		}
	      r4(j,i)=rnew;
	      l4=(n23+n22)*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
		(n33+n32+n24+n21+n13+n12)*log(rnew-(rnew*rnew))+
		2.0*(n34+n11)*log(rnew)+
		2.0*(n31+n14)*log(1.0-rnew);
	      r4(i,j)=(l4-l01)/log(10.0); /*transforming to base 10 logarithm*/	     
	    }
	}
    }
  List z = List::create(wrap(r1), wrap(r2), wrap(r3), wrap(r4));
  return(z);





}



 





