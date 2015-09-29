#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.00001


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
  int ct=0, mis=0, n_ind=250;
  double q, rold, rnew, r0, r1, r2;
  NumericMatrix r(((int)x.size()/n_ind), ((int)x.size()/n_ind));
  int n[5][5] = {};
  for(int i=0; i < (int)(x.size()/n_ind)-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < (int)x.size()/n_ind; j++)
        {
          rold=0, rnew=0.01;
          for (int row=0; row<5; row++)
	    {
	      for(int columns=0; columns<5; columns++)
		n[row][columns]=0;
	    }
	  mis=0;
          for(int k=0; k < n_ind; k++)
            {
              if(x(i*n_ind+k)==1) 
                {
		  if (x(j*n_ind+k)==1) n[0][0]++;
		  else if (x(j*n_ind+k)==2) n[0][1]++;
		  else if (x(j*n_ind+k)==3) n[0][2]++;
		  else if (x(j*n_ind+k)==4) n[0][3]++;
		  else if (x(j*n_ind+k)==5) n[0][4]++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==2) 
                {
		  if (x(j*n_ind+k)==1) n[1][0]++;
		  else if (x(j*n_ind+k)==2) n[1][1]++;
		  else if (x(j*n_ind+k)==3) n[1][2]++;
		  else if (x(j*n_ind+k)==4) n[1][3]++;
		  else if (x(j*n_ind+k)==5) n[1][4]++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==3) 
                {
		  if (x(j*n_ind+k)==1) n[2][0]++;
		  else if (x(j*n_ind+k)==2) n[2][1]++;
		  else if (x(j*n_ind+k)==3) n[2][2]++;
		  else if (x(j*n_ind+k)==4) n[2][3]++;
		  else if (x(j*n_ind+k)==5) n[2][4]++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==4) 
                {
		  if (x(j*n_ind+k)==1) n[3][0]++;
		  else if (x(j*n_ind+k)==2) n[3][1]++;
		  else if (x(j*n_ind+k)==3) n[3][2]++;
		  else if (x(j*n_ind+k)==4) n[3][3]++;
		  else if (x(j*n_ind+k)==5) n[3][4]++;
		  else mis++;
                }
              else if(x(i*n_ind+k)==5) 
                {
		  if (x(j*n_ind+k)==1) n[4][0]++;
		  else if (x(j*n_ind+k)==2) n[4][1]++;
		  else if (x(j*n_ind+k)==3) n[4][2]++;
		  else if (x(j*n_ind+k)==4) n[4][3]++;
		  else if (x(j*n_ind+k)==5) n[4][4]++;
		  else mis++;
                }
	      else mis++;
	    }
	  
	  /*No dominant markers*/
    	  if((n[0][3] + n[0][4] + n[1][3] + n[1][4] + n[2][3] + n[2][4] + 
	      n[3][0] + n[3][1] + n[3][2] + n[3][3] + n[3][4] + 
	      n[4][0] + n[4][1] + n[4][2] + n[4][3] + n[4][4]) == 0) 
	    {
	      while(abs(rold-rnew) > TOL)
	        {
		  rold=rnew;
		  q=(rold*rold)/((rold*rold)+((1-rold)*(1-rold)));
		  rnew=(2*(q*n[1][1]+n[2][0]+n[0][2])+n[2][1]+n[1][2]+n[1][0]+n[0][1])/(2*(n_ind-mis));
	        }
	      r(j,i)=r(i,j)=rnew;
	    }
	  /*One dominant marker 5-1*/
	  //else if()
	  /*One dominant marker 4-3*/
	  //else if()
	  /*Two dominant marker 5-1*/
	  //else if()
	  /*Two dominant marker 4-3*/
	    //  else if()
	  /*Two dominant marker 5-1 and 4-3*/
	  //else if()
	  else
	    {
	      while(abs(rold-rnew) > TOL)
	        {
		  rold=rnew;
		  r0=(1-rold)*(1-rold);
		  r1=(1-rold)*rold;
		  r2=rold*rold;
		  rnew=(n[0][1] + n[1][0] + n[1][2] + n[2][1] + //OK 1
			(n[0][4] + n[2][3] + n[3][2] + n[4][0])*((2*r1)/((r2+2*r1))) + //OK 2
			(n[1][3] + n[1][4] + n[3][1] + n[4][1])*(r1/(r1+r0+r2)) + //OK 3
			(n[0][3] + n[2][4] + n[3][0] + n[4][2])*((2*r1)/((r0+2*r1))) + //OK 4
			(n[3][3] + n[4][4])*(4*r1)/(3*r0 + 4*r1 + 2*r2)+ //OK 5
			(n[3][4] + n[4][3])*(4*r1)/(2*r0 + 4*r1 + 3*r2)+ //OK 6
			2*(n[0][2]+n[2][0]+ //OK 8
			   (n[0][4] + n[2][3] + n[3][2] + n[4][0])*(r2/((r2+2*r1))) +  //OK 2
			   n[1][1]*(r2/(r0+r2)) + //OK 7 
			   (n[1][3] + n[1][4] + n[3][1] + n[4][1]) *(r2/(r1+r0+r2)) +  //OK 3
			   (n[3][3] + n[4][4])*(2*r2)/(3*r0 + 4*r1 + 2*r2)+ // OK 5
			   (n[3][4] + n[4][3])*(3*r2)/(2*r0 + 4*r1 + 3*r2)   //OK 6
			   ))/(2*(n_ind-mis));
		}
	      r(j,i)=r(i,j)=rnew;
	     }

	  /*Mixed markers 1-2-3, 5-1 and 4-3*/
	  //else if()
	  /*Should not get here*/
	      // else 
	      //return(-1);
	  /*

	      for (int row=0; row<5; row++)
	      {
	        for(int columns=0; columns<5; columns++)
	          printf("%d \t", n[row][columns]);
	        printf("\n");
	      }
	      printf("\n----------------------------------------------------\n");
	  */
	}
    }
  return wrap(r);
}

