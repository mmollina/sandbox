/* Written by Marcelo Mollinari
   Adapted from hmm_main.c, hmm_f2.c and util.c (found in the R package qtl)
   copyright (c) 2001-10, Karl W Broman    
   These codes are under the GNU General Public License, version 3
   A copy of the GPL 3 is available at http://www.r-project.org/Licenses/GPL-3
*/


#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define THRESH 200.0

double emitf(int obs_gen, int true_gen, double error_prob)
{
  switch(obs_gen) {
  case 0: return(1.0);
  case 1:
    switch(true_gen) {
    case 1: return(1.0-error_prob);
    case 2: case 3: case 4: return(error_prob/2.0);
    }
  case 2:
    switch(true_gen) {
    case 1: case 4: return(error_prob/2.0);
    case 2: case 3: return(1.0-error_prob);
    }
  case 3:
    switch(true_gen) {
    case 4: return(1.0-error_prob);
    case 1: case 2: case 3: return(error_prob/2.0);
    }
  case 4: /* AA or AB (not BB) */
    if(true_gen != 4) return(1.0-0.5*error_prob);
    else return(error_prob);
  case 5: /* AB or BB (not AA) */
    if(true_gen != 1) return(1.0-0.5*error_prob);
    else return(error_prob);
  }
  return(1.0); /* shouldn't get here */
}


double stepf(int gen1, int gen2, double rf)
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return((1.0-rf)*(1.0-rf));
    case 2: case 3: return((1.0-rf)*rf);
    case 4: return(rf*rf);
    }
  case 2:
    switch(gen2) {
    case 1: case 4: return((1.0-rf)*rf);
    case 2: return((1.0-rf)*(1.0-rf));
    case 3: return(rf*rf);
    }
  case 3:
    switch(gen2) {
    case 1: case 4: return((1.0-rf)*rf);
    case 3: return((1.0-rf)*(1.0-rf));
    case 2: return(rf*rf);
    }
  case 4:
    switch(gen2) {
    case 1: return(rf*rf);
    case 2: case 3: return((1.0-rf)*rf);
    case 4: return((1.0-rf)*(1.0-rf));
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double nrecf(int gen1, int gen2, double rf)
{
  switch(gen1) {
  case 1:
    switch(gen2) {
    case 1: return(0.0);
    case 2: case 3: return(0.5);
    case 4: return(1.0);
    }
  case 2:
    switch(gen2) {
    case 1: case 4: return(0.5);
    case 2: return(0.0);
    case 3: return(1.0);
    }
  case 3:
    switch(gen2) {
    case 1: case 4: return(0.5);
    case 3: return(0.0);
    case 2: return(1.0);
    }
  case 4:
    switch(gen2) {
    case 4: return(0.0);
    case 2: case 3: return(0.5);
    case 1: return(1.0);
    }
  }
  return(log(-1.0)); /* shouldn't get here */
}

double addlog(double a, double b)
{
  if(b > a + THRESH) return(b);
  else if(a > b + THRESH) return(a);
  else return(a + log1p(exp(b-a)));
}

// [[Rcpp::export]]
SEXP est_hmm_f2_brand_new(NumericMatrix Geno, NumericVector rf, int verbose) {
  int n_mar = Geno.nrow();
  int n_ind = Geno.ncol();
  int n_gen = 4;
  int it, i, v, v2, j, j2, flag=0, maxit=1000;
  double error_prob = 0.001;
  double tol=1e-6, s=0.0; 
  double loglik, curloglik; 
  NumericMatrix alpha(n_gen, n_mar);
  NumericMatrix beta(n_gen, n_mar);
  NumericMatrix gamma(n_gen, n_gen);
  NumericVector cur_rf(n_mar-1);
  NumericVector initf(4,0.25);

  NumericMatrix tr(n_gen, (n_mar-1)*n_gen);
  NumericMatrix em(5, (n_mar-1)*n_gen);

  for(it=0; it<maxit; it++) {
    //Rcpp::Rcout << "it: " << it << "\n";
    for(j=0; j<n_mar-1; j++) {
      cur_rf[j] = rf[j];
      rf[j] = 0.0;
    }
    /*
      3. se for testar ordens parecidas nao precisa alocar tudo outra vez. 
      Isso é complicado de programar mas salva um bom tempo de processamento
    */
    for(j=0; j < ((n_mar-1)*n_gen); j++){
      for(i=0; i<n_gen; i++){
	tr(i,j)= stepf(i+1, (j%4)+1, cur_rf(j/4));
      }
    }
    for(j=0; j < n_gen; j++){
      for(i=0; i < 5; i++){
	em(i,j)= emitf(i+1, j+1, error_prob);
      }
    }
    std::fill(rf.begin(), rf.end(), 0);
    for(i=0; i<n_ind; i++) { /* i = individual */
      R_CheckUserInterrupt(); /* check for ^C */
      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	alpha(v,0) = initf(v)  * em(Geno(0,i)-1, v);
	beta(v,n_mar-1) = 1.0;
      }
     
      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {

	for(v=0; v<n_gen; v++) {
	  alpha(v,j) = alpha(0,j-1) * tr(0, (j-1)*n_gen+v);
	  beta(v,j2) = beta(0,j2+1) * tr(v, j2*4) * 
	    em(Geno(j2+1,i)-1,0);

	  for(v2=1; v2<n_gen; v2++) {
	    alpha(v,j) = alpha(v,j) + alpha(v2,j-1) * tr(v2,(j-1)*n_gen+v);
	    beta(v,j2) = beta(v,j2) + beta(v2,j2+1) * tr(v, j2*n_gen+v2)  * em(Geno(j2+1,i)-1,v2);
	  }

	  alpha(v,j) *= em(Geno(j,i)-1,v);
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma(v,v2) = alpha(v,j) * beta(v2,j+1) *
	      em(Geno(j+1,i)-1, v2)*
	      tr(v, j*n_gen+v2);

	    if(v==0 && v2==0) s = gamma(v,v2);
	    else s = s + gamma(v,v2);
	  }
	}

	for(v=0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    rf(j) += nrecf(v+1,v2+1, cur_rf(j)) * gamma(v,v2)/s;
	  }
	}
      }

    } /* loop over individuals */
    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf(j) /= (double)n_ind;
      if(rf(j) < tol/1000.0) rf(j) = tol/1000.0;
      else if(rf(j) > 0.5-tol/1000.0) rf(j) = 0.5-tol/1000.0;
    }
    /* check convergence */    
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(rf(j) - cur_rf(j)) > tol*(cur_rf(j)+tol*100.0)) {
	flag = 1;
	break;
      }
    }
    if(!flag) break;

  } /* end EM algorithm */

  //if(flag) warning("Didn't converge!\n");

  // calculate log likelihood
  loglik = 0.0;
  
  for(i=0; i<n_ind; i++) { // i = individual 
    // initialize alpha 
    for(v=0; v<n_gen; v++) {
      alpha(v,0) = initf(v) * emitf(Geno(0,i), v+1, error_prob);
    }
    // forward equations 
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha(v,j) = alpha(0,j-1) *
	  stepf(1, v+1, rf(j-1));

	for(v2=1; v2<n_gen; v2++)
	  alpha(v,j) = alpha(v,j) + alpha(v2,j-1) *
	    stepf(v2+1,v+1,rf(j-1));
	alpha(v,j) *= emitf(Geno(j,i),v+1,error_prob);
      }
    }
    curloglik = alpha(0,n_mar-1);
    for(v=1; v<n_gen; v++)
      curloglik = curloglik + alpha(v,n_mar-1);
    loglik += log(curloglik);
  }
  
  List z = List::create(wrap(rf), wrap(loglik));
  return(z);
}
