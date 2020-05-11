#include <Rinternals.h>

/*
How to compile:
R CMD SHLIB skewkt_Call.c
How to run a program:
Rscript run_skewkt_Call.R
*/

SEXP skewkt(SEXP e_in){
  double *e; 
  int n; int p;
  SEXP SUM;

  if(!isReal(e_in))
	error("[ERROR] First argument must be a double vector");

  SEXP Rdim = getAttrib(e_in, R_DimSymbol);
  n = INTEGER(Rdim)[0];
  p = INTEGER(Rdim)[1];

  e = REAL(e_in);
  
  double sum = 0;
  for(int i = 0; i < n; i++){
	for(int j = 0; j < n; j++){
	  double sum1 = 0;
	  for(int k = 0; k < p; k++){
		sum1 += e[i + k * n] * e[j+ k * n];
	  }
	  sum += sum1 * sum1 * sum1;
	}
  }
  sum /= (n * n);

  PROTECT(SUM = allocVector(REALSXP, 1));
  REAL(SUM)[0] = sum;
  UNPROTECT(1);
  return(SUM);
}


