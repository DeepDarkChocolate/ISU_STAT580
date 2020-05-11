/*  Calculate the SVD of a A using LAPACK */
/* 
gcc -c mat_vec.c -Wall -pedantic
gcc -c eval.c -Wall -pedantic
gcc -o testeval testeval.c eval.o mat_vec.o -Wall -pedantic -lm -llapack
*/

#include<stdio.h> 
#include "mat_vec.h"
#include "array.h"

int evd(double **a, int n, double *lambda, double **p);
void print_dvector(double *a, int rows, const char *format);
void print_dmatrix(double **a, int rows, int cols, const char *format);

int main(void)
{
	double **a,**p, *lambda;
	int i,j,m,n,flag;

	m = 9;
	n = 9;
	if(m != n){
	  fprintf(stderr, "ERROR: not symmetric matrix\n");
	  exit(1);
	}

	MAKE_MATRIX(a,m,n);
	for(i = 0; i < n; i++){
	  for(j = 0; j < m; j++){
		if(i == j) a[i][j] = i+2;
		else if ((i == 0) | (j == 0)) a[i][j] = 1;
		else a[i][j] = 0;
	  }
	}
	MAKE_MATRIX(p,m,n);
	MAKE_VECTOR(lambda,n);

	flag = evd(a, n, lambda, p);
	if (flag!=0) {
	  fprintf(stderr, "The algorithm failed to compute eigenvalues.\n");
	  exit(1);
	}	
	else {
	  printf("Original matrix:\n");
	  print_dmatrix(a, m, n, "%g ");
	  printf("\nEigen Values:\n");
	  print_dvector(lambda, n, "%g "); printf("\n"); 
	  printf("\nEigen Vectors:\n");
	  print_dmatrix(p,m,m,"%g ");		
	  FREE_MATRIX(a);		
	  FREE_MATRIX(p);
	  FREE_VECTOR(lambda);
	}
	return EXIT_SUCCESS;
}	
  
