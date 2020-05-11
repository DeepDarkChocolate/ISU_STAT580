#include <stdio.h> 
#include <stdlib.h> 
#include "array.h"

#define MIN(i,j) ((i)<(j) ? (i) : (j))
#define MAX(i,j) ((i)>(j) ? (i) : (j))

void dsyevd_( char* jobz, char* uplo, int* n, double* a, int* lda,
	                double* w, double* work, int* lwork, int* iwork, int* liwork, int* info );

int evd(double **a, int n, double *lambda, double **p){
  int lda = n, info, lwork, liwork;
  int iwkopt;
  int* iwork;
  double wkopt;
  double* work;
  int i, j, k;
  double *AT;

  MAKE_VECTOR(AT, n*n);
  for (j=0, k=0; j<n; j++)
	for (i=0; i<n; i++)
	  AT[k++] = a[i][j];
  
  /* Query and allocate the optimal workspace */
  lwork = -1;
  liwork = -1;
  dsyevd_( "Vectors", "Upper", &n, AT, &lda, lambda, &wkopt, &lwork, &iwkopt,
	                          &liwork, &info );
  lwork = (int)wkopt;								          
  liwork = iwkopt;												         

  MAKE_VECTOR(work, lwork);
  MAKE_VECTOR(iwork, liwork);

  dsyevd_( "Vectors", "Upper", &n, AT, &lda, lambda, work, &lwork, iwork,
	                          &liwork, &info );
  
  for (j=0, k=0; j<n; j++){
	for (i=0; i<n; i++){
	  p[i][j] = AT[k++];
	}
  }
  FREE_VECTOR(AT);
  FREE_VECTOR(work);
  FREE_VECTOR(iwork);
  
  return info;
}
