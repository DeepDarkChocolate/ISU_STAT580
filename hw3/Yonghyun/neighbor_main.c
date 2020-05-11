#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"array.h"

typedef double(*sMeasure)(double*, double*, int p);

int n_k(int i, int j, sMeasure s, int k, double ** mat, int nrow, int ncol);
int n_k2(int i, int j, sMeasure s, int k, double ** mat, int nrow, int ncol);
int n_k_wp(int i, int j, double *pi, sMeasure s, int k, double ** mat, int nrow, int ncol);
void iso_point(sMeasure s, int k, double **mat, int nrow, int ncol);

double corrCoef(double *x, double *y, int n);

int main(void){
  double **mat2;
  int nrow2 = 10;
  int ncol2 = 5;
  MAKE_MATRIX(mat2, nrow2, ncol2);
  int i; int j;
  for(i = 0; i < nrow2; i++){
    for(j = 0; j < ncol2; j++){
	  mat2[i][j] = 1/(fabs((double)i-(double)j) + 1);
    }
  }

  int k = 1;
  sMeasure sptr = &corrCoef;
  

  for(int i3 = 0; i3 < nrow2; i3++){
	for(int j3 = 0; j3 < nrow2; j3++){
	  printf("%d ", n_k(i3, j3, sptr, k, mat2, nrow2, ncol2));
	}printf("\n");
  }


  printf("isolated points:\n");
  iso_point(sptr, k, mat2, nrow2, ncol2);
  FREE_MATRIX(mat2);
  return(0);
}
