/* Provide eigenvalue and eigenvector using the power method.
 * 04/04/2020 
 * author: Yonghyun Kwon
 * Disclaimer: The parameters of epsilon and maximum iteration number are taken from powerFunction of matlib package in R.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"array.h"

int matxvec(double **a, int arows, int acols,double *x, int xrows, double *y);
double dEnorm(double *x, int n);
int cpy(double **a,int nrows,int ncols,double **b);
void print_dmatrix(double **a, int rows, int cols, const char *format);
void print_dvector(double *a, int rows, const char *format);

typedef void(*prod)(double *, double *, double *, int *, int *, int lenc);

void ccsprod (double *x, double *y, double a[], int r[], int c[], 
		int lena){ /* Size of an array should be passed in a function */ 
  int j, k;
  for(j = 0, k = 0; k < lena; k++){
    while(c[j+1] <= k){
      j++;
    }
	/*printf("y[%d] = %f\n", r[k], y[r[k]]);
	printf("a[%d]*x[%d] = %f\n",k, j,  a[k]*x[j]);*/
    y[r[k]] += a[k] * x[j];
  }
}
void ccsprod_sym (double *x, double *y, double a[], int r[], int c[], 
		int lena){ /* Size of an array should be passed in a function */ 
  int j, k;
  for(j = 0, k = 0; k < lena; k++){
    while(c[j+1] <= k){
      j++;
    }
    y[r[k]] += a[k] * x[j];
	if(j != r[k]) y[j] += a[k] * x[r[k]];
	/*printf("y[%d] = %f\n", r[k], y[r[k]]);
	printf("a[%d]*x[%d] = %f\n",k, j,  a[k]*x[j]);*/
  }
}



void eigen(int m, double *evals, double **evecs, prod ccs_ptr, double *a, int *r, int *c, int lena, int lenc){
  int i; int j; double eval; int mi; int steps; double lenv; double tmpa;
  double *v;
  double *v2;
  double *vtmp; MAKE_VECTOR(vtmp, lenc);
  double *v_old; MAKE_VECTOR(v_old, lenc);
  double *v_oldtmp; MAKE_VECTOR(v_oldtmp, lenc);
  int maxiter = 1000;
  for(j = 0; j < lenc; j++){
	v_oldtmp[j] = 1/sqrt(lenc);
  }
  for(mi = 0; mi < m; mi++){
	v = calloc(sizeof(double), lenc);
	memcpy(v_old, v_oldtmp, sizeof(double) * lenc);
	steps = 1;
	while(steps < maxiter){
	  memcpy(vtmp, v, sizeof(double) * lenc);
	  (*ccs_ptr)(v_old, v, a, r, c, lena);
	  for(i = 0; i < mi; i++){
		tmpa = 0;
		for(j = 0; j < lenc; j++){
		  tmpa += evecs[j][i] * vtmp[j];
		}
		/*print_dvector(v, lenc, "%f ");printf("\nblabla\n");*/
		for(j = 0; j < lenc; j++){
		  v[j] -= evals[i] * tmpa * evecs[j][i];
		}
	  }
	  lenv = dEnorm(v, lenc);
	  if(lenv == 0){
		fprintf(stderr, "lenv = 0 ERROR\n");
		exit(EXIT_FAILURE);
	  }
	  for(i = 0; i < lenc; i++){
		/*printf("evec[%d] = %f\n", i, v[i]);*/  /* to be removed */
		v[i] /= lenv;
	  }
	  memcpy(v_old, v, sizeof(double) * lenc);
	  steps += 1;
	}
	v2 = calloc(sizeof(double), lenc);
	(*ccs_ptr)(v, v2, a, r, c, lena);
	eval = 0;
	for(i = 0; i < lenc; i++){
	  eval += v[i] * v2[i];
	}
	for(i = 0; i < lenc; i++){
	  evecs[i][mi] = v[i];
	}
	evals[mi] = eval;
	free(v2);
	free(v);
  }
  FREE_VECTOR(vtmp);
  FREE_VECTOR(v_old);
  FREE_VECTOR(v_oldtmp);
}





