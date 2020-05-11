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
int multiply(double **a, int arows, int acols,
	        double **b, int brows, int bcols, double **c);
double dEnorm(double *x, int n);
int cpy(double **a,int nrows,int ncols,double **b);
void print_dmatrix(double **a, int rows, int cols, const char *format);
void print_dvector(double *a, int rows, const char *format);
void matrpose(double **a,int rows,int cols,double **aT);

typedef void(*prod)(double *, double *, double *, int *, int *, int lenc);

void ccsprod (double *x, double *y, double a[], int r[], int c[], 
		int lena); /* Size of an array should be passed in a function */ 
void ccsprod_sym (double *x, double *y, double a[], int r[], int c[], 
		int lena); /* Size of an array should be passed in a function */ 

void eigen(int m, double *evals, double **evecs, prod ccs_ptr, double *a, int *r, int *c, int lena, int lenc);

int main(void){
  int m = 5; int i;
  prod pptr = &ccsprod;
  
  int dim = 10;
  int lena = 2 * dim - 1;
  int lenc = dim + 1;
  double *a = malloc(sizeof(double) * lena);
  a[0] = 1;
  for(i = 0; i < dim - 1; i++){
	a[2*i + 1] = 1;
	a[2*i + 2] = i + 2;
  }
  int *r = malloc(sizeof(double) * lena);
  for(i = 0; i < dim - 1; i++){
	r[2*i + 2] = i + 1;
  }
  int * c = malloc(sizeof(int) * lenc);
  for(i = 0; i < dim; i++){
	c[i + 1] = 2*i + 1;
  }
  lenc -= 1;

  double *evals; MAKE_VECTOR(evals, m);
  double **evecs; MAKE_MATRIX(evecs, lenc, m);

  eigen(m, evals, evecs, pptr, a, r, c, lena, lenc);
  printf("eigen-vector:\n");
  print_dmatrix(evecs, lenc, m, "%f  ");
  printf("eigen-value:\n");
  print_dvector(evals, m, "%f  ");
  printf("\n");

  double ** mat; MAKE_MATRIX(mat, m, m);
  double ** mat2; MAKE_MATRIX(mat2, lenc, lenc);
  double ** tevecs; MAKE_MATRIX(tevecs, m, lenc);

  matrpose(evecs, lenc, m, tevecs);
  
  multiply(tevecs, m, lenc,
		      evecs, lenc, m, mat);
  printf("P^TP = \n");
  print_dmatrix(mat, m, m, "%f  ");

  FREE_VECTOR(evals);
  FREE_MATRIX(evecs);
}






