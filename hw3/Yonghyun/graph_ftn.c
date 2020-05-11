#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "array.h"

void getfield(char* line, int num, double *ptr){
  const char* tok;
  char *ptr0;
  int i = 0;
  for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",\n")){
	if(i == 0) {i++; continue;}
	ptr[i - 1] = strtod(tok, &ptr0);
	if(i == num) break;
	i++;
  }
}

typedef double(*sMeasure)(double*, double*, int p);
typedef void(*prod)(double *, double *, double *, int *, int *, int lenc);
int n_k(int i, int j, sMeasure s, int k, double ** mat, int nrow, int ncol);
double corr(double *x, double *y, int n);
double corrCoef(double *x, double *y, int n);
void ccsprod_sym (double *x, double *y, double a[], int r[], int c[], int lena);
void print_dvector(double *a, int rows, const char *format);
void print_dmatrix(double **a, int rows, int cols, const char *format);
void eigen(int m, double *evals, double **evecs, prod ccs_ptr, double *a, int *r, int *c, int lena, int lenc);


void graph(int maxrow,  int k, int m, double rho){
  FILE* stream = fopen("diurnaldata.csv", "r");
  char line[1024];
  int row = 0;
  int col = 11;
  int i; int j;
  double * tmpcol = malloc(sizeof(double)* col);
  int n_knum;

  double * mat[maxrow];
  for(i = 0; i < maxrow; i++) mat[i] = (double *) malloc(sizeof(double) * col);
  fgets(line, 1024, stream);
  while (fgets(line, 1024, stream)){
	char * tmp = strdup(line);
	getfield(tmp, col, tmpcol);
	memcpy(mat[row], tmpcol, sizeof(double) * col);
	free(tmp);	
	row++;
	if(row ==  maxrow) break;
  }
  free(tmpcol);
  printf("FILE READ COMPLETED\n");

  int lena = maxrow * 3;
  int lenc = maxrow + 1;
  double *a = malloc(sizeof(double) * lena);
  int *r = malloc(sizeof(double) * lena);
  int *c = calloc(sizeof(int),  lenc);
  double *g = calloc(sizeof(double),  maxrow);
  double *one = malloc(sizeof(double) * maxrow);
  for(i = 0; i < maxrow; i++){
	one[i] = 1;
  }
  int cnt = 0;
  sMeasure sptr = &corrCoef;
  for(j= 0; j < maxrow; j++){
	for(i = j; i < maxrow; i++){
	  if(i == j){
		a[cnt] = 1; r[cnt] = i; cnt += 1;continue;
	  }
	  if((*sptr)(mat[i], mat[j], col) <= rho){
	    continue;
	  }
	  n_knum = n_k(i, j, sptr, k, mat, maxrow, col);
	  if(n_knum == 1){
		a[cnt] = 1;
		r[cnt] = i;
		cnt += 1;
	  }
	}
	c[j+1] = cnt; 
  }
  lena = cnt;

  printf("W calculation completed\n");

  ccsprod_sym(one, g, a, r, c, lena);
  cnt = 0;
  for(i = 0; i < lena; i++){
	a[i] /= sqrt(g[r[i]]);
	if(c[cnt + 1] <= i) cnt++;
	a[i] /= sqrt(g[cnt]);
  }
  double *evals; MAKE_VECTOR(evals, m);
  double **evecs; MAKE_MATRIX(evecs, lenc, m);

  prod pptr = &ccsprod_sym;
  
  /*print_dvector(a, lena, "%f "); printf("\n");
  for(i = 0; i < lenc; i++){printf("%d ", c[i]);} printf("\n");*/

  eigen(m, evals, evecs, pptr, a, r, c, lena, lenc - 1);
/*
  printf("eigen-vector:\n");
  print_dmatrix(evecs, lenc - 1, m, "%f  ");
  printf("eigen-value:\n");
  print_dvector(evals, m, "%f  ");
  printf("\n");
*/
  FILE *f = fopen("evals.txt", "wb");
  for(i = 0; i < m; i++){
	fprintf(f, "%f", evals[i]);
	if(i + 1 != m) fprintf(f, ", ");
  }
  fclose(f);

  FILE *f2 = fopen("evecs.txt", "wb");
  for(i = 0; i < maxrow; i++){
	for(j = 0; j < m; j++){
	  fprintf(f2, "%f", evecs[i][j]);
	  if(j + 1 != m) fprintf(f2, ", ");
	}fprintf(f2, "\n");
  }
  fclose(f2);


  free(a);
  free(r);
  free(c);
  free(one);
  free(g);
  FREE_VECTOR(evals);
  FREE_MATRIX(evecs);

  /*
  for(i = 0; i < maxrow; i++){
	for(int j = 0; j < col; j++){
	  printf("%f ", mat[i][j]);
	}
	printf("\n");
  }
  */
}
