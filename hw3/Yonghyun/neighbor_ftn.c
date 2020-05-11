#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"array.h"

typedef double(*sMeasure)(double*, double*, int p);

int n_k(int i, int j, sMeasure s, int k, double ** mat, int nrow, int ncol){
  double *pi = mat[i];
  double *pj = mat[j];
  double sij = (*s)(pi, pj, ncol);
  int cnt = 0;
  double *pj2;
  for(int j2 = 0; j2 < nrow; j2++){
	if((j2 == i) | (j2 == j)) continue;
	pj2 = mat[j2];
	if( (*s)(pi, pj2, ncol) > sij ) cnt += 1;
	if(cnt >= k) return (0);
  }
  
    
  int cnt2 = 0;
  double *pi2;
  for(int i2 = 0; i2 < nrow; i2++){
	if((i2 == i)| (i2 == j)) continue;
	pi2 = mat[i2];
	if( (*s)(pi2, pj, ncol) > sij) cnt2 += 1;
  	if(cnt2 >= k) return(0);
  }
  /*printf("cnt2 = %d\n", cnt2);*/ /* ... */
  return(1);
}

int n_k2(int i, int j, sMeasure s, int k, double ** mat, int nrow, int ncol){
  double *pi = mat[i];
  double *pj = mat[j];
  double sij = (*s)(pi, pj, ncol);
  int cnt = 0;
  double *pj2;
  for(int j2 = 0; j2 < nrow; j2++){
	if((j2 == i) | (j2 == j)) continue;
	pj2 = mat[j2];
	if( (*s)(pi, pj2, ncol) > sij ) cnt += 1;
  }
  if(cnt >= k) return (0);
    
  return(1);
}

void iso_point(sMeasure s, int k, double **mat, int nrow, int ncol){
  int cnt;
  for(int j = 0; j < nrow; j++){
	cnt = 0;
	for(int i = 0; i < nrow; i++){
	  cnt += n_k2(i, j, s, k, mat, nrow, ncol);
	}
	if(cnt == 1) printf("%d ",j);
  }
  printf("\n");
}


double sum(double *x, int n){
  double res = 0;
  for(int i = 0; i < n; i++){
	res += x[i];
  }
  return(res);
}

double mean(double *x, int n){
  double res = sum(x, n) / n;
  return(res);
}

double cov(double *x, double *y, int n){
  double res = 0;
  double meanx = mean(x, n);
  double meany = mean(y, n);
  for(int i = 0; i < n; i++){
	res += x[i] * y[i];
  }
  res -= n*meanx*meany;
  res /= (n-1);
  return(res);
}

double corr(double *x, double *y, int n){
  double res = cov(x, y, n) / sqrt(cov(x, x, n) * cov(y, y, n));
  return(res);
}

double corrCoef(double *X, double *Y, int n) 
{ 
  
    double sum_X = 0, sum_Y = 0, sum_XY = 0; 
    double squareSum_X = 0, squareSum_Y = 0; 
  
    for (int i = 0; i < n; i++) 
    { 
        // sum of elements of array X. 
        sum_X = sum_X + X[i]; 
  
        // sum of elements of array Y. 
        sum_Y = sum_Y + Y[i]; 
  
        // sum of X[i] * Y[i]. 
        sum_XY = sum_XY + X[i] * Y[i]; 
  
        // sum of square of array elements. 
        squareSum_X = squareSum_X + X[i] * X[i]; 
        squareSum_Y = squareSum_Y + Y[i] * Y[i]; 
    } 
  
    // use formula for calculating correlation coefficient. 
    double corr = (n * sum_XY - sum_X * sum_Y)  
                  / sqrt((n * squareSum_X - sum_X * sum_X)  
                      * (n * squareSum_Y - sum_Y * sum_Y)); 
  
    return corr; 
} 
