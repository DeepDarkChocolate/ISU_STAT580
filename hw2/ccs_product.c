/* @ file name = ccs_product.c
 * @ author = Yonghyun Kwon
 * 
 * How to complie: gcc -o ccs_product ccs_product.c -ansi -Wall -pedantic
 * No warnings!
 * 
 * How to execute: ./ccs_product
 *
 */

#include<stdio.h>
#include<stdlib.h>

void cssprod (double x[], double y[], double a[], int r[], int c[], 
		size_t lena){ /* Size of an array should be passed in a function */
  int j, k;
  for(j = 0, k = 0; k < lena; k++){
    while(c[j+1] <= k){
      j++;
    }
    y[r[k]] += a[k] * x[j];
  }
}


int main(void){
  double x[] = {1, 2, 3, 4, 5.5, 6, 7};
  double a[] = {4, 2, 1, 13, 10, -4, 3, 4, 5};
  int r[] = {1, 3, 0, 2, 1, 2, 3, 0, 1};
  int c[] = {0, 2, 4, 6, 6, 7, 7, 9};
  double y[6] = {0};
  int i;
  size_t lena = sizeof(a) / sizeof(a[0]); /* Get size of each element*/
  size_t lenx = sizeof(x) / sizeof(x[0]);
  size_t leny = sizeof(y) / sizeof(y[0]);
  if(c[0] != 0){
    for(i = 0; i < lenx + 1; i++){
      c[i] =c[i] -c[0];
    }
  } /* Standardize array c so that the first element of c becomes 0*/
  cssprod(x, y, a, r, c, lena);
  for(i = 0; i < leny; i++){
    printf("y[%d] = %f\n", i, y[i]);
  }
  return 0;
}


