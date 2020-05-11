#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "array.h"
void graph(int maxrow,  int k, int m, double rho);

int main(void){
  int maxrow = 22810; /* 22810*/
  int k = 3;
  int m = 20;
  double rho = 0.99;
  graph(maxrow, k, m, rho);
}

