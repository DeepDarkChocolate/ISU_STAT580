#include<stdio.h>
#include<stdlib.h>
#define MATHLIB_STANDALONE  /*It is essential to have this before the call 
							to the Rmath’s header file because this determines
							the definitions to be set.*/
#include <Rmath.h>

/* gcc -Wall -pedantic -o Rmathex Rmathex.c -lRmath  */

double dgammac(double x, double shape, double rate){
  return dgamma(x, shape, rate, 1);
}

int main(void) {

  double x[] = {0.5, 1, 2};
  double shape[] = {1, 2, 3};
  double rate[] = {0.5, 1, 2};

  for (int i = 0; i < 3; i++){
	for(int j = 0; j < 3; j++){
	  for(int k = 0; k < 3; k++){
		printf("log(f(%g | shape = %g, rate = %g)) = %g\n", x[k], shape[j], rate[i], dgammac(x[k], shape[j], rate[i]));
		/* Return loglikelihood of gamma density:logf(x) with shape para = shape, rate para = rate*/
	  }
	}
  }

  return EXIT_SUCCESS;
}
