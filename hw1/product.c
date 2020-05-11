/* @file product.c
 * @author Yonghyun Kwon
 *
 * How to compile
 * gcc -o product product.c -ansi -Wall -pedantic
 *
 * How to execute
 * ./product
 *
 * Results
 * Type short: 16384
 * Type int: 81920
*/ 

# include <stdio.h>

int main(void){
  int i = 320; int j = 256;
  short res1 = i * j;
  int res2 = i * j;
  printf("Type short:%d\n", res1);
  printf("Type int:%d\n", res2);
  return 0;
}

/* We can find out that the type short output is 16384
 * and the type int output is 81920, which is expected.
 * Storing as short type, arithmetic overflow occurs,
 * since the maximum value of short type is 2^15 - 1 = 32767
 * Short type consists of 2 Bytes(=16 Bits), so 82920(mod 2^16) = 16384 is obtained.
 * Storing as int type, however, artihmetic overflow does not occur,
 * since the maximum value of int type is 2^31 - 1.
 * Int type consists of 4 Bytes(=32 Bits) so 82920(mod 2^32) = 82920 is obtained.
*/
