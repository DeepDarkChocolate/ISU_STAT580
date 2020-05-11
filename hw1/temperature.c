/* 
@file temperature.c
@author Yonghyun Kwon
  - converts temperature from Faherheit to the Celsius scale and vice-versa.
  - If the input value is 1, convert from F to C.
  - If the input value if 2, convert form C to F.

  -How to compile
  gcc -o temperature temperature.c -ansi -Wall -pedantic

  -How to excute
  ./temperature  

  -Results(Example)
  First)
  From Fahernheit to Celsius, enter 1
  From Celsius to Fahernheit, enter 2
  1
  Enter a value: 0
  0.00 F = -17.78 C

  Second)
  From Fahernheit to Celsius, enter 1
  From Celsius to Fahernheit, enter 2
  2
  Enter a value: 0
  0.00 C = 32.00 F

  Third)
  From Fahernheit to Celsius, enter 1
  From Celsius to Fahernheit, enter 2
  123
  Please enter suitable input value(1 or 2)
*/

#include <stdio.h>

int main(void){
  int i; /* i is the indicator variable. If i == 1, F -> C; If i ==2, C -> F */
  double i2, res; /* i2 is the input value; res is the output value */
  printf("From Fahernheit to Celsius, enter 1\nFrom Celsius to Fahernheit, enter 2\n");
  scanf("%d", &i);
  if((i != 1) & (i != 2)){
    printf("Please enter suitable input value(1 or 2)\n");
    return 0;
  }
  printf("Enter a value: ");
  scanf("%lf", &i2);
  if(i == 1){
    res = (i2 -32) *5/9;
    printf("%.2f F = %.2f C\n", i2, res);
  }else if (i == 2){
    res = (i2 * 9/5) + 32;
    printf("%.2f C = %.2f F\n", i2, res);
  }
  return 0;
}
