/* @file pass_value
 * @author Yonghyun Kwon
 * */

#include <stdio.h>
#include <stdlib.h>

void ftn(int i, int *pi){
  printf("----Inside the function before increment----\n");
  printf("address of the integer is &i = %p\n", (void *)&i);
  printf("value of the pointer is pi = %p\n", (void *)pi);
  printf("address of the pointer to the integer is &pi = %p\n", (void *)&pi);
  printf("---------------------------------------------\n");
  i = i + 1;
  pi = pi + 1;
  printf("----Inside the function after increment----\n");
  printf("address of the integer is &i = %p\n", (void *)&i);
  printf("value of the pointer is pi = %p\n", (void *)pi);
  printf("address of the pointer to the integer is &pi = %p\n", (void *)&pi);
  printf("--------------------------------------------\n");
}

int main(void){
  int d = 0;
  int *pd = &d;
  printf("---Outside the function---\n");
  printf("address of the integer is &d = %p\n", (void *)&d);
  printf("value of the pointer is pd = %p\n", (void *)pd);
  printf("address of the pointer to the integer is &pd = %p\n", (void *)&pd);
  printf("--------------------------\n");
  ftn(d, pd);
  printf("---Outside the function---\n");
  printf("address of the integer is &d = %p\n", (void *)&d);
  printf("value of the pointer is pd = %p\n", (void *)pd);
  printf("address of the pointer to the integer is &pd = %p\n", (void *)&pd);
  printf("--------------------------\n");
  return 0;
}
