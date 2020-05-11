#include<stdlib.h>
#include<stdio.h>
#include<string.h>

int main(void){
  
  char c[4] = "ACGT";
  int lenc = strlen(c);
  int i1; int i2; int i3; int i4; int i5; int i6; int i7; int i8;
  for(i1 = 0; i1 < lenc; i1++){
	for(i2 = 0; i2 < lenc; i2++){
	  for(i3 = 0; i3 < lenc; i3++){
		for(i4 = 0; i4 < lenc; i4++){
		  for(i5 = 0; i5  < lenc; i5++){
			for(i6 = 0; i6  < lenc; i6++){
			  for(i7 = 0; i7  < lenc; i7++){
				for(i8 = 0; i8  < lenc; i8++){
				 /* printf("%c%c%c%c%c%c%c%c\n", c[i1], c[i2], c[i3], c[i4], c[i5], c[i6], c[i7], c[i8]);*/
				}
			  }
			}
		  }
		}
	  }
	}
  }
  int cnt = 0;
  int strlen = 8; /* change to 16 if you want*/
  int max = 65536;  /* change to 2^strlen if you want */
  while(cnt < max){
	int tmp = cnt;
	for(int i = 0; i < strlen; i++){
	  printf("%c", c[tmp & 3]);
	  tmp = tmp >>2;
	}
	printf("\n");
	cnt++;
  }
}
