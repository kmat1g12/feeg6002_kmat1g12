/* Laboratory 6, SESG6025, 2013/2014, Template */

#include <stdio.h>
#include <string.h>
#define MAXLINE 1000

/* Function void rstrip(char s[])
modifies the string s: if at the end of the string s there are one or more spaces,
then remove these from the string.

The name rstrip stands for Right STRIP, trying to indicate that spaces at the 'right'
end of the string should be removed.
*/

void reverse(char source[], char target[]){
  long slen = strlen(source);
  long x = 0;
  long i;

  for(i = slen-1; i>=0; i--){
    target[x] = source[i];
    x = x+1;
  }
}

void rstrip(char s[]) {
    long i;
    i = 0;

    while(s[i]!='\0'){
    	i++;
    }

    while(s[i-1] == ' '){
    	s[i-1] = '\0';
    	i--;
    }
   
}	

void lstrip(char s[]) {
    long len;

    char reversed[MAXLINE];

    reverse(s,reversed);
    rstrip(reversed);
    len = strlen(reversed);
    reverse(reversed,s);
    s[len] = '\0';
} 


int main(void) {
  char test1[] = "Hello World   ";
  char test2[] = "   Hello World";

  printf("Original string reads  : |%s|\n", test1);
  rstrip(test1);
  printf("r-stripped string reads: |%s|\n", test1);
  printf("Original string reads  : |%s|\n", test2);
  lstrip(test2);
  printf("l-stripped string reads: |%s|\n", test2);


  return 0;
}