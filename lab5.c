#include <stdio.h>
#include <string.h>
#define MAXLINE 1000 /* maximum length of string */

/* function prototype */
void reverse(char source[], char target[]);

long string_length(char s[]){
	long length;
	length = strlen(s);
	return length;
}

void reverse(char source[], char target[]){
  long slen = string_length(source);
  long x = 0;
  long i;

  for(i = slen-1; i>=0; i--){
    target[x] = source[i];
    x = x+1;
  }
}

int main(void) {
  char original[] = "This is a test: can you print me in reverse character order?";
  char reversed[MAXLINE];

  printf("%s\n", original);
  reverse(original, reversed);
  printf("%s\n", reversed);
  return 0;
}
