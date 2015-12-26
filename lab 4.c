#include <stdio.h>
#include <limits.h>  /* limits for integers */
#include <float.h>   /* limits for floats */
#include <math.h>

double upper_bound(long n);

long maxlong(void) {
       return LONG_MAX;
    }

long factorial(long x){

  long i;
  long val = 1;
  long fin;

  if(x<0){
    fin = -2;
  }

  else{

    for (i = 1; i <= x; i++){
      val = val*i; 
    }

    if (val<maxlong()){
      fin = val;
    }
    if (val>=maxlong()){
      fin = -1;
    }
  }

  return fin;
}

double upper_bound(long n){
    long nmin = 0;
    long nsix = 6;
    long ans=0;

    if(ans==0){
      if(n>= nmin && n< nsix){
          ans = factorial(n);
      }
      if(n>=nsix){
          ans =  pow((n/2.),n);
      }
      if(n<0){
          ans = factorial(n);
      }
    }

    if(ans<maxlong() && ans>=-2)
      return ans;

    if(ans>=maxlong())
      return factorial(n);

    return -1;
}

int main(void){

	long i;
	printf("maxlong()=%ld\n", maxlong());
	
	for (i=-5; i<15; i++) {
	printf("upper_bound(%ld)=%g\n", i, upper_bound(i));
	}

	return 0;
}