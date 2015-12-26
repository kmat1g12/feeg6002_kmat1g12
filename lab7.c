#include <stdio.h>
#include <stdlib.h>
#include <string.h>

long* make_fib_array(long n){

    /* Given an integer number n, return a list with
    the first n fibbonnaci numbers. Assume that n>=2""" 

    # create list fibs with n elements */
    long i;
    long *fibs;
    fibs = (long*)malloc(sizeof(long)*n);

    /* populate with data*/

    fibs[0] = 0;
    fibs[1] = 1;

    for(i=2;i<n;i++){

        fibs[i] = fibs[i - 1] + fibs[i - 2];
    }

    return fibs;

}
    
long* make_long_array(long n){

    long *a;
    a= (long*)malloc(sizeof(long)*n);
    
    return a;
    
    }
    

void use_make_long_array(long n) {
  int i;
  long *p = make_long_array(n);

  printf("In use_make_long_array(%ld)\n", n);

  /* if p is not NULL, we could allocate memory, and we proceed
     with testing: */
     
  if (p != NULL) {

    /* write some data to the array -- if the allocated memory
       is too short, this might trigger a segfault */
       
       for (i=0; i<n; i++) {
       p[i]=i+42;           /* just write some data */
       printf("%d ",p[i]);
       }
       printf("\n");

    /* free array -- if the allocated array is too short, we may
       have corrupted malloc/free metadata when writing the i+42 data
       above, and this may show when we call the free command: */

    free(p);

    /* if the program does not crash, it is a good sign [but
       no proof for correctness]. The other way round: if the
       program segfaults or crashes, this is not a good sign. */
  }
  
  else { /* we get here if memory allocation didn't work for
            some reason. */
            
        printf("Error - it seems that the memory allocation failed.\n");
  }
  
}


void use_fib_array(long N) {
  /* N is the maximum number for fibarray length */
  long n;      /* counter for fibarray length */
  long i;      /* counter for printing all elements of fibarray */
  long *fibarray;  /* pointer to long -- pointer to the fibarray itself*/

  /* Print one line for each fibarray length n*/
  for (n=2; n<=N; n++) {
    /* Obtain an array of longs with data */
    fibarray = make_fib_array(n);

    /* Print all elements in array */
    printf("fib(%2ld) : [",n);
    for (i=0; i<n; i++) {
      printf(" %ld", fibarray[i]);
    }
    printf(" ]\n");

    /* free array memory */
    free(fibarray);
  }
}


int main(void) {
  use_fib_array(30);
  return 0;
}
