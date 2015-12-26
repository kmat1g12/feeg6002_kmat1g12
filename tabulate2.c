#include <stdio.h>
#include <math.h>

#define N 100
#define XMIN 1
#define XMAX 10

double x;
double y;
double z;
int i;

int main(void){
	for (i=0; i<N; i++){
	x = XMIN + (XMAX - XMIN) / (N - 1.) * i;
	y = sin(x);
	z = cos(x);
	printf("%f %f %f\n", x,y,z);
	}

	return 0;
}
