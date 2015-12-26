#include <stdio.h>

int main(void){

	int s = 1000;
	float debt = s;
	float rate = 0.03;
	int N = 25;
	int month;
	float interest;
	float interest_sum = 0;
	float frac;

	for (month = 1; month < N; month++){
		interest = debt * rate;
		debt = debt + interest;
		interest_sum = interest_sum + interest;
		frac = interest_sum/s *100;
		printf("month %2d: debt = %.2f, interest= %.2f, total_interest= %.2f, frac= %.2f%% \n", month, debt, interest, interest_sum, frac);

		}
	return 0;
}
