#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{

double multiplier;
int window_size;

window_size = atoi(argv[1]);

for (int i = 0; i <= window_size; i++) 
{
     	multiplier = 0.5 * (1 - cos(2*M_PI*i/window_size));
	printf("%f\n",multiplier);
}

printf("\n");
return 0;

}
