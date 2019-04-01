#include <stdio.h>
#include <stdlib.h>
void reall(double ** data)
{
	int i = 0;
	data[0] = (double *)malloc(4*sizeof(double));
	for(i = 0; i < 4; i++)
		data[0][i] = i;
	data[0] = realloc(data[0], sizeof(double)*10);
}
int main(void)
{
	int i;
	double * data;
	reall(&data);
	for(i = 0; i < 10; i++)
		printf("%lf ", data[i]);
	free(data);
	return 0;
}
		 
