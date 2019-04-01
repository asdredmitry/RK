#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const double EPS = 1e-10;
double f(double t)
{
	return cos(t);
}
double xi(double x1, double x2)
{
	return (x1 - f(x1)*(x2 - x1)/(f(x2) - f(x1)));
}
double find_zero(double t1, double t2, double eps)
{
	double tmp = t2;
	while(fabs(t1 - t2) > eps)
	{
		tmp = xi(t1, t2);
		t2 = t1;
		t1 = tmp;
	}
	return tmp;
}
int main()
{
	printf("%lf \n", find_zero(4.5, 6, EPS));
	return 0;
}