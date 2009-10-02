#include "alea.h"
#include <stdio.h>
#include <sys/types.h>
#include <time.h>

double alea( double a, double b)
{
	//random number (uniform distribution) in [a b]
	double r ;
	r = (double)rand();
	r = r / RAND_MAX;
	return a + r * (b-a);
}


int alea_integer(int a, int b)
{
	// Integer random number in [a b]
	int ir;
	double r;
	r = alea(0, 1);
	ir= (int)(a + r*(b+1-a));
	if (ir>b) ir = b;
        return ir;
}
                                  



