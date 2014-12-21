
#include "random_gen.h"
#include <cmath>
#include <cstdio>


park_miller::park_miller(long idum_ini)
{
	initialize(idum_ini);
}

park_miller::~park_miller()
{
}

void park_miller::initialize(long idum_ini) 
{
	long k;
	if (idum_ini <= 0) 
		idum = 1;
	else
		idum = idum_ini;
	// load shuffle table
	for (int j=NTAB+7; j>=0; j--) {
		k = idum/IQ;
		idum = IA*(idum-k*IQ)-IR*k;
		if (idum < 0) 
			idum += IM;
		if (j<NTAB)
			iv[j] = idum;
	}
	iy = iv[0];
}
	
double park_miller::random()
{
	//printf("idum = %Li\n", idum);
	int j;
	long k;
	double temp;
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	//printf("idum = %Li\n", idum);
	if (idum<0)
		idum += IM;
	//printf("idum = %Li\n", idum);
	j = iy/NDIV;
	iy = iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX)
		return RNMX;
	else
		return temp;
}

box_muller::box_muller(long idum_ini)
{
	initialize(idum_ini);
}

box_muller::~box_muller()
{}

void box_muller::initialize(long idum_ini)
{
	pm.initialize(idum_ini);
	iset = 0;
}

double box_muller::random()
{
	double fac, rsq, v1, v2;
	if (iset == 0) {
		do {
			v1 = 2.0*pm.random() - 1.0;
			v2 = 2.0*pm.random() - 1.0;
			rsq = v1*v1+v2*v2;
		} while ( rsq >= 1.0 || rsq == 0.0 );
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	} else {
		iset = 0;
		return gset;
	}
}



