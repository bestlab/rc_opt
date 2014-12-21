
#ifndef _RANDOM_GEN_H
#define _RANDOM_GEN_H

using namespace std;

const int NTAB = 32;
const int IA = 16807;
const long IM =  2147483647;
const double AM =  1.0/2147483647.0;
const long IQ = 127773;
const int IR = 2836;
const long NDIV = ( 1+(2147483647-1)/32);
const double EPS = 1.2e-7;
const double RNMX = 1.0 - 1.2e-7;

class park_miller {
	/* from numerical recipes p280 */
	/* generates a random number from a uniform distribution on interval [0,1] */
	private:
		long idum;
		long iy;
		long iv[NTAB];
	public: 
		park_miller(long idum_ini=1);
		~park_miller();
		void initialize(long idum_ini);
		double random();
};

class box_muller {
	/* generates a random number from a gaussian distribution with zero mean and
	 * unit variance */
	private:
		park_miller pm;
		int iset;
		double gset;
	public:
		box_muller(long idum_ini=1);
		~box_muller();
		void initialize(long idum_ini);
		double random();
};


#endif	// #ifdef _RANDOM_GEN_H
