
#ifndef __CONFIG_H

#define __CONFIG_H

#define SINGLE_PRECISION		

#ifdef SINGLE_PRECISION		// to save memory & poss. increase performance
typedef float crd_type;
typedef float proj_type;
#else
typedef double crd_type;
typedef double proj_type;
#endif

#endif 		// __CONFIG_H
