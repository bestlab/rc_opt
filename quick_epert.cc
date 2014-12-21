/*
 * QUICK ENERGY PERTURBATION
 *
 * input: dcd format trajectory
 * 	  file with a list in the following format:
 * 	  	i	j	delta_eij	rij
 * 	  	....
 * 
 * output: a single column list of energies
 *
 * usage:
 * 	quick_epert parm.dat traj1.dcd ... trajN.dcd > epert.dat
 */

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "TrajFile.h"
#include "random_gen.h"
#include <sys/time.h>

const char *Usage = "\n\nquick_epert parm.dat traj.dcd > epert.dat\n\n";

using namespace std;

int main(int argc, char **argv)
{
	const int maxij = 100;
	DCDITrajFile inptraj;
	vector<string> dcdfiles;
	int potential_type = 0;		// 0:KB; 1: ...
	double morse_alpha = 1.2;
	double fsw, riju, rijl, rul3;
	int c,natom, nframe, ii,jj, argind;
	int I[maxij],J[maxij];
	double A[maxij],B[maxij],C[maxij];
	double *X,*Y,*Z,rij,deij;
	int npair = 0;
	double dx,dy,dz,dr2,dr4;
	double E, EE, tmpr, r_on2, r_off2;
	double r_on = 25.0;
	double r_off = 27.5;

	while (1) {
		c=getopt(argc,argv,"hktlm:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'k':
				potential_type = 0;	// Karanicolas
				break;			// 12-10-6
			case 't':
				potential_type = 1;	// 12-10
				break;
			case 'l':
				potential_type = 2;	// 12-6
				break;
			case 'm':
				potential_type = 3;	// morse
				morse_alpha = atof(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	argind = optind;
	fstream ijinp(argv[argind++]);
	inptraj.open(argv[argind]);
	natom = inptraj.num_atoms();
	inptraj.close();
	X = new double[natom];
	Y = new double[natom];
	Z = new double[natom];
	npair = 0;
	while (ijinp.good()) {
		ijinp >> I[npair] >> J[npair] >> deij >> rij;
		A[npair] = 13.0*deij*pow(rij,double(12));
		B[npair] = -18.0*deij*pow(rij,double(10));
		C[npair] = 4.0*deij*pow(rij,double(6));
		npair++;
	}
	npair--;
	//printf("Number of dcd files= %i\n",argc-argind);
	//exit(0);
	dcdfiles.resize(argc-argind);
	for (int i=0; i< argc-argind; i++) {
		dcdfiles[i] = argv[argind+i];
	}
	r_off2 = r_off*r_off;
	r_on2 = r_on*r_on;
	tmpr = r_off2-r_on2;
	rul3 = 1.0/(tmpr*tmpr*tmpr);
	for (int file = 0; file<dcdfiles.size(); file++) {
		inptraj.open(dcdfiles[file].c_str());
		nframe = inptraj.total_frames();
		for (int f=0; f<nframe; f++) {
			inptraj.read_frame(X,Y,Z,natom); 
			E = 0.0;
			EE = 0.0;
			for (int p=0; p<npair; p++) {
				ii = I[p]-1;
				jj = J[p]-1;
				dx = X[ii]-X[jj];
				dy = Y[ii]-Y[jj];
				dz = Z[ii]-Z[jj];
				dr2 = dx*dx+dy*dy+dz*dz;
				dr4 = dr2*dr2;
				tmpr = A[p]/(dr4*dr4*dr4) + B[p]/(dr4*dr4*dr2)
					+ C[p]/(dr4*dr2);
				E += tmpr;
				if (sqrt(dr2)>r_off) {
					fsw = 0.0;
				} else if (sqrt(dr2)<r_on) {
					fsw = 1.0;
				} else {
					rijl = r_on2-dr2;
					riju = r_off2-dr2;
					fsw = riju*riju*(riju-3.0*rijl)*rul3;
				}
				EE += tmpr*fsw;
			}
			cout << EE << endl;
		}
		inptraj.close();
	}
	delete [] X;
	delete [] Y;
	delete [] Z;

	return 0;
}
