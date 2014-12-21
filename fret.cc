/* 
 * given a pdb file, a set of trajectories, and a list of "native"
 * contacts and cutoffs, compute the fraction of native contacts
 * over the trajectories.
 */

//#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <unistd.h>

#include "traj.h"

string usage = "\n\n	Usage:\n"
"	distance -D donor1 -d donor2 -A acceptor1 -a acceptor2 -R R0...\n"
"		-C donor_centre -c acceptor_centre go1.dcd ... goN.dcd\n"
"               calculate kappa2 and kET over trajectory\n\n";

double vec_norm(double &x, double &y, double &z)
{
	double mag;
	mag = sqrt(x*x+y*y+z*z);
	x/=mag; y/=mag; z/=mag;
	return mag;
}

double vec_dotp(double &x1,double &y1, double &z1, 
		double &x2,double &y2, double &z2)
{
	double dot;
	dot = x1*x2+y1*y2+z1*z2;
	return dot;
}

int main( int argc, char ** argv )
{
	vector<string> dcd_names;
	int c,natom;
	double com[3], mat[3][3], dr,dx,dy,dz;
	int donor1, donor2, acceptor1, acceptor2, centre1, centre2;
	double donor_life, R0, prefac;

	while (1) {
		c=getopt(argc,argv,"hD:d:A:a:R:C:c:t:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage.c_str());
				exit(0);
				break;
			case 'D':
				donor1 = atoi(optarg);
				break;
			case 'd':
				donor2 = atoi(optarg);
				break;
			case 'A':
				acceptor1 = atoi(optarg);
				break;
			case 'a':
				acceptor2 = atoi(optarg);
				break;
			case 'R':
				R0 = atof(optarg);
				break;
			case 'C':
				centre1 = atoi(optarg);
				break;
			case 'c':
				centre2 = atoi(optarg);
				break;
			case 't':
				donor_life = atof(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage.c_str());
				exit(1);
		}
	}

	if (argc - optind < 1) {
		fprintf(stdout,"%s\n",usage.c_str());
		exit(0);
	}
	dcd_names.resize(argc-optind);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[optind+t]);
	}

	//DCDITrajFile inptraj;
	BaseITrajFile *inptraj;
	int slen = dcd_names[0].size();
	if (slen-dcd_names[0].rfind(string(".dcd"))==4) {
		fprintf(stderr,"dcd file\n");
		inptraj = new DCDITrajFile(dcd_names[0].c_str());
	} else if (slen-dcd_names[0].rfind(string(".xtc"))==4) {
		fprintf(stderr,"xtc file\n");
		inptraj = new XTCITrajFile(dcd_names[0].c_str());
	} else {
		fprintf(stderr,"unknown trajectory file type:\n");
		fprintf(stderr,"%s\n",dcd_names[0].c_str());
		exit(1);
	}
	//inptraj.open(dcd_names[0].c_str());
	natom = inptraj->num_atoms();
	inptraj->close();
	//
	float *X = new float[natom];
	float *Y = new float[natom];
	float *Z = new float[natom];
	int frames = 0;
	donor1 -= 1; 
	donor2 -= 1; 
	acceptor1 -= 1; 
	acceptor2 -= 1; 
	centre1 -= 1; 
	centre2 -= 1; 
	prefac = 1.5 * pow(R0,6.0) / donor_life;
	// first pass to compute mean coordinates
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		int slen = dcd_names[trajfile].size();
		if (slen-dcd_names[trajfile].rfind(string(".dcd"))==4) {
			fprintf(stderr,"dcd file\n");
			inptraj = new DCDITrajFile(dcd_names[trajfile].c_str());
		} else if (slen-dcd_names[trajfile].rfind(string(".xtc"))==4) {
			fprintf(stderr,"xtc file\n");
			inptraj = new XTCITrajFile(dcd_names[trajfile].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",dcd_names[trajfile].c_str());
			exit(1);
		}
		while (inptraj->frames_left()) {
			double kET, kappa2, R, r_dot_muD, r_dot_muA, muD_dot_muA;
			double dx_D, dy_D, dz_D;
			double dx_A, dy_A, dz_A;
			double dx_R, dy_R, dz_R;
			//
			frames++;
			inptraj->read_frame(X,Y,Z,natom);
			dx_D = X[donor1]-X[donor2];
			dy_D = Y[donor1]-Y[donor2];
			dz_D = Z[donor1]-Z[donor2];
			dx_A = X[acceptor1]-X[acceptor2];
			dy_A = Y[acceptor1]-Y[acceptor2];
			dz_A = Z[acceptor1]-Z[acceptor2];
			dx_R = X[centre1]-X[centre2];
			dy_R = Y[centre1]-Y[centre2];
			dz_R = Z[centre1]-Z[centre2];
			vec_norm(dx_D,dy_D,dz_D);
			vec_norm(dx_A,dy_A,dz_A);
			R = vec_norm(dx_R,dy_R,dz_R);
			r_dot_muD = vec_dotp(dx_R,dy_R,dz_R,dx_D,dy_D,dz_D);
			r_dot_muA = vec_dotp(dx_R,dy_R,dz_R,dx_A,dy_A,dz_A);
			muD_dot_muA = vec_dotp(dx_D,dy_D,dz_D,dx_A,dy_A,dz_A);
			kappa2 = muD_dot_muA - 3.*r_dot_muD*r_dot_muA;
			kappa2 = kappa2*kappa2;
			kET = prefac * kappa2 / pow(R,6.0);
			fprintf(stdout,"%i %12.6e %12.6e %12.6e\n",frames, R, kappa2, kET);
		}
		inptraj->close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

