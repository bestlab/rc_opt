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
"	distance -i idxi -j idxj go1.dcd ... goN.dcd\n"
"               Just measures the distance in Angstroms between"
"               atoms/beads with indices idxi and idxj over the trajectories\n\n";


int main( int argc, char ** argv )
{
	string pdbfile, qlist, output;
	vector<string> dcd_names;
	int c,nselect,resi,resj,natom;
	int * selection;
	vector<int> i, j;
	vector<double> rij;
	double com[3], mat[3][3], dr,dx,dy,dz;
	double beta,gamma;

	beta = 5.0;
	gamma = 1.0;

	while (1) {
		c=getopt(argc,argv,"hi:j:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage.c_str());
				exit(0);
				break;
			case 'i':
				resi = atoi(optarg);
				break;
			case 'j':
				resj = atoi(optarg);
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
	resi -= 1;	// index from 0
	resj -= 1;
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
		//inptraj->open(dcd_names[trajfile].c_str());
		// iterate over frames ...
		//inptraj.open( dcd_names[trajfile].c_str() );
		//int no_frames = inptraj.total_frames();
		//int no_frames = inptraj.actual_frames();
		//for (int frame = 0; frame < no_frames; frame++) {
		while (inptraj->frames_left()) {
			frames++;
			inptraj->read_frame(X,Y,Z,natom);
			dx = X[resi]-X[resj];
			dy = Y[resi]-Y[resj];
			dz = Z[resi]-Z[resj];
			dr = sqrt(dx*dx+dy*dy+dz*dz);
			//fprintf(stdout,"%8.3f %8.3f %8.3f %8.3f\n",Z[resi],Z[resj],dz,dr);
			fprintf(stdout,"%12.6f %12.6f %12.6f %12.6f\n",dr,dx,dy,dz);
		}
		inptraj->close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

