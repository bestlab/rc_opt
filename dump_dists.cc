/*
 * dump contact map for each frame of trajectory
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <fstream>
#include "traj.h"

using namespace std;

const char *Usage = "\n\n\
	dump_dists -c cut -o outp.txt -L boxL [ -m mth ] [ -s skip ] trj.dcd\n\
	cut is a contact cutoff in A\n\
	output goes to outp.txt\n\
	option '-s' causes only every skip'th frame to be used\n\
	option '-m' causes only every m'th residue/atom to be considered\n\
	\n\n";

inline int index(int i, int j) 
{
	return j>=i ? -1: i*(i-1)/2+j;
}

void get_ij(int ind, int &i, int &j) 
{
	i = int(floor((1.0+sqrt(1.0+8.0*ind))/2.0));
	j = ind - i*(i-1)/2;
	return;
}


void dump_dists(vector<string> &dcd_names, double rcut, FILE *outp, int skip, int m, double boxL)
{
	BaseITrajFile *inptraj;
	int natom, nfile, dim, idx, istp;
	int nu;
	double dx, dy, dz, dr;
	nfile = dcd_names.size();
	int slen = dcd_names[0].size();
	if (slen-dcd_names[0].rfind(string(".dcd"))==4) {
		fprintf(stdout,"dcd file\n");
		inptraj = new DCDITrajFile(dcd_names[0].c_str());
	} else if (slen-dcd_names[0].rfind(string(".xtc"))==4) {
		fprintf(stdout,"xtc file\n");
		inptraj = new XTCITrajFile(dcd_names[0].c_str());
	} else {
		fprintf(stderr,"unknown trajectory file type:\n");
		fprintf(stderr,"%s\n",dcd_names[0].c_str());
		exit(1);
	}
	natom = inptraj->num_atoms();
	fprintf(outp,"# ");
	for (int i=m; i<natom; i+=m) {
		for (int j=0; j<i; j+=m) {
			fprintf(outp," %i:%i",i+1,j+1);
		}
	}
	fprintf(outp,"\n");

	inptraj->close();
	float *X = new float[natom];
	float *Y = new float[natom];
	float *Z = new float[natom];
	int frame = 0;
	nu = 0;
	for (int f=0; f<nfile; f++) {
		int slen = dcd_names[f].size();
		if (slen-dcd_names[f].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			inptraj = new DCDITrajFile(dcd_names[f].c_str());
		} else if (slen-dcd_names[f].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			inptraj = new XTCITrajFile(dcd_names[f].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",dcd_names[f].c_str());
			exit(1);
		}
		inptraj->open(dcd_names[f].c_str());
		while (inptraj->frames_left()) {
			if (frame %1000 == 0) {
				fprintf(stderr,"frame %i\n", frame);
			}
			frame++;
			inptraj->read_frame(X,Y,Z,natom);
			if (frame % skip == 0) {
				for (int i=m; i<natom; i+=m) {
					for (int j=0; j<i; j+=m) {
						dx = X[i] - X[j];
						dy = Y[i] - Y[j];
						dz = Z[i] - Z[j];
						if (boxL>0.) {
							dx-=boxL*round(dx/boxL);
							dy-=boxL*round(dy/boxL);
							dz-=boxL*round(dz/boxL);
						}
						dr = sqrt(dx*dx+dy*dy+dz*dz);
						fprintf(outp," %8.3f",dr);
					}
				}
				fprintf(outp,"\n");
			}
		}
		inptraj->close();
	}
	delete [] X;
	delete [] Y;
	delete [] Z;
	return;
}

int main(int argc, char **argv)
{
	vector<string> dcd_names, sel_names;
	vector<double> con;
	double cutoff;
	string outp_name;
	int nu;
	int c, natom, ntraj, dim;
	// defaults
	cutoff = 6.0;
	float boxL = -1;
	if (argc == 1) {
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	int skip = 1;
	int m=1;
	while (1) {
		c=getopt(argc,argv,"hc:o:s:L:m:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stderr,"%s\n",Usage);
				exit(0);
				break;
			case 'L':
				boxL = atof(optarg);
				break;
			case 'c':
				cutoff = atof(optarg);
				break;
			case 's':
				skip = atoi(optarg);
				break;
			case 'm':
				m = atoi(optarg);
				break;
			case 'o':
				outp_name = string(optarg);
				break;
			//case 'p':
			//	pdb_name = optarg;
			//	break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	FILE *outp = fopen(outp_name.c_str(),"w");

	ntraj = (argc-optind);	
	dcd_names.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		dcd_names[i] = argv[optind+2*i];
	}
	dump_dists(dcd_names, cutoff, outp,skip,m,boxL);
	return 0;
}


