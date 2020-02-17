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
#include <fstream>
#include <iostream>

#include "traj.h"

using namespace std;

string usage = "\n\n	Usage:\n"
"	distance -i idxi -j idxj go1.dcd ... goN.dcd\n"
"               Just measures the distance in Angstroms between"
"               atoms/beads with indices idxi and idxj over the trajectories\n\n";


void read_idx(const string &file, vector<int> &i, vector<int> &j)
{
	ifstream inp(file.c_str());
	const int bufsz = 1024;
	char buf[bufsz];
	char * junk;
	int p,q,nc,nfield;
	nc = 0;
	nfield = 0;
	while (inp.good()) {
		inp.getline(buf,bufsz,'\n');
		if (nc==0) {
			junk = strtok(buf," ");
			while (junk != NULL) {
				junk = strtok(NULL," ");
				nfield++;
			}

		}
		nc++;
	}
	fprintf(stderr,"number of fields = %i\n",nfield);
	nc--;
	cerr << nc << endl;
	inp.close();
	if (nfield!=2) {
		fprintf(stderr,
				"unknown number of fields (%i) in idx file\n",
				nfield);
		exit(1);
	}
	i.resize(nc);
	j.resize(nc);
	ifstream inp2(file.c_str());
	for (int t=0; t<nc; t++) {
		inp2 >> i[t] >> j[t];
		cerr << i[t] << " " << j[t]  << endl;
	}
	inp2.close();
	return;
}

int main( int argc, char ** argv )
{
	string pdbfile, qlist, output;
	vector<string> dcd_names;
	int c,nselect,atom_i,atom_j,natom;
	int * selection;
	vector<int> i, j;
	vector<double> rij;
	bool use_idx = false;
	double com[3], mat[3][3], dr,dx,dy,dz;
	string idx_file = "NULL";

	while (1) {
		c=getopt(argc,argv,"hi:j:f:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage.c_str());
				exit(0);
				break;
			case 'i':
				atom_i = atoi(optarg);
				break;
			case 'j':
				atom_j = atoi(optarg);
				break;
			case 'f':
				idx_file = optarg;
				use_idx = true;
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
	//fprintf(stdout,"%s\n",idx_file.c_str()); fflush(stdout);
	if (use_idx) {
		fprintf(stderr,"Using index file %s\n",idx_file.c_str());
		read_idx(idx_file, i, j);
	} else {
		i.resize(1);
		j.resize(1);
		i[0] = atom_i;
		j[0] = atom_j;
	}
	int ndat = i.size();
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
	natom = inptraj->num_atoms();
	inptraj->close();
	//
	float *X = new float[natom];
	float *Y = new float[natom];
	float *Z = new float[natom];
	int frames = 0;
	//resi -= 1;	// index from 0
	//resj -= 1;
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
			frames++;
			inptraj->read_frame(X,Y,Z,natom);
			for (int l=0; l<ndat; l++) {
				atom_i = i[l]-1; 
				atom_j = j[l]-1; 
				dx = X[atom_i]-X[atom_j];
				dy = Y[atom_i]-Y[atom_j];
				dz = Z[atom_i]-Z[atom_j];
				dr = sqrt(dx*dx+dy*dy+dz*dz);
				fprintf(stdout," %12.6f %12.6f %12.6f %12.6f",dr,dx,dy,dz);
			}
			fprintf(stdout,"\n");
		}
		inptraj->close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

