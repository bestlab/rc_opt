/*
 * run through trajectory with defined folded and unfolded states,
 * and identify those contacts that DIFFER most between the two
 * i.e. this should be useful when there are substantial non-native
 * contacts in the unfolded state
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <fstream>
#include "TrajFile.h"

using namespace std;

const char *Usage = "\n\n\
	ave_contacts -c cut trj.dcd selection.dat\n\
	cut is a contact cutoff in A\n\
	output goes to stdout\n\
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

void write_conts(vector<double> &con)
{
	int i,j;
	int dim = con.size();
	fprintf(stdout,"%-6s %-6s %-8s\n", "I","J","Q");
	for (int p=0; p<dim; p++) {
		get_ij(p,i,j);
		fprintf(stdout,"%6i %6i %8.5f\n", i,j,con[p]);
	}
	return;
}


void ave_conts(vector<string> &dcd_names, vector<string> &istp_names,
		vector<double> &con, double rcut)
{
	//DCDITrajFile inptraj;
	BaseITrajFile *inptraj;
	ifstream inptp;
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
	//inptraj.open(dcd_names[0].c_str());
	natom = inptraj->num_atoms();
	inptraj->close();
	dim = natom*(natom-1)/2;
	con.resize(dim);
	int *icon = new int[dim];
	for (int i=0; i<dim; i++) {
		icon[i] = 0;
	}
	float *X = new float[natom];
	float *Y = new float[natom];
	float *Z = new float[natom];
	int frame = 0;
	nu = 0;
	for (int f=0; f<nfile; f++) {
		int slen = dcd_names[i].size();
		if (slen-dcd_names[i].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			inptraj = new DCDITrajFile(dcd_names[i].c_str());
		} else if (slen-dcd_names[i].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			inptraj = new XTCITrajFile(dcd_names[i].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",dcd_names[i].c_str());
			exit(1);
		}
		inptraj->open(dcd_names[f].c_str());
		inptp.open(istp_names[f].c_str());
		while (inptraj->frames_left()) {
			if (frame %1000 == 0) {
				fprintf(stderr,"frame %i\n", frame);
			}
			frame++;
			inptraj->read_frame(X,Y,Z,natom);
			inptp >> istp;
			for (int i=1; i<natom; i++) {
				for (int j=0; j<i; j++) {
					dx = X[i] - X[j];
					dy = Y[i] - Y[j];
					dz = Z[i] - Z[j];
					dr = sqrt(dx*dx+dy*dy+dz*dz);
					if (dr<rcut) {
						idx = index(i,j);
						if (istp==1) {
							icon[idx]++;
						} 
					}
				}
			}
			if (istp==1) {
				nu++;
			} 
		}
		inptraj.close();
		inptp.close();
		inptp.clear();
	}
	for (int i=0; i<dim; i++) {
		nu > 0 ? con[i] = double(icon[i])/double(nu) : 0.0;
	}
	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] icon;
	return;
}

int main(int argc, char **argv)
{
	vector<string> dcd_names, sel_names;
	vector<double> con;
	double cutoff;
	int nu;
	int c, natom, ntraj, dim;
	// defaults
	cutoff = 6.0;
	while (1) {
		c=getopt(argc,argv,"hc:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stderr,"%s\n",Usage);
				exit(0);
				break;
			case 'c':
				cutoff = atof(optarg);
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
	ntraj = (argc-optind)/2;	// half dcd, half sel
	dcd_names.resize(ntraj);
	sel_names.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		dcd_names[i] = argv[optind+2*i];
		sel_names[i] = argv[optind+2*i+1];
	}
	ave_conts(dcd_names, sel_names, con, cutoff);
	write_conts(con);
	return 0;
}


