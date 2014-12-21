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

const char *Usage = "\n\n\
	ncdiff_contacts -c cut trj1.dcd istp1.dcd ... trjN.dcd istpN.dcd\n\
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

void write_conts(vector<double> &ucon,vector<double> &fcon,vector<double> &tpcon)
{
	int i,j;
	int dim = ucon.size();
	fprintf(stdout,"%-6s %-6s %-8s %-8s %-8s\n", "I","J","A","B","TP");
	for (int p=0; p<dim; p++) {
		get_ij(p,i,j);
		fprintf(stdout,"%6i %6i %8.5f %8.5f %8.5f\n", i,j,ucon[p],fcon[p],tpcon[p]);
	}
	return;
}


void ave_conts(vector<string> &dcd_names, vector<string> &istp_names,
		vector<double> &ucon, vector<double> &fcon, vector<double> &tpcon, double rcut)
{
	DCDITrajFile inptraj;
	ifstream inptp;
	int natom, nfile, dim, idx, istp;
	int nf,nu,ntp;
	double dx, dy, dz, dr;
	nfile = dcd_names.size();
	inptraj.open(dcd_names[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	dim = natom*(natom-1)/2;
	ucon.resize(dim);
	fcon.resize(dim);
	tpcon.resize(dim);
	int *iucon = new int[dim];
	int *ifcon = new int[dim];
	int *itpcon = new int[dim];
	for (int i=0; i<dim; i++) {
		iucon[i] = 0;
		ifcon[i] = 0;
		itpcon[i] = 0;
	}
	float *X = new float[natom];
	float *Y = new float[natom];
	float *Z = new float[natom];
	int frame = 0;
	nu = nf = ntp = 0;
	for (int f=0; f<nfile; f++) {
		inptraj.open(dcd_names[f].c_str());
		inptp.open(istp_names[f].c_str());
		while (inptraj.frames_left()) {
			if (frame %1000 == 0) {
				fprintf(stderr,"frame %i\n", frame);
			}
			frame++;
			inptraj.read_frame(X,Y,Z,natom);
			inptp >> istp;
			for (int i=1; i<natom; i++) {
				for (int j=0; j<i; j++) {
					dx = X[i] - X[j];
					dy = Y[i] - Y[j];
					dz = Z[i] - Z[j];
					dr = sqrt(dx*dx+dy*dy+dz*dz);
					if (dr<rcut) {
						idx = index(i,j);
						if (istp==0) {
							iucon[idx]++;
						} else if (istp==2) {
							ifcon[idx]++;
						} else {
							itpcon[idx]++;
						}
					}
				}
			}
			if (istp==0) {
				nu++;
			} else if (istp==2) {
				nf++;
			} else {
				ntp++;
			}
		}
		inptraj.close();
		inptp.close();
		inptp.clear();
	}
	for (int i=0; i<dim; i++) {
		nu > 0 ? ucon[i] = double(iucon[i])/double(nu) : 0.0;
		nf > 0 ? fcon[i] = double(ifcon[i])/double(nf) : 0.0;
		ntp > 0 ? tpcon[i] = double(itpcon[i])/double(ntp) : 0.0;
	}
	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] iucon;
	delete [] ifcon;
	delete [] itpcon;
	return;
}

int main(int argc, char **argv)
{
	vector<string> dcd_names, istp_names;
	vector<double> ucon, fcon, tpcon;
	double cutoff;
	int nf,nu,ntp;
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
	ntraj = (argc-optind)/2;	// half dcd, half istp
	dcd_names.resize(ntraj);
	istp_names.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		dcd_names[i] = argv[optind+2*i];
		istp_names[i] = argv[optind+2*i+1];
	}
	ave_conts(dcd_names, istp_names, ucon, fcon, tpcon, cutoff);
	write_conts(ucon,fcon,tpcon);
	return 0;
}


