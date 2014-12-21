/*
 * count number of "correct" or incorrect residues from
 * Go-model simulation
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "Coorcharmm.h"
#include "TrajFile.h"

const char * Usage = "\n\n\
	Usage:\n\
\n\
	go_correct -h \n\
		for help \n\
	go_correct -s sigma_phi -p xyz.pdb [-q conts.dat] -o outp.dat f1.dcd ...fN.dcd\n\
	where:\n\
		-s sigma_phi is the width of the gaussian for counting native dihedrals\n\
		-p xyz.pdb is the pdb file to take the native dihedral values from\n\
		-o outp.dat is the file the reaction coord will be written to\n\
		f1...fN.dcd are the trajectory files\n\
		[-q conts.dat] (optional) contact list for computing intersegment contacts\n\
		\n\
		The coordinate is sum_i A*exp(-((phi_nat(i)-phi(i))/sigma)^2)\n\n";

double delta_dihe(double dihe1, double dihe2) 
{
	double max, min;
	double xx, yy;
	if (dihe1 > dihe2) {
		max = dihe1;
		min = dihe2;
	} else {
		max = dihe2;
		min = dihe1;
	}
	xx = abs(max-min-360);
	yy = abs(max-min);
	return xx < yy ? xx : yy;
}

int main(int argc, char **argv)
{
	const int BUF_LEN = 1024;
	char pdb_file[BUF_LEN], outp_file[BUF_LEN], tmps[BUF_LEN];
	vector<string> traj_files;
	vector<double> phi_nat;
	double S, sigma, phi,tmpf;
	int c, ntraj, natom, ndihe, cum_fr;
	FILE *outp;
	DCDITrajFile inptraj;

	sigma = 15.0;
	strcpy(outp_file,"def.out");
	strcpy(pdb_file,"NULL");
	while (1) {
		c=getopt(argc,argv,"hs:p:o:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 's':
				sigma = atof(optarg);
				break;
			case 'p':
				strcpy(pdb_file,optarg);
				break;
			case 'o':
				strcpy(outp_file,optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	if (strcmp(pdb_file, "NULL")==0) {
		fprintf(stderr, "Must specify pdb file with -p\n");
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	ntraj=argc-optind;
	if (ntraj==0) {
		fprintf(stderr, "Must specify at least one trajectory\n");
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	traj_files.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		traj_files[i] = string(argv[optind+i]);
	}

	Coorcharmm crd(pdb_file,"pdb");
	natom = crd.Get_N();
	ndihe = natom-3;
	phi_nat.resize(ndihe);
	for (int i=0; i<ndihe; i++) {
		phi_nat[i] = crd.dihedral_angle(i,i+1,i+2,i+3);
	}
	outp = fopen(outp_file,"w");
	if (outp == NULL) {
		fprintf(stderr,"Could not open output file %s\n",outp_file);
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	cum_fr = 0;
	for (int t=0; t<ntraj; t++) {
		inptraj.open(traj_files[t].c_str());
		while (inptraj.frames_left()) {
			inptraj >> crd;
			S = 0;
			for (int i=0; i<ndihe; i++) {
				phi = crd.dihedral_angle(i,i+1,i+2,i+3);
				tmpf = delta_dihe(phi,phi_nat[i])/sigma;
				S += exp(-tmpf*tmpf);
			}
			cum_fr++;
			fprintf(outp,"%12i %12.5f\n", cum_fr, S);
		}
		inptraj.close();
	}
}
