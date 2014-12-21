/* 
 * given a pdb file, a set of trajectories, and a list of "native"
 * contacts and cutoffs, compute the fraction of native contacts
 * over the trajectories.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TrajFile.h"

string usage = "\n\n	Usage:\n"
"	qcalc go.qlist ncout.dat go1.dcd ... goN.dcd\n"
"		where\n"
"		* go.qlist is the input native contact list (pairs of\n"
"			atom numbers and 'native' distances)\n"
"		* ncout.dat is the file to write the fraction of native\n"
"			contacts during simulations\n"
"			including 'native' distances\n"
"		* goi.dcd are the simulation trajectories\n\n";

void ReadContacts(const string &file, vector<int> &i, vector<int> &j,
		vector<double> &rij)
{
	ifstream inp(file.c_str());
	const int bufsz = 1024;
	char buf[bufsz];
	int p,q,nc;
	nc = 0;
	while (inp.good()) {
		inp.getline(buf,bufsz,'\n');
		nc++;
	}
	nc--;
	inp.close();
	i.resize(nc);
	j.resize(nc);
	rij.resize(nc);
	ifstream inp2(file.c_str());
	for (int t=0; t<nc; t++) {
		inp2 >> i[t] >> j[t] >> rij[t];
		//cout << i[t] << j[t] << endl;
	}
	inp2.close();
	return;
}

int main( int argc, char ** argv )
{
	string pdbfile, qlist, output;
	vector<string> dcd_names;
	int c,nselect,ii,jj,natom;
	int * selection;
	vector<int> i, j;
	vector<double> rij;
	double com[3], mat[3][3], dr,dx,dy,dz;
	double gamma;

	gamma = 5.0;

	while (1) {
		c=getopt(argc,argv,"hg:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage.c_str());
				exit(0);
				break;
			case 'g':
				gamma = atof(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage.c_str());
				exit(1);
		}
	}

	if (argc - optind < 3) {
		cout << usage << endl;
		exit(0);
	}
	qlist = string(argv[optind]);
	output = string(argv[optind+1]);
	dcd_names.resize(argc-optind-2);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[optind+2+t]);
	}
	ReadContacts(qlist,i,j,rij);
	// holds reference coordinates

	DCDITrajFile inptraj;
	inptraj.open(dcd_names[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	double *X = new double[natom];
	double *Y = new double[natom];
	double *Z = new double[natom];
	int frames = 0;
	// first pass to compute mean coordinates
	ofstream outp(output.c_str());
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		//int no_frames = inptraj.total_frames();
		int no_frames = inptraj.actual_frames();
		for (int frame = 0; frame < no_frames; frame++) {
		//for (int frame = 0; frame < 1; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			float Q = 0.0;
			int N = 0;
			for (int t=0; t<i.size(); t++) {
				ii = i[t]-1;
				jj = j[t]-1;
				dx = X[ii]-X[jj];
				dy = Y[ii]-Y[jj];
				dz = Z[ii]-Z[jj];
				dr = sqrt(dx*dx+dy*dy+dz*dz);
				double q = 1.0/(1.0+exp(gamma*(dr-rij[t])));
				//fprintf(stdout,"%5i %5i %5i %8.3f %8.3f %8.3f %8.3f %8.3f\n",t,ii,jj,dr,rij[t],
				//		q,X[ii],X[jj]);
				if ( dr < rij[t] ) 
					N++;
			}
			Q = float(N)/float(i.size());
			frames++;
			outp << frames << "\t" << Q << "\t" << N << endl;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

