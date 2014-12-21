/* 
 * given a pdb file, a set of trajectories, and a list of "native"
 * contacts and cutoffs, compute the fraction of native contacts
 * over the trajectories.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>

#include "TrajFile.h"

string usage = "\n\n	Usage:\n		qcalc go.qlist ncout.dat go1.dcd ... goN.dcd\n	where\n		* go.qlist is the input native contact list (pairs of\n		  atom numbers and 'native' distances)\n		* ncout.dat is the file to write the fraction of native\n		  contacts during simulations\n		  including 'native' distances\n		* goi.dcd are the simulation trajectories\n\n";

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
	int nselect,ii,jj,natom;
	int * selection;
	int pulli, pullj;
	vector<int> i, j;
	vector<double> rij;
	double com[3], mat[3][3], dr,dx,dy,dz;

	if (argc < 4) {
		cout << usage << endl;
		exit(0);
	}
	qlist = string(argv[1]);
	pulli = atoi(argv[2]);
	pullj = atoi(argv[3]);

	output = string(argv[4]);
	//cout << qlist.c_str() << endl;
	//cout << pulli << '\t' << pullj << endl;
	//cout << output.c_str() << endl;
	dcd_names.resize(argc-5);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[5+t]);
		//cout << dcd_names[t].c_str() << endl;
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
	FILE * outp = fopen(output.c_str(),"w");
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		int no_frames = inptraj.total_frames();
		for (int frame = 0; frame < no_frames; frame++) {
		//for (int frame = 0; frame < 1; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			double Q = 0.0;
			double CO = 0.0;
			double Q2 = 0.0;
			double CO2 = 0.0;
			for (int t=0; t<i.size(); t++) {
				ii = i[t]-1;
				jj = j[t]-1;
				dx = X[ii]-X[jj];
				dy = Y[ii]-Y[jj];
				dz = Z[ii]-Z[jj];
				dr = sqrt(dx*dx+dy*dy+dz*dz);
				double q = 1.0/(1.0+exp(5.0*(dr-rij[t])));
				Q += q;
				CO += q*float(abs(i[t]-j[t]));
			}
			dx = X[pulli]-X[pullj];
			dy = Y[pulli]-Y[pullj];
			dz = Z[pulli]-Z[pullj];
			dr = sqrt(dx*dx+dy*dy+dz*dz);
			CO /= Q;
			Q /= float(i.size());
			frames++;
			fprintf(outp,"%10i %12.6f %8.3f\n",frames,Q,dr);
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

