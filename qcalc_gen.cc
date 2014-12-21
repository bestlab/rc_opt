/* 
 * Compute a general contact based reaction coordinate
 * based on the following definition of Q:
 *
 * Q = sum_(i,j) w_ij H(r_ij,rcut_ij)
 *
 * where w_ij is a weight factor
 * H(r,rc) is a step function of form H(r,rc) = 1.0/(1.0+exp(beta*(r-rc)))
 * r_ij is the distance between atoms i and j
 * rcut_ij is a cutoff distance specific to the pair (i,j)
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TrajFile.h"

#define DEBUG

string usage = "\n\n	Usage:\n		qcalc go.qlist ncout.dat go1.dcd ... goN.dcd\n	where\n		* go.qlist is the input native contact list (pairs of\n		  atom numbers and 'native' distances)\n		* ncout.dat is the file to write the fraction of native\n		  contacts during simulations\n		  including 'native' distances\n		* goi.dcd are the simulation trajectories\n\n";

/*
string usage = "\n\n	Usage:\n		qcalc_gen go.qlist ncout.dat go1.dcd ... goN.dcd\n	where\n		* go.qlist is the input native contact list (pairs of\n		  atom numbers, 'native' distances and weights)\		* ncout.dat is the file to write the fraction of native		  contacts during simulations\n		  including 'native' distances\n		* goi.dcd are the simulation trajectories\n\n";
*/

void ReadContacts(const string &file, vector<int> &i, vector<int> &j,
		vector<double> &rcut_ij, vector<double> &w_ij)
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
	rcut_ij.resize(nc);
	w_ij.resize(nc);
	ifstream inp2(file.c_str());
	for (int t=0; t<nc; t++) {
		inp2 >> i[t] >> j[t] >> rcut_ij[t] >> w_ij[t];
	}
	inp2.close();
	return;
}

int main( int argc, char ** argv )
{
	string qlist, output;
	vector<string> dcd_names;
	int nselect,natom,ii,jj;
	int * selection;
	vector<int> i, j;
	vector<double> rcut_ij, w_ij;
	double com[3], mat[3][3], dr,dx,dy,dz,beta;

	beta = 5.0;
	if (argc < 4) {
		cout << usage << endl;
		exit(0);
	}
	qlist = string(argv[1]);
	output = string(argv[2]);
	dcd_names.resize(argc-3);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[3+t]);
	}
	ReadContacts(qlist,i,j,rcut_ij,w_ij);
	// holds reference coordinates
	//Coorcharmm coor_main (pdbfile.c_str(),"pdb");
	//coor_main.set_uniform_weights();	// = 1

	//coor_main.center_of_mass(com);
	//coor_main.translate( -com[0], -com[1], -com[2] );
	//Coorcharmm coor_work = coor_main;
	DCDITrajFile inptraj;
	inptraj.open(dcd_names[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	double *X = new double[natom];
	double *Y = new double[natom];
	double *Z = new double[natom];
	int frames = 0;
	ofstream outp(output.c_str());
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		int no_frames = inptraj.total_frames();
		for (int frame = 0; frame < no_frames; frame++) {
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
				double q = w_ij[t]/(1.0+exp(beta*(dr-rcut_ij[t])));
				Q += q;
			}
			frames++;
			outp << Q << endl;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

