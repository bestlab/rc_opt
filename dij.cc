/* 
 * calculate drms
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TrajFile.h"

string usage = "\n\n\
	Usage:\n\
		drms contacts.dat ncout.dat trj1.dcd ... trjN.dcd\n\
	where\n\
		* contacts.dat is a reference contact list (i,j,rij) \n\
		* ncout.dat is the file to write drms to\n\
		* trji.dcd are the simulation trajectories\n\n";

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

double dRMS(const vector<int> &i, const vector<int> &j, const vector<double> &rij, 
		double *X, double *Y, double *Z)
{
	int ncontact = i.size();
	int ii,jj;
	double dX, dY, dZ, dR, ddR;
	double drms;
	drms = 0.0;
	//cout << i[0] << "\t" << j[0]  << "\t" << rij[0] << X[0]  << endl;
	for (int p=0; p<ncontact; p++) {
		ii = i[p];
		jj = j[p];
		dX = X[ii]-X[jj];
		dY = Y[ii]-Y[jj];
		dZ = Z[ii]-Z[jj];
		dR = sqrt(dX*dX+dY*dY+dZ*dZ);
		ddR = dR - rij[p];
		drms += ddR*ddR;
	}
	//cout << sqrt(drms/double(ncontact))  << endl;
	return sqrt(drms/double(ncontact));
}


int main( int argc, char ** argv )
{
	string pdbfile, qlist, output;
	double DRMS;
	vector<string> dcd_names;
	int nselect,ncontact,ii,jj,natom;
	int * selection;
	int I,J;
	vector<int> i, j;
	vector<double> rij;
	double com[3], mat[3][3], dr,dx,dy,dz;

	if (argc < 4) {
		cout << usage << endl;
		exit(0);
	}
	I = atoi(argv[1]);
	J = atoi(argv[2]);
	output = string(argv[3]);
	dcd_names.resize(argc-4);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[4+t]);
	}
	//ReadContacts(qlist,i,j,rij);
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
	double dR,dX,dY,dZ,ddR;
	ncontact=i.size();
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		int no_frames = inptraj.total_frames();
		for (int frame = 0; frame < no_frames; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			//DRMS = dRMS(i,j,rij,X,Y,Z);
			DRMS = 0.0;
			dX = X[I-1]-X[J-1];
			dY = Y[I-1]-Y[J-1];
			dZ = Z[I-1]-Z[J-1];
			dR = sqrt(dX*dX+dY*dY+dZ*dZ);
			frames++;
			outp << frames << "\t" << dR << endl;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

