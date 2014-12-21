/* 
 * calculate drms
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>

#include "TrajFile.h"

const char *usage = "\n\n"
"	Usage:\n"
"		drms [-L boxlen] -q contacts.dat -o ncout.dat trj1.dcd ... trjN.dcd\n"
"		drms [-L boxlen] [-b|-a] [-s i1-j1:i2-j2...] -p xyz.pdb -c cut \n"
"				-o ncout.dat trj1.dcd ... trjN.dcd\n"
"\n"
"		In the first example, contacts are read from contacts.dat and\n"
"	used to compute dRMS from the trajectories, writing output\n"
"	to ncout.dat\n"
"		In the second example, contacts are computed from the\n"
"	coordinates in xyz.pdb using cut-off cut and the rest is as above;\n"
"      	-s 11-20:33-45:56-72 will only use the selected residues\n"
"      	-b will only use backbone atoms (CA,C,N,O)\n"
"      	-a [default] will use all atoms (note that hydrogens (H*) are\n"
"			always ignored).\n"
"		In either case an optional box length argument (\"-L\") will calculate\n"
"	distances using a cubic box.\n\n"
"\n"
"	Old Usage (still works):\n"
"		drms contacts.dat ncout.dat trj1.dcd ... trjN.dcd\n"
"	where\n"
"		* contacts.dat is a reference contact list (i,j,rij) \n"
"		* ncout.dat is the file to write drms to\n"
"		* trji.dcd are the simulation trajectories\n\n";


void calc_contacts_pdb(const string &pdb_file_name, vector<int> &i, vector<int> &j, 
		vector<double> &rij, double cut, const string &selection, bool bb_only)
{
	// TODO
}

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
	string pdb_file_name, qlist_file_name, output_name, selection;
	double DRMS;
	int c, dcd_off;
	bool bb_only;
	vector<string> dcd_names;
	int nselect,ncontact,ii,jj,natom;
	int * sel;
	vector<int> i, j;
	vector<double> rij;
	double com[3], mat[3][3], dr,dx,dy,dz;
	double box_len = -1.;
	double cut = 6.0;
	bb_only = false;
	qlist_file_name = "NONE";
	output_name = "drms.dat";
	pdb_file_name = "NONE";
	selection = "ALL";


	// first try parse options
	while (1) {
		c=getopt(argc,argv,"hL:q:o:p:c:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage);
				exit(0);
				break;
			case 'L':
				box_len = atof(optarg);
				break;
			case 'q':
				qlist_file_name = optarg;
				break;
			case 'o':
				output_name = optarg;
				break;
			case 'p':
				pdb_file_name = optarg;
				break;
			case 'c':
				cut = atof(optarg);
				break;
			case 's':
				selection = optarg;
				break;
			case 'b':
				bb_only = true;
				break;
			case 'a':
				bb_only = false;
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}
	
	// then fall back on legacy mode!
	if ( qlist_file_name == string("NONE") ) {
		if (argc < 4) {
			cout << usage << endl;
			exit(0);
		}
		qlist_file_name = string(argv[1]);
		output_name = string(argv[2]);
		dcd_off = 3;
		ReadContacts(qlist_file_name,i,j,rij);
	} else { 
		dcd_off = optind;
		// compute contact list from pdb file
		calc_contacts_pdb(pdb_file_name, i, j, rij, cut, selection, bb_only);
	}
	dcd_names.resize(argc-dcd_off);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[dcd_off+t]);
	}
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
	ofstream outp(output_name.c_str());
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
			for (int t=0; t<i.size(); t++) {
				ii = i[t]-1;
				jj = j[t]-1;
				dX = X[ii]-X[jj];
				dY = Y[ii]-Y[jj];
				dZ = Z[ii]-Z[jj];
				dR = sqrt(dX*dX+dY*dY+dZ*dZ);
				ddR = dR - rij[t];
				DRMS += ddR*ddR;
			}
			frames++;
			outp << frames << "\t" << sqrt(DRMS/float(i.size())) << endl;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

