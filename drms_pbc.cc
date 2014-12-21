/* 
 * calculate drms
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>

#include "TrajFile.h"

string usage = "\n\n"
"Usage:"
"\n		drms_pbc boxlen beta go.qlist ncout.dat go1.dcd ... goN.dcd\n"
"where\n"
"* go.qlist is the input native contact list (pairs of\n"
"atom numbers and 'native' distances)\n"
"* ncout.dat is the file to write the fraction of native\n"
"contacts during simulations\n"
"including 'native' distances\n"
"* goi.dcd are the simulation trajectories\n\n";

inline double anint( double x )
{
	return ceil( x - 0.50 );
}

void read_pdb(double *X, double *Y, double *Z, int nres, const char *f)
{
        const int buf_len = 1024;
        char buf[buf_len], buf2[buf_len], buf3[buf_len];
	int i;
	float x,y,z;

	FILE *pdb = fopen(f, "r");
	i = 0;
	while (feof(pdb) == 0 && i<nres) {
                fgets(buf,buf_len,pdb);
		if (strncmp(buf,"ATOM",4) == 0) {
			//fprintf(stdout,"%s\n",buf);
			sscanf(buf,"%30c%8f%8f%8f",
					buf2,&x,&y,&z);
			//sscanf(buf,"%30c\n",
			//		buf2);
			X[i] = x;
			Y[i] = y;
			Z[i] = z;
			i++;
		}
	}
}

double compute_dmat(double *X,double *Y, double *Z, int n, int na,
		vector<vector<double> > &dmat, double boxlen)
/* returns shortest pair distance between A and B */
{
	double dx,dy,dz,r,aint_x_l,aint_y_l,aint_z_l,rmin;
	int nb = n-na;
	dmat.resize(na);
	for (int q=0;q<na;q++)
		dmat[q].resize(nb);
	rmin = -1.;
	for (int i=0; i<na;i++) {
		for (int j=0;j<nb;j++) {
			int jj = na+j;
			dx = X[i] - X[jj];
			dy = Y[i] - Y[jj];
			dz = Z[i] - Z[jj];
			aint_x_l = anint( (dx)/boxlen );
			aint_y_l = anint( (dy)/boxlen );
			aint_z_l = anint( (dz)/boxlen );
			dx -=  boxlen*aint_x_l;
			dy -=  boxlen*aint_y_l;
			dz -=  boxlen*aint_z_l;
			r = sqrt(dx*dx+dy*dy+dz*dz);
			if (r<rmin || rmin < 0.)
				rmin = r;
			dmat[i][j] = r;
		}
	}
	return rmin;
}

double drms(vector<vector<double> > &a, vector<vector<double> > &b)
{
	double tot  = 0.;
	int nx = a.size();
	int ny = a[0].size();
	for (int i=0;i<nx;i++)
		for (int j=0; j<ny; j++) {
			tot += abs(a[i][j]-b[i][j]);
		}
	return tot/float(nx*ny);
}
	
int main( int argc, char ** argv )
{
	string pdbfile, qlist, output;
	vector<string> dcd_names;
	int nselect,ii,jj,natom;
	int * selection;
	vector<int> i, j;
	vector<double> rij;
	vector<vector<double> > dmat_nat, dmat_cur;
	double beta;
	double aint_x_l,aint_y_l,aint_z_l;
	double com[3], mat[3][3], dr,dx,dy,dz,boxlen;
	string pdb;
	int piv;

	if (argc < 4) {
		cout << usage << endl;
		exit(0);
	}
	boxlen = atof(argv[1]);
	pdb = string(argv[2]);
	piv = atoi(argv[3]);
	output = string(argv[4]);
	dcd_names.resize(argc-5);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[5+t]);
	}
	DCDITrajFile inptraj;
	inptraj.open(dcd_names[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	double *X = new double[natom];
	double *Y = new double[natom];
	double *Z = new double[natom];
	read_pdb(X, Y, Z, natom, pdb.c_str());
	compute_dmat(X,Y,Z,natom,piv,dmat_nat,boxlen);
	
	int frames = 0;
	double DRMS,RMIN;

	// first pass to compute mean coordinates
	ofstream outp(output.c_str());
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		int no_frames = inptraj.total_frames();
		double xtal[6];
		for (int frame = 0; frame < no_frames; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			if (inptraj.crystal()) {
				inptraj.get_crystal_data(xtal);
				boxlen = xtal[0];
			}
			RMIN = compute_dmat(X,Y,Z,natom,piv,dmat_cur,boxlen);
			DRMS = drms(dmat_nat, dmat_cur);
			frames++;
			outp << frames << "\t" << DRMS << "\t" << RMIN << endl;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

