/* 
 * apply kmt reduction to determine if there is a knot in a chain
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TrajFile.h"

string usage = "\n\n	Usage:\n		kmt skip go1.dcd ... goN.dcd\n"
"	where\n"
"                  skip is an integer value (only analyze every skip frames)\n"
"		   go1.dcd ... goN.dcd are trajectory files\n\n";


double dist(vector<double> &p, vector<double>&q )
{
	double dx, dy, dz;
	dx = p[0]-q[0];
	dy = p[1]-q[1];
	dz = p[2]-q[2];
	return sqrt(dx*dx+dy*dy+dz*dz);
}

bool intriangle(vector<double> &I, vector<double> &J, vector<double> &K,
		vector<double> &P, vector<double> &Q)
{
	vector<double> a(3),b(3),N(3),S(3),L(3);
	double num,den;
	double big,small;
	// (i) find S, the point on line PQ and in plane IJK
	a[0] = J[0]-I[0];
	a[1] = J[1]-I[1];
	a[2] = J[2]-I[2];
	b[0] = J[0]-K[0];
	b[1] = J[1]-K[1];
	b[2] = J[2]-K[2];
	N[0] = a[1]*b[2]-a[2]*b[1];
	N[1] = a[2]*b[0]-a[0]*b[2];
	N[2] = a[0]*b[1]-a[1]*b[0];
	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",N[0],N[1],N[2]); 
	num = N[0]*J[0]+N[1]*J[1]+N[2]*J[2]-N[1]*P[1]-N[2]*P[2]
		+ N[1]*(Q[1]-P[1])/(Q[0]-P[0])*P[0] + N[2]*(Q[2]-P[2])/(Q[0]-P[0])*P[0];
	den = N[0]+N[1]*(Q[1]-P[1])/(Q[0]-P[0])+N[2]*(Q[2]-P[2])/(Q[0]-P[0]);
	//fprintf(stdout,"%8.3f %8.3f\n",num,den); 
	S[0] = num/den;
	S[1] = P[1]+(Q[1]-P[1])*(S[0]-P[0])/(Q[0]-P[0]);
	S[2] = P[2]+(Q[2]-P[2])*(S[0]-P[0])/(Q[0]-P[0]);
	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",S[0],S[1],S[2]); 
	//
	// (ii) is S between P and Q!
	if (abs(S[0]-P[0])-abs(Q[0]-P[0])>=0) //>-1.e-15) // S not even between P and Q
		return false;
	//
	// (iii) is S in the triangle IJK?
	// find intersection L of IK and JS
	num = J[0]*(S[1]-J[1])/(S[0]-J[0])-I[0]*(K[1]-I[1])/(K[0]-I[0])-J[1]+I[1];
	den = (S[1]-J[1])/(S[0]-J[0])-(K[1]-I[1])/(K[0]-I[0]);
	//fprintf(stdout,"%8.3f %8.3f\n",num,den); 
	L[0] = num/den;
	L[1] = J[1] + (L[0]-J[0])/(S[0]-J[0])*(S[1]-J[1]);
	L[2] = J[2] + (L[0]-J[0])/(S[0]-J[0])*(S[2]-J[2]);
	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",L[0],L[1],L[2]); 
	big = dist(J,L);
	small = dist(J,S);
	//fprintf(stdout,"small, big = %8.3f %8.3f\n",small,big); 
	//if (small+1.e-15>=big)
	if (small>=big)
		return false;
	// find intersection L of JK and IS
	num = I[0]*(S[1]-I[1])/(S[0]-I[0])-J[0]*(K[1]-J[1])/(K[0]-J[0])-I[1]+J[1];
	den = (S[1]-I[1])/(S[0]-I[0])-(K[1]-J[1])/(K[0]-J[0]);
	//fprintf(stdout,"%8.3f %8.3f\n",num,den); 
	L[0] = num/den;
	L[1] = I[1] + (L[0]-I[0])/(S[0]-I[0])*(S[1]-I[1]);
	L[2] = I[2] + (L[0]-I[0])/(S[0]-I[0])*(S[2]-I[2]);
	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",L[0],L[1],L[2]); 
	big = dist(I,L);
	small = dist(I,S);
	//fprintf(stdout,"small, big = %8.3f %8.3f\n",small,big); 
	//if (small+1.e-15>=big)
	if (small>=big)
		return false;
	// find intersection L of IJ and KS
	num = K[0]*(S[1]-K[1])/(S[0]-K[0])-I[0]*(J[1]-I[1])/(J[0]-I[0])-K[1]+I[1];
	den = (S[1]-K[1])/(S[0]-K[0])-(J[1]-I[1])/(J[0]-I[0]);
	//fprintf(stdout,"%8.3f %8.3f\n",num,den); 
	L[0] = num/den;
	L[1] = K[1] + (L[0]-K[0])/(S[0]-K[0])*(S[1]-K[1]);
	L[2] = K[2] + (L[0]-K[0])/(S[0]-K[0])*(S[2]-K[2]);
	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",L[0],L[1],L[2]); 
	big = dist(K,L);
	small = dist(K,S);
	//fprintf(stdout,"small, big = %8.3f %8.3f\n",small,big); 
	if (small>=big)
		return false;
	return true;
}

void kmt(double *X, double *Y, double *Z, int natom, 
		double *X_kmt, double *Y_kmt, double *Z_kmt, int &natom_kmt)
{
	int prevkill,kill;
	bool *alive = new bool[natom];
	int *alive_idx = new int[natom];
	natom_kmt = natom;
	vector<double> A(3),B(3),C(3),D(3),E(3);
	for (int i=0; i<natom; i++) {
		alive[i] = true;
		alive_idx[i] = i;
	}
	//prevkill = kill =0;
	prevkill = -1;
	kill = 100;
	int firsti = 1;
	int failedkills = 0;
	//while (! (kill == 0 and prevkill==0) ) {
	while ( failedkills < 4 ) {
	        //fprintf(stdout,"prevkill, kill = %5i %5i\n",prevkill,kill);
		prevkill = kill;
		kill = 0;
		if (firsti == 1) {
			firsti = 2;
		}else {
			firsti = 1;
		}
		for (int i=firsti; i<natom_kmt-1; i+=2) {
			int aa = alive_idx[i-1];
			int bb = alive_idx[i];
			int cc = alive_idx[i+1];
			A[0] = X[aa];
			A[1] = Y[aa];
			A[2] = Z[aa];
			B[0] = X[bb];
			B[1] = Y[bb];
			B[2] = Z[bb];
			C[0] = X[cc];
			C[1] = Y[cc];
			C[2] = Z[cc];
			bool cankill = true;
			for (int j=0; j<natom_kmt-1; j++) {
				if (j==i || j==i-1) 
					continue;
				int dd = alive_idx[j];
				int ee = alive_idx[j+1];
				D[0] = X[dd];
				D[1] = Y[dd];
				D[2] = Z[dd];
				E[0] = X[ee];
				E[1] = Y[ee];
				E[2] = Z[ee];
				if (intriangle(A,B,C,D,E)) {
					cankill = false;
					break;
				}
			}
			if (cankill) {
				alive[bb] = false;
				kill++;
				//fprintf(stdout,"killing atom %i\n",bb);
			}

		}
		natom_kmt = 0;
		int idx = 0;
		for (int i=0; i<natom; i++) {
			if (alive[i]) {
				alive_idx[idx++]=i;
				natom_kmt++;
			}
		}
		if (kill==0)
			failedkills++;
		else
			failedkills==0;
	}
	for (int i=0; i<natom_kmt; i++) {
		X_kmt[i] = X[alive_idx[i]];
		Y_kmt[i] = Y[alive_idx[i]];
		Z_kmt[i] = Z[alive_idx[i]];
	}

	delete [] alive;
	delete [] alive_idx;
}


int main( int argc, char ** argv )
{
	/*
	 * STUFF BELOW TEST CASES FOR WHETHER LINE DE 
	 * PASSES THROUGH TRIANGLE ABC!
	 *
	vector<double> A(3),B(3),C(3),D(3),E(3);
	// an example in the triangle 
	A[0] = 1.01;
	A[1] = -0.001;
	A[2] = 0.;
	B[0] = 0.001;
	B[1] = 0.02;
	B[2] = 0.005;
	C[0] = 0.03;
	C[1] = 1.1;
	C[2] = 0.006;
	D[0] = 0.20;
	D[1] = 0.25;
	D[2] = -1.;
	E[0] = 0.30;
	E[1] = 0.26;
	E[2] = 0.99;
	// an example not in the triangle 
	A[0] = 1.01;
	A[1] = -0.001;
	A[2] = 0.;
	B[0] = 0.001;
	B[1] = 0.02;
	B[2] = 0.005;
	C[0] = 0.03;
	C[1] = 1.1;
	C[2] = 0.006;
	D[0] = 5.20;
	D[1] = 0.25;
	D[2] = -1.;
	E[0] = 6.30;
	E[1] = 0.26;
	E[2] = 0.99;
	if (intriangle(A,B,C,D,E)) {
		fprintf(stdout,"IN TRIANGLE!\n");
	} else {
		fprintf(stdout,"NOT IN TRIANGLE!\n");
	}
	*/


	string output;
	vector<string> dcd_names;
	int natom,natom_kmt;

	
	if (argc ==1) {
		cout << usage << endl;
		exit(0);
	}
	//output = string(argv[2]);
	int skip = atoi(argv[1]);
	dcd_names.resize(argc-2);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[2+t]);
	}

	DCDITrajFile inptraj;
	inptraj.open(dcd_names[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	double *X = new double[natom];
	double *Y = new double[natom];
	double *Z = new double[natom];
	double *X_kmt = new double[natom];
	double *Y_kmt = new double[natom];
	double *Z_kmt = new double[natom];
	int frames = 0;
	// first pass to compute mean coordinates
	ofstream outp(output.c_str());
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		int no_frames = inptraj.total_frames();
		for (int frame = 0; frame < no_frames; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			if (frames%skip == 0) {
				kmt(X, Y, Z, natom, X_kmt, Y_kmt, Z_kmt, natom_kmt);
				fprintf(stdout,"%5i\n",natom_kmt);
			}
			frames++;
			//outp << frames << "\t" << Q << "\t" << CO << endl;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] X_kmt;
	delete [] Y_kmt;
	delete [] Z_kmt;
	return 0;
}

