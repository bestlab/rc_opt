/*
 * write vanilla native dcd files
 */

#include <cstdio>
#include <cstdlib>
#include "TrajFile.h"

int main(int argc, char **argv)
{
	DCDITrajFile inp;
	DCDOTrajFile out;
	int natom,nframes;
	float *X,*Y,*Z;
	double x_O,y_O,z_O;

	inp.open(argv[1]);
	out.setup(inp);
	out.open(argv[2]);
	out.write_header();

	natom = inp.num_atoms();
	nframes = inp.total_frames();
	X = new float[natom];
	Y = new float[natom];
	Z = new float[natom];

	for (int i=0; i<nframes; i++) {
		inp.read_frame(X,Y,Z,natom);

		x_O = X[179]; y_O = Y[179]; z_O = Z[179];
		for (int k=179; k>176;k++) {
			X[k] = X[k-1];
			Y[k] = Y[k-1];
			Z[k] = Z[k-1];
		}
		X[176] = x_O;
		Y[176] = y_O;
		Z[176] = z_O;

		out.write_frame(X,Y,Z,natom);
	}

	inp.close();
	out.close();

	delete [] X;
	delete [] Y;
	delete [] Z;

	return 0;
}

