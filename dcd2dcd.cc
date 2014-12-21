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
		out.write_frame(X,Y,Z,natom);
	}

	inp.close();
	out.close();

	delete [] X;
	delete [] Y;
	delete [] Z;

	return 0;
}

