/*
   file - trunctraj.cc
   Author - Robert Best

   Purpose - simple program to extract a single frame from a trajectory
   Usage - pull_frame <n> <crd_file> <dcd_file> 
 */
	
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstdlib>
#include "TrajFile.h"
#include "config.h"

#define DEBUG

string usage = "Usage: trunctraj <first> <last> <dcdin> <dcdout>:\n\
The input trajectory will be written to unit dcdout, but\n\
with frames selected between first and last (inclusive).\n";

int main ( int argc, char * argv [] )
{
	string idcd, odcd, pdb;
	double shapemat[6];
	DCDITrajFile inptraj;
	int natom;
	int firstframe, lastframe;
	
	// parse command line and input file
	if ( argc < 5 )
	{
		cerr << "\nIncorrect number of command line args\n\n"
			<< usage << endl;
		exit(1);
	}

	int arg = 1;
	
	firstframe = atoi( argv[1] );
	lastframe = atoi( argv[2] );
	idcd = string( argv[3] );
	odcd = string( argv[4] );

	inptraj.open( idcd.c_str() );
	natom = inptraj.num_atoms(); 
	int nframe = inptraj.total_frames();
	if (lastframe > nframe) {
		cerr << "Even if this trajectory is complete, it would" <<endl;
		cerr << "only containt " << nframe << " frames;" <<endl;
		cerr << "your last frame is " << lastframe  <<endl;
	}
	float *xx = new float[natom];
	float *yy = new float[natom];
	float *zz = new float[natom];
	DCDOTrajFile outtraj;
	outtraj.setup(inptraj);
	outtraj.set_frames(lastframe-firstframe+1);
	outtraj.open(odcd.c_str());
	outtraj.write_header();
	inptraj.skip(firstframe-1);
	for (int i=firstframe; i<=lastframe; i++) {
		inptraj.read_frame(xx,yy,zz,natom);
		outtraj.write_frame(xx,yy,zz,natom);
	}
	inptraj.close();
	outtraj.close();
	delete [] xx;
	delete [] yy;
	delete [] zz;
		
	return 0;
}
