/*
 * dcdinfo.cc
 * Robert Best
 * Tue Sep 25 17:10:49 BST 2001
 *
 * give basic info about dcd file
 */

#include <iostream>
#include "TrajFile.h"
#include "config.h"

int main(int argc, char ** argv)
{
	bool verb = false;
	DCDITrajFile tmptraj;
	for (int i=1; i< argc; i++) {
		if (!strcmp(argv[i], "-v")) {
			verb=true;
			continue;
		}
		try {
			tmptraj.open(argv[i]);
		} catch ( DCDITrajFile::OpenErr err ) {
			cerr << "\n\nCould not open file "
				<< argv[i] << endl;
		}

		if (tmptraj.actual_frames() != tmptraj.total_frames()) {
			cout << "====================================================" << endl;
			cout << "Actual number of frames != number stated in header!" << endl;
			cout << "Number in header = " << tmptraj.total_frames() << endl;
			cout << "Number in file   = " << tmptraj.actual_frames() << endl;
		}

		cout << "===================================" << endl;
		cout << "INFO FOR: " << argv[i] << endl;
		cout << "===================================" << endl;
		string s; tmptraj.read_title(s);
		cout << "Title: 	" << endl << s << endl;
		cout << "First step: 	" << tmptraj.initial_step() << endl;
		cout << "Total frames: 	" << tmptraj.total_frames() << endl;
		cout << "Step size:	" << tmptraj.step_size() << endl;
		cout << "Total # atoms:	" << tmptraj.num_atoms() << endl;
		cout << "# Fixed atoms	" << tmptraj.fixed_atoms() << endl;
		cout << "CHARMM version	" << tmptraj.charmm_version() << endl;
		cout << "Traj. type	";
		switch(tmptraj.traj_type()) {
			case 'c': cout << "Coordinates" << endl;
				  break;
			case 'v': cout << "Velocities" << endl;
				  break;
			default: cout << "Unknown" << endl;
		}
		cout << "Byte-swapping?	";
		if ( tmptraj.swapping_bytes() ) {
			cout << "Yes" << endl;
		} else {
			cout << "No" << endl;
		}
		cout << "Crystal?	";
		if ( tmptraj.crystal() ) {
			cout << "Yes" << endl;
		} else {
			cout << "No" << endl;
		}
		cout << "# DOF:		" << tmptraj.deg_free() << endl;
		cout << "===================================" << endl
			<< endl;

		if (verb) {	// print info for every frame!
			int N=tmptraj.num_atoms();
			double xtal[6];
			double * X = new double[tmptraj.num_atoms()];
			double * Y = new double[tmptraj.num_atoms()];
			double * Z = new double[tmptraj.num_atoms()];
			int frame=0;
			while ( tmptraj.frames_left() ) {
				tmptraj.read_frame(X, Y, Z, N);
				tmptraj.get_crystal_data(xtal);
				cout << "Frame " << frame << endl;
				cout << "Crystal data:" << endl;
				cout << xtal[0] << '\n' << xtal[1] <<  '\t'
					<< xtal[2] << '\n' << xtal[3]
					<< '\t' << xtal[4] << '\t' 
					<< xtal[5] << '\n' << endl;
				frame++;
			}
		}
		tmptraj.close();
	}

	return 0;
}

