/*
 * extract transition path DCD's from equilibrium traj...
 */

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <fstream>
#include "TrajFile.h"

void read_istp(vector<string> &istp_names, int *isTP, int skip, int nsave)
{
	ifstream inp;
	int nline,idx;
	int tmpi;
	int nfile = istp_names.size();
	idx=0;
	for (int f=0; f<nfile; f++) {
		inp.open(istp_names[f].c_str());
		nline = 0;
		while (!inp.eof()) {
			inp >> tmpi;
			nline++;
		}
		nline--;
		fprintf(stdout,"nline = %i\n",nline);
		inp.close();
		inp.clear();
		inp.open(istp_names[f].c_str());
		for (int l=0; l<nline; l++) {
			inp >> tmpi;
			if (l%skip == 0) {
				isTP[idx++] = tmpi;
			}
		}
		inp.close();
		inp.clear();
	}
	if (idx != nsave) {
		fprintf(stderr,"istp len [%i] != traj len [%i]\n",idx,nsave);
		exit(1);
	} else {
		fprintf(stderr,"GOT HERE: istp len [%i] = traj len [%i]\n",idx,nsave);
	}
}

void write_tp(vector<string> &dcd_names, const char *oname, 
		int tp_start, int tp_end)
{
	fprintf(stdout,"Transition path: %i - %i\n",tp_start,tp_end);
	DCDITrajFile inptraj;
	DCDOTrajFile outtraj;
	inptraj.open(dcd_names[0].c_str());
	inptraj.close();
	int ltraj = inptraj.total_frames();
	int file = tp_start/ltraj;
	inptraj.open(dcd_names[file].c_str());
	inptraj.skip(tp_start-file*ltraj);
	outtraj.setup(inptraj);
	outtraj.open(oname);
	outtraj.set_frames(tp_end-tp_start+1);
	outtraj.write_header();
	int natom = inptraj.num_atoms();
	float *X = new float[natom];
	float *Y = new float[natom];
	float *Z = new float[natom];
	for (int i=tp_start; i<tp_end; i++) {
		if (!inptraj.frames_left()) {
			inptraj.close();
			file++;
			inptraj.open(dcd_names[file].c_str());
		}
		inptraj.read_frame(X,Y,Z,natom);
		outtraj.write_frame(X,Y,Z,natom);
	}
	outtraj.close();
}

int main(int argc, char **argv)
{
	int c,max_save;
	int *isTP;
	int nstep, ntraj;
	const char *Usage = 
		"\nUsage: extract_tp -o oroot dcd1 istp1 ... dcdN istpN\n";
	const int buf_len=1024;
	char outp_root[buf_len] = "default";
	char buf[buf_len];
	vector<string> dcd_names, istp_names;
	int L_eq, L_tp;
	int nc;		// # native conts
	int natom, ncontact;
	DCDITrajFile inptraj;
	int nframe, tot_frame, f_ind;
	
	// parse cmd line options
	while (1) {
		c=getopt(argc,argv,"ho:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'o':
				strcpy(outp_root,optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	// now read dcd names
	ntraj = (argc-optind)/2;	// half dcd, half istp
	dcd_names.resize(ntraj);
	istp_names.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		dcd_names[i] = argv[optind+2*i];
		istp_names[i] = argv[optind+2*i+1];
	}
	max_save = 0;
	fprintf(stdout,"=========================================================\n");
	fprintf(stdout,"Trajectory files:\n");
	for (int i=0; i<ntraj; i++) {
		fprintf(stdout,"%s\n",dcd_names[i].c_str());
		inptraj.open(dcd_names[i].c_str());
		max_save += inptraj.total_frames();
		inptraj.close();
	}
	fprintf(stdout,"Output root = %s\n", outp_root);
	isTP = new int[max_save];
	read_istp(istp_names,isTP,1,max_save);
	// read trajectories and fill 
	char state = 'i';
	int ltp, leq, tt;
	ltp = 0;
	leq = 0;
	inptraj.open(dcd_names[0].c_str());
	nframe = inptraj.total_frames();
	inptraj.close();
	tot_frame = nframe*dcd_names.size();
	f_ind = -1;
	int tp_start, tp_end;
	int tp_ind = 0;
	int i=0;
	while (i<max_save) {
		if (isTP[i] == 1) {
			tp_start = i++;
			while (i<max_save && (isTP[i++]==1)) ;
			tp_end = i;
			sprintf(buf,"%s_tp%i.dcd",outp_root,tp_ind++);
			write_tp(dcd_names, buf, tp_start, tp_end);
		} else {
			i++;
		}
	}

	return 0;
}


