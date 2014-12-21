/*
 * extract a particular weight matrix from a trajectory
 */

#include <cstdlib>
#include <cstdio>
#include <string>
#include <unistd.h>
#include "config.h"

using namespace std;

void print_usage(FILE *outp)
{
	fprintf(outp,"\n\nUsage:\n");
	fprintf(outp,"\txtr_pvec -f frame -o outp_file weight_traj\n");
	fprintf(outp,"\twhere:\n");
	fprintf(outp,"\t\t-f frame is the frame number (1-nframe) to extract [default: 1]\n");
	fprintf(outp,"\t\t-o outp_file is the name of the file to write the rate matrix to\n");
	fprintf(outp,"\t\t\t[default: junk.dat]\n");
	fprintf(outp,"\t\tweight_traj is the binary trajectory of weight matrices\n");
	fprintf(outp,"\n\n");
}

int main(int argc, char **argv)
{
	int c, frame, nframe, dim;
	string outp_name, traj_name;
	outp_name = "junk.dat";
	proj_type *w;
	proj_type *r;
	int *a;
	FILE *traj, *outp;
	int active_contacts;
	int change_cutoffs, noswap;
	frame = 1;
	while (1) {
		c=getopt(argc,argv,"hf:o:c");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				print_usage(stdout);
				exit(0);
				break;
			case 'f':
				frame = atoi(optarg);
				break;
			case 'o':
				outp_name = optarg;
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				print_usage(stderr);
				exit(1);
		}
	}
	if (argc == optind) {
		fprintf(stderr,"\nERROR: No trajectory specified\n");
		print_usage(stderr);
		exit(1);
	} else {
		traj_name = argv[optind];
	}
	traj = fopen(traj_name.c_str(),"r");
	fread(&nframe,sizeof(int),1,traj);
	fread(&dim,sizeof(int),1,traj);
	fread(&active_contacts,sizeof(int),1,traj);
	if (frame > nframe) {
		fprintf(stderr,"\nFrame number (%i) is greater than number of frames (%i)\n",
				frame, nframe);
		print_usage(stderr);
		exit(1);
	}
	fprintf(stdout,"Number of frames = %i\n",nframe);
	fprintf(stdout,"Frame to extract = %i\n",frame);
	fprintf(stdout,"Contact vector dimension = %i\n",dim);
	fprintf(stdout,"Number of active contacts = %i\n",active_contacts);
	w = new proj_type[dim];
	a = new int[dim];
	r = new proj_type[dim];
	fread(&change_cutoffs,sizeof(int),1,traj);
	if (change_cutoffs) {
		fread(r,sizeof(proj_type),active_contacts,traj);
	}
	fread(&noswap,sizeof(int),1,traj);
	if (noswap) { 
		fread(a,sizeof(int),active_contacts,traj);
	}
	for (int i=0; i<frame; i++) {
		fread(w,sizeof(proj_type),active_contacts,traj);
		if (!noswap) 
			fread(a,sizeof(int),active_contacts,traj);
		if (change_cutoffs)
			fread(r,sizeof(proj_type),active_contacts,traj);
	}
	fclose(traj);
	outp = fopen(outp_name.c_str(),"w");
	for (int i=0; i<active_contacts; i++) {
		fprintf(outp,"%6i %12.5e %12.5e\n", a[i],r[i],w[i]);
	}
	fclose(outp);
	delete [] w;
	delete [] a;
	delete [] r;
}

