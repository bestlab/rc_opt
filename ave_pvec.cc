/*
 * calculate average and standard deviation for weight matrices
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <unistd.h>
#include "config.h"

using namespace std;

void print_usage(FILE *outp)
{
	fprintf(outp,"\n\nUsage:\n");
	fprintf(outp,"\tave_pvec -g -o ave_out weight_traj\n");
	fprintf(outp,"\twhere:\n");
	fprintf(outp,"\t\t-o outp is the file to write average/sdev w,r to [default: junk.dat]\n");
	fprintf(outp,"\t\tweight_traj is the binary trajectory of weight matrices\n");
	fprintf(outp,"\t\t-g (optional) specifies a gcmc trajectory\n");
	fprintf(outp,"\n\n");
}

int main(int argc, char **argv)
{
	int change_cutoffs,noswap;
	int c, frame, nframe, dim, natom;
	string oname, traj_name;
	oname = "junk.dat";
	proj_type *w, *w_sum, *w2_sum, *r, *r_sum, *r2_sum;
	int *r_cnt;
	double wtmp, mean_w, sdev_w, mean_r, sdev_r;
	int *a;
	bool gcmc = false;
	FILE *traj, *outp;
	int active_contacts;
	frame = 1;
	while (1) {
		c=getopt(argc,argv,"hgo:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				print_usage(stdout);
				exit(0);
				break;
			case 'o':
				oname = optarg;
				break;
			case 'g':
				gcmc = true;
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
	fprintf(stdout,"nframe =  %i\n", nframe);
	fread(&dim,sizeof(int),1,traj);
	fprintf(stdout,"dim =  %i\n", dim);
	fread(&active_contacts,sizeof(int),1,traj);
	fprintf(stdout,"active contacts =  %i\n", active_contacts);
	if (frame > nframe) {
		fprintf(stderr,"\nFrame number (%i) is greater than number of frames (%i)\n",
				frame, nframe);
		print_usage(stderr);
		exit(1);
	}
	natom = int(floor(1.0+sqrt(1.0+8.0*float(dim))))/2;
	fprintf(stdout,"Number of frames = %i\n",nframe);
	fprintf(stdout,"Contact vector dimension = %i\n",dim);
	fprintf(stdout,"Number of active contacts = %i\n",active_contacts);
	fprintf(stdout,"Number of atoms = %i\n",natom);
	fprintf(stdout,"Averaging weight matrices...\n");
	a = new int[active_contacts];
	w = new proj_type[dim];
	r = new proj_type[dim];
	r_sum = new proj_type[dim];
	r2_sum = new proj_type[dim];
	w_sum = new proj_type[dim];
	w2_sum = new proj_type[dim];
	r_cnt = new int[dim];
	fread(&change_cutoffs,sizeof(int),1,traj);
	fprintf(stdout,"change_cutoffs = %i\n", change_cutoffs);
	if (!change_cutoffs) {
		fprintf(stdout,"reading cutoff list\n");
		fread(r,sizeof(proj_type),active_contacts,traj);
	}
	fprintf(stdout,"GOT HERE 0\n");
	fread(&noswap,sizeof(int),1,traj);
	fprintf(stdout,"noswap = %i\n", noswap);
	if (noswap) { 
		fprintf(stdout,"reading active contact list\n");
		fread(a,sizeof(int),active_contacts,traj);
	}
	fprintf(stdout,"GOT HERE 1\n");

	for (int j=0; j<dim; j++) {
		r_sum[j] = 0.0;
		r2_sum[j] = 0.0;
		w_sum[j] = 0.0;
		w2_sum[j] = 0.0;
		r_cnt[j] = 0;
	}
	fprintf(stdout,"GOT HERE 2\n");
	for (int i=0; i<nframe; i++) {
		fread(w,sizeof(proj_type),active_contacts,traj);
		fprintf(stdout,"GOT HERE 2.0\n");
		if (!noswap) 
			fread(a,sizeof(int),active_contacts,traj);
		fprintf(stdout,"GOT HERE 2.1\n");
		if (change_cutoffs)
			fread(r,sizeof(proj_type),active_contacts,traj);
		fprintf(stdout,"GOT HERE 2.2\n");
		for (int j=0; j<active_contacts; j++) {
			//fprintf(stdout,"j=%i; a[j]=%i\n",j,a[j]);
			wtmp = w[j];
			w_sum[a[j]] += wtmp;
			w2_sum[a[j]] += wtmp*wtmp;
			if (!noswap) {
				wtmp = r[j];
				r_sum[a[j]] += wtmp;
				r2_sum[a[j]] += wtmp*wtmp;
				r_cnt[a[j]] += 1;
			}
		}
		fprintf(stdout,"GOT HERE 2.3\n");
	}
	fprintf(stdout,"GOT HERE 3\n");
	outp = fopen(oname.c_str(),"w");
	fprintf(outp,"#%5s %5s %12s %12s %12s %12s\n","I","J","MEAN_W","SDEV_W",
			"MEAN_R","SDEV_R");
	for (int i=1; i<natom; i++) {
		for (int j=0; j<i; j++) {
			int ind = i*(i-1)/2+j;
			mean_w = w_sum[ind] / double(nframe);
			sdev_w = sqrt(w2_sum[ind]/double(nframe)-mean_w*mean_w);
			if (r_cnt[ind]==0) {
				mean_r = 0.0;
				sdev_r = 0.0;
			} else {
				mean_r = r_sum[ind] / double(r_cnt[ind]);
				if (r_cnt[ind]==1) {
					sdev_r = 0.0;
				} else {
					double tmpf = r2_sum[ind]/double(r_cnt[ind])-mean_r*mean_r;
					if (tmpf>0.0) {
						sdev_r = sqrt(tmpf);
					} else {
						sdev_r = 0.0;
					}
				}
			}
			fprintf(outp," %5i %5i %12.5e %12.5e %12.5e %12.5e\n",i+1,j+1,mean_w,sdev_w,mean_r,sdev_r);
		}
	}
	fclose(traj);
	fclose(outp);
	delete [] a;
	delete [] w;
	delete [] r;
	delete [] w_sum;
	delete [] w2_sum;
	delete [] r_sum;
	delete [] r2_sum;
	delete [] r_cnt;
}

