/*
 * pTP_crd.cc
 *
 * assign points to transition paths
 */

#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;

const char * Usage = "\n\n\
	Usage:\n\
\n\
	assign_tp -h \n\
		for help \n\
	assign_tp -a a_max -b b_min -r r.dat -o outp.dat [-A amax2 -B bmin2 -R r2.dat]\n\
\n\
	Will write a sequence of numbers, 1 per line\n\
		0 for state A\n\
		1 for TP\n\
		2 for state B\n\
	where:\n\
		-a a_max is the maximum value of the coord defining state A;\n\
		-b b_min is the minimum value of the coord defining state B;\n\
		-r r.dat is the file containing values of the \n\
			reaction coordinate\n\
		-A amax2 is the maximum value of the 2nd coord defining state A;\n\
		-B bmin2 is the minimum value of the 2nd coord defining state B;\n\
		-R r2.dat is the file containing values of the \n\
			(optional) second reaction coordinate\n\
		-o outp.dat is the file for writing output\n\n";

double next_float(FILE *f)
{
	const int BUF_LEN = 1024;
	char *tok, buf[BUF_LEN];
	fgets(buf,BUF_LEN,f);
	tok = strtok(buf," \t");
	return atof(tok);
}

void bin_data(vector<double> &bins, double bin_lo, double delta, double x)
{
	int bin, nbin;
	nbin = bins.size();
	bin = int(floor((x-bin_lo)/delta));
	if (bin >= 0 && bin<nbin)
		bins[bin]++;
	return;
}

void add_data(vector<double> &a, vector<double> &b)
{
	/* add a to b */
	for (int i=0; i<a.size(); i++) {
		b[i] += a[i];
	}
	return;
}

int main(int argc, char *argv[]) 
{
	int c, tt;
	char state;
	const int BUF_LEN = 1024;
	char eq_crd[BUF_LEN], eq_crd2[BUF_LEN], outp_name[BUF_LEN];
	char tok[BUF_LEN], buf[BUF_LEN];
	int ndat, nbin;
	vector<double> pxTP, pTPx, peqx, ptmp;
	FILE *eq_rxf, *eq_rxf2, *outp;
	double eqx, eqx2, tpx, amax, bmin, amax2, bmin2;
	// parse cmd line options
	amax = 0.5;
	bmin = 0.8;
	amax2 = 0.5;
	bmin2 = 0.8;

	while (1) {
		c=getopt(argc,argv,"ha:b:r:A:B:R:o:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'a':
				amax = atof(optarg);
				break;
			case 'b':
				bmin = atof(optarg);
				break;
			case 'r':
				strcpy(eq_crd,optarg);
				break;
			case 'A':
				amax2 = atof(optarg);
				break;
			case 'B':
				bmin2 = atof(optarg);
				break;
			case 'R':
				strcpy(eq_crd2,optarg);
				break;
			case 'o':
				strcpy(outp_name,optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	fprintf(stdout,"state A max = %8.3f\n",amax);
	fprintf(stdout,"state B min = %8.3f\n",bmin);
	eq_rxf = fopen(eq_crd,"r");
	fprintf(stdout,"[rc2] state A max = %8.3f\n",amax2);
	fprintf(stdout,"[rc2] state B min = %8.3f\n",bmin2);
	eq_rxf2 = fopen(eq_crd2,"r");
	if (eq_rxf == NULL) {
		fprintf(stderr,"Could not read rxn coord file: %s\n",eq_crd);
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	if (eq_rxf2 == NULL) {
		fprintf(stderr,"Could not read rxn coord file: %s\n",eq_crd2);
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	// count equilibrium data
	fgets(buf,BUF_LEN,eq_rxf);
	ndat = 0;
	while(feof(eq_rxf)==0) {
		ndat++;
		fgets(buf,BUF_LEN,eq_rxf);
	}
	fclose(eq_rxf);
	eq_rxf = fopen(eq_crd,"r");

	state = 'i';
	outp = fopen(outp_name,"w");
	if (outp == NULL) {
		fprintf(stderr,"Could not open file %s for output\n", outp);
		exit(1);
	}
	fprintf(stdout,"ndat = %i\n",ndat);
	tt = 0;
	for (int i=0; i<ndat; i++) {
		eqx = next_float(eq_rxf);
		eqx2 = next_float(eq_rxf2);
		//fprintf(stdout,"%12.6f %12.6f\n",eqx,eqx2);
		//fprintf(outp,"%1i\t",0);
		//continue;
		if (state == 'i') {
			if (eqx < amax && eqx2 < amax2) { // in A
				state = 'A';
				tt++;
				for (int j=0; j<tt; j++) {
					fprintf(outp,"%1i\n",0); 	// in A
				}
			} else if (eqx > bmin && eqx2 > bmin2) { // in B
				state = 'B';
				tt++;
				for (int j=0; j<tt; j++) {
					fprintf(outp,"%1i\n",2); 	// in B
				}
			} else { // not in A or B yet
				tt++;
		//		fprintf(stdout,"GOT HERE!\n");
				continue;		
			}
		} else if (state == 'A' && (eqx>amax || eqx2>amax2)) { // possible tp
			tt = 1;
			state = 'a';
		} else if (state == 'B' && (eqx<bmin || eqx2<bmin2)) {
			tt = 1;
			state = 'b';
		} else if (state == 'a') {
			if (eqx < amax && eqx2 < amax2) {
				tt++;
				for (int j=0; j<tt; j++) {
					fprintf(outp,"%1i\n",0); 	// in A
				}
				state = 'A';
			} else if (eqx > bmin && eqx2 > bmin2) {
				for (int j=0; j<tt; j++) {
					fprintf(outp,"%1i\n",1); 	// is tp
				}
				fprintf(outp,"%1i\n",2); 	// now in B
				state = 'B';
			} else {
				tt++;
			}
		} else if (state == 'b') {
			if (eqx > bmin && eqx2 > bmin2) {
				tt++;
				for (int j=0; j<tt; j++) {
					fprintf(outp,"%1i\n",2); 	//  in B
				}
				state = 'B';
			} else if (eqx < amax && eqx2 < amax2) {
				for (int j=0; j<tt; j++) {
					fprintf(outp,"%1i\n",1); 	// is tp
				}
				fprintf(outp,"%1i\n",0); 	// now in A
				state = 'A';
			} else {
				tt++;
			}
		} else {
			if (state == 'A') {
				fprintf(outp,"%1i\n",0);
			} else {
				fprintf(outp,"%1i\n",2);
			}
		}
	}
	if (state == 'a') {
		for (int j=0; j<tt; j++) {
			fprintf(outp,"%1i\n",0); 	// in A
		}
	} else if (state == 'b') {
		for (int j=0; j<tt; j++) {
			fprintf(outp,"%1i\n",2); 	// in B
		}
	}

	fclose(eq_rxf);
	fclose(eq_rxf2);
	fclose(outp);
	return 0;
}
