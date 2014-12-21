/*
 * pTP_crd.cc
 *
 * calculate p(TP|x) for some putative reaction coordinate x, given:
 * (i) a file with a series of equilibrium x-values
 * (ii) a file with a series of values for some reaction coordinate
 * 	used to define equilibrium states & transition paths
 * Needless to say these must be of the same length!!
 *
 * This version assumes state A for values of the equilibrium reaction
 * coordinate less than the -a argument and state B for values greater than -b
 * argument 
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
	pTP_crd -h \n\
		for help \n\
	pTP_crd -a a_max -b b_min -L bin_lo -H bin_hi -n nbin -e eq.dat \n\
		-t tp.dat -o outp.dat\n\
	where:\n\
		-a a_max is the maximum value of the eq. coord defining state A;\n\
		-b b_min is the minimum value of the eq. coord defining state B;\n\
		-L bin_lo is the lower limit for binning the TP coord;\n\
		-H bin_hi is the upper limit for binning the TP coord;\n\
		-n nbin is the number of bins for binning the TP coord;\n\
		-e eq.dat is the file containing values of the equilibrium\n\
			reaction coordinate\n\
		-t eq.dat is the file containing values of the transition path\n\
			reaction coordinate\n\
		-o outp.dat is the file for writing output\n\n";

double next_float(FILE *f)
{
	const int BUF_LEN = 1024;
	char *happyhappy, *tok, buf[BUF_LEN];
	happyhappy = fgets(buf,BUF_LEN,f);
	tok = strtok(buf," \t");
	return atof(tok);
}

//void bin_data(vector<double> &bins, double bin_lo, double delta, double x)
void bin_data(vector<int> &bins, double bin_lo, double delta, double x)
{
	int bin, nbin;
	nbin = bins.size();
	bin = int(floor((x-bin_lo)/delta));
	if (bin >= 0 && bin<nbin)
		bins[bin]++;
	return;
}

//void add_data(vector<double> &a, vector<double> &b)
void add_data(vector<int> &a, vector<int> &b)
{
	/* add a to b */
	for (int i=0; i<a.size(); i++) {
		b[i] += a[i];
	}
	return;
}

int main(int argc, char *argv[]) 
{
	int c, leq, ltp, ntp, tt;
	char state;
	double bin_lo, bin_hi, delta;
	const int BUF_LEN = 1024;
	char eq_crd[BUF_LEN], tp_crd[BUF_LEN], outp_name[BUF_LEN];
	char tplen_name[BUF_LEN];
	char tok[BUF_LEN], buf[BUF_LEN];
	int ndat, nbin;
	int nblock;
	//vector<double> pxTP, pTPx, peqx, ptmp;
	vector<int> pxTP, pTPx, peqx, ptmp;
	vector<int> ncross_eq,ncross_fwd, ncross_rev, ncross_tmp;
	FILE *eq_rxf, *tp_rxf, *outp, *tplen_out;
	double eqx, tpx, prev_tpx, amax, bmin;
	double p_eq_x, p_TP_x, p_x_TP, p_TP;
	// parse cmd line options
	bin_lo = -50000.0;
	bin_hi = -50000.0;
	nbin = 20;
	amax = 0.5;
	bmin = 0.8;
	nblock = 1;

	while (1) {
		c=getopt(argc,argv,"ha:b:L:H:n:e:t:o:T:");
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
			case 'L':
				bin_lo = atof(optarg);
				break;
			case 'H':
				bin_hi = atof(optarg);
				break;
			case 'n':
				nbin = atoi(optarg);
				break;
			case 'e':
				strcpy(eq_crd,optarg);
				break;
			case 't':
				strcpy(tp_crd,optarg);
				break;
			case 'o':
				strcpy(outp_name,optarg);
				break;
			case 'T':
				strcpy(tplen_name,optarg);
				break;
			case 'B':
				nblock = atoi(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	pTPx.resize(nbin);
	pxTP.resize(nbin);
	peqx.resize(nbin);
	ptmp.resize(nbin);
	ncross_eq.resize(nbin);
	ncross_fwd.resize(nbin);
	ncross_rev.resize(nbin);
	ncross_tmp.resize(nbin);
	delta = (bin_hi-bin_lo)/float(nbin);
	for (int i=0; i<nbin; i++) {
		//pTPx[i] = 0.0;
		//pxTP[i] = 0.0;
		//peqx[i] = 0.0;
		//ptmp[i] = 0.0;
		pTPx[i] = 0;
		pxTP[i] = 0;
		peqx[i] = 0;
		ptmp[i] = 0;
		ncross_eq[i] = 0;
		ncross_fwd[i] = 0;
		ncross_rev[i] = 0;
		ncross_tmp[i] = 0;
	}
	eq_rxf = fopen(eq_crd,"r");
	if (eq_rxf == NULL) {
		fprintf(stderr,"Could not read equilibrium coord file: %s\n",eq_crd);
		fprintf(stderr,"%s\n",Usage);
		exit(1);
	}
	tp_rxf = fopen(tp_crd,"r");
	if (tp_rxf == NULL) {
		fprintf(stderr,"Could not read transition coord file: %s\n",tp_crd);
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

	tplen_out = fopen(tplen_name,"w");
	state = 'i';
	leq = 0;
	ltp = 0;
	ntp = 0;
	tpx = -1; // needed?
	for (int i=0; i<ndat; i++) {
		prev_tpx = tpx;
		eqx = next_float(eq_rxf);
		tpx = next_float(tp_rxf);
		//fprintf(stdout,"%12.6f %12.6f\n",eqx,tpx);
		if (state == 'i') {
			if (eqx < amax) { // in A
				state = 'A';
				bin_data(peqx, bin_lo, delta, tpx);
				leq += 1;
			} else if (eqx > bmin) { // in B
				state = 'B';
				bin_data(peqx, bin_lo, delta, tpx);
				leq += 1;
			} else { // not in A or B yet
				continue;		
			}
		} else if (state == 'A' && eqx > amax) { // possible tp
			tt = 1;
			for (int p=0; p<nbin; p++) 
				ptmp[p] = 0.;
			bin_data(ptmp,bin_lo,delta,tpx);
			state = 'a';
		} else if (state == 'B' && eqx < bmin) {
			tt = 1;
			for (int p=0; p<nbin; p++) 
				ptmp[p] = 0.;
			bin_data(ptmp,bin_lo,delta,tpx);
			state = 'b';
		} else if (state == 'a') {
			//fprintf(stdout,"%8.3f\n",eqx);
			if (eqx < amax) {
				bin_data(peqx,bin_lo,delta,tpx);
				add_data(ptmp,peqx);
				leq += tt+1;
				state = 'A';
			} else if (eqx > bmin) {
				bin_data(peqx,bin_lo,delta,tpx);
				//for (int k=0; k<ptmp.size(); k++) {
				//	fprintf(stdout,"%5.3f ",ptmp[k]);
				//}
				//fprintf(stdout,"\n");
				add_data(ptmp,pxTP);
				add_data(ptmp,peqx);
				ltp += tt;
				fprintf(tplen_out,"%i\n", tt);
				leq += tt+1;
				ntp += 1;
				state = 'B';
			} else {
				bin_data(ptmp,bin_lo,delta,tpx);
				tt++;
			}
		} else if (state == 'b') {
			if (eqx > bmin) {
				bin_data(peqx,bin_lo,delta,tpx);
				add_data(ptmp,peqx);
				leq += tt+1;
				state = 'B';
			} else if (eqx < amax) {
				bin_data(peqx,bin_lo,delta,tpx);
				add_data(ptmp,pxTP);
				add_data(ptmp,peqx);
				ltp += tt;
				fprintf(tplen_out,"%i\n", tt);
				leq += tt+1;
				ntp += 1;
				state = 'A';
			} else {
				bin_data(ptmp,bin_lo,delta,tpx);
				tt++;
			}
		} else {
			bin_data(peqx,bin_lo,delta,tpx);
			leq++;
		}
	}
	fclose(eq_rxf);
	fclose(tp_rxf);
	//sprintf(buf,"%s_pTPx.dat",outp_name);
	outp = fopen(outp_name,"w");
	if (outp == NULL) {
		fprintf(stderr,"Could not open file %s for output\n", outp);
		exit(1);
	}
	p_TP = float(ltp)/float(leq);
	fprintf(outp,"# p(TP|x) analysis\n");
	fprintf(outp,"# File with equilibrium coordinate: %s\n", eq_crd);
	fprintf(outp,"# File with transition path coordinate: %s\n", tp_crd);
	fprintf(outp,"# Total length of equilibrium trajectory = %i\n", leq);
	fprintf(outp,"# Total transition path length = %i; p(TP) = %12.5e\n", ltp, p_TP);
	fprintf(outp,"# Number of transition paths = %i\n", ntp);
	fprintf(outp,"#\n");
	fprintf(outp,"# %10s %12s %12s %12s\n", "x","Peq(x)","P(x|TP)","P(TP|x)");
	for (int i=0; i<nbin; i++) {
		double xx = bin_lo + delta*(float(i)+0.5);
		p_eq_x = float(peqx[i])/(float(leq)*delta);
		p_x_TP = float(pxTP[i])/(float(ltp)*delta);
		if (peqx[i] == 0) {
			p_TP_x = 0.0;
		} else {
			p_TP_x = p_x_TP*p_TP/p_eq_x;
		}
		fprintf(outp,"%12.6f %12.6f %12.6f %12.6f\n", xx, p_eq_x, p_x_TP, p_TP_x);
	}
	fclose(outp);
	fclose(tplen_out);
	return 0;
}
