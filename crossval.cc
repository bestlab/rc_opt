/*
 * calculate p(TP|r) given a projection and assignment of states
 */



#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "TrajFile.h"
#include "config.h"

using namespace std;

const char *Usage = "\n\
	Usage: \n\
		crossval -o ptp.dat -L proj_lo -H proj_hi -n proj_nbin \n\
			trj1_proj.dat istp1.dcd ... trjN_proj.dat istpN.dcd\n\n";

template<class T>
inline T min(T a, T b) {
	return a < b ? a : b;
}

template<class T>
inline T sign(T x) {
	return x < 0.0 ? -1.0 : 1.0;
}

template<class T>
inline T contact(T dr, T cut) 
{
	return 1.0/(1.0+exp(5.0*(dr-cut)));
}

template<class T>
void read_files(vector<string> &file_names, vector<T> &data)
{
	ifstream inp;
	int ndat,idx,ndatf;
	T tmpi;
	int nfile = file_names.size();
	vector<int> file_len;
	ndat = 0;
	file_len.resize(file_names.size());
	for (int f=0; f<nfile; f++) {
		ndatf = 0;
		inp.open(file_names[f].c_str());
		if (!inp.good()) {
			fprintf(stderr,"Could not open file %s\n",file_names[f].c_str());
			exit(1);
		}
		while (!inp.eof()) {
			inp >> tmpi;
			ndatf++;
		}
		ndatf--;
		file_len[f] = ndatf;
		ndat += ndatf;
		inp.close();
		inp.clear();
	}
	data.resize(ndat);
	idx=0;
	for (int f=0; f<nfile; f++) {
		inp.open(file_names[f].c_str());
		for (int l=0; l<file_len[f]; l++) {
			inp >> tmpi;
			data[idx++] = tmpi;
		}
		inp.close();
		inp.clear();
	}
}

void write_ptp(const char *fname, double bin_lo, double bin_hi, int nbin, 
		vector<vector<proj_type> > &p_r, vector<vector<proj_type> > &p_rTP, 
		vector<vector<proj_type> > &p_TPr,
		vector<proj_type> &p_r_glob, vector<proj_type> &p_rTP_glob, vector<proj_type> &p_TPr_glob,
		int leq, int ltp, int ntp, double pTP_mu, double pTP_sig,
		double rate_mean, double rate_sdev)
{
	double dr, pri; 
	FILE * ptp_out = fopen(fname,"w"); 
	if (ptp_out == NULL) {
		fprintf(stderr,"Could not open file %s for output\n",fname);
		exit(1);
	}
	int nblock = p_r.size();
	double dnblock = double(nblock);
	fprintf(ptp_out,"# Total number of points in equilibrium histogram = %-8i\n",leq);
	fprintf(ptp_out,"# Total number of points in transition path histogram = %-8i\n",ltp);
	fprintf(ptp_out,"# Global p(TP) = %-12.8f\n",float(ltp)/float(leq));
	fprintf(ptp_out,"# Block average p(TP) = %-12.8f\n",pTP_mu);
	fprintf(ptp_out,"# Block p(TP) st.dev. = %-12.8f\n",pTP_sig);
	fprintf(ptp_out,"# Block p(TP) st.err. = %-12.8f\n",pTP_sig/sqrt(float(nblock)));
	fprintf(ptp_out,"# Number of transition paths = %-8i\n",ntp);
	fprintf(ptp_out,"# Average transition path length = %-8i (units = save interval)\n",
			float(ltp)/float(ntp));
	fprintf(ptp_out,"# Global rate =2/(tf+tu)=p(TP)/<t_TP>= %-12.6e (units = 1/(save interval))\n",
			float(ntp)/float(leq));
	fprintf(ptp_out,"# Block average rate = %-12.6e (units = 1/(save interval))\n",rate_mean);
	fprintf(ptp_out,"# Rate std.err = %-12.6e (units = 1/(save interval))\n",rate_sdev);
	dr = (bin_hi-bin_lo)/double(nbin);
	for (int i=0; i<nbin; i++) {
		pri = bin_lo+(double(i)+0.5)*dr;
		double pr_sum, pr_sum2,prtp_sum,prtp_sum2,ptpr_sum,ptpr_sum2;
		double pr,prtp,ptpr;
		pr_sum= pr_sum2=prtp_sum=prtp_sum2=ptpr_sum=ptpr_sum2 = 0.0;
		for (int b=0; b<nblock; b++) {
			pr = p_r[b][i];
			pr_sum += pr; pr_sum2 += pr*pr;
			prtp = p_rTP[b][i];
			prtp_sum += prtp; prtp_sum2 += prtp*prtp;
			ptpr = p_TPr[b][i];
			ptpr_sum += ptpr; ptpr_sum2 += ptpr*ptpr;
		}
		pr_sum /= dnblock;
		pr_sum2 = sqrt(pr_sum2/dnblock-pr_sum*pr_sum)/sqrt(dnblock);
		prtp_sum /= dnblock;
		prtp_sum2 = sqrt(prtp_sum2/dnblock-prtp_sum*prtp_sum)/sqrt(dnblock);
		ptpr_sum /= dnblock;
		ptpr_sum2 = sqrt(ptpr_sum2/dnblock-ptpr_sum*ptpr_sum)/sqrt(dnblock);
		fprintf(ptp_out,"%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",pri,p_r_glob[i],pr_sum2,
				p_rTP_glob[i],prtp_sum2,p_TPr_glob[i],ptpr_sum2);
	}
	fclose(ptp_out);
}

void pTPr_calc(vector<int> &isTP, vector<proj_type> &proj, double bin_lo, double bin_hi, int nbin,
		vector<proj_type> &p_r, vector<proj_type> &p_rTP, vector<proj_type> &p_TPr,
		int &leq, int &ltp, int &ntp, int lo, int hi)
{
	// only use data points with idx%base=mod
	int nsave;
	int bin;
	double r, dr;
	nsave = isTP.size();
	p_r.resize(nbin);
	p_rTP.resize(nbin);
	p_TPr.resize(nbin);
	dr = (bin_hi-bin_lo)/double(nbin);
	for (int i=0; i<nbin; i++) {
		p_r[i] = 0.0;
		p_rTP[i] = 0.0;
		p_TPr[i] = 0.0;
	}
	leq = 0; ltp = 0;
	ntp = 0;
	bool ontp = false;
	for (int s=lo; s<hi; s++) {
		r = proj[s];
		bin = int(floor((r-bin_lo)/dr));
		if (bin>=0 && bin<nbin) {
			p_r[bin] += 1.0;
			leq++;
			if (isTP[s]==1) {
				if (!ontp) {
					ontp = true;
					ntp++;
				}
				p_rTP[bin] += 1.0;
				ltp++;
			} else {
				ontp = false;
			}
		} else {
			fprintf(stdout,"r = %12.5f, bin = %i\n", r, bin);
		}
	}
	for (int i=0; i<nbin; i++) {
		if (p_r[i] > 10.0) {
			p_TPr[i] = p_rTP[i]/p_r[i];
		} else {
			p_TPr[i] = 0.0;
		}
		p_r[i] /= double(leq);
		p_rTP[i] /= double(ltp);
	}
	return;
}


int main(int argc, char **argv)
{
	vector<proj_type> p_r_glob, p_rTP_glob, p_TPr_glob;
	vector<vector<proj_type> > p_r, p_rTP, p_TPr;
	vector<string> proj_names, istp_names;
	vector<int> isTP;
	vector<proj_type> proj;
	string oname;
	int nproj, nbin, c, nblock;
	int leq,ltp,ntp,leq_glob,ltp_glob,ntp_glob,ndat;
	double ptp, pTP_sum, pTP_sum2;
	double bin_lo, bin_hi;
	bin_lo =0.0;
	bin_hi =1.0;
	nbin = 100;
	nblock = 10;
	/*
	double ptot, ptot_trial, psdev, psdev_trial, pmean, pmean_trial, crit, crit_trial, best_crit;
	string pdb_name, oroot, qfile,pvec_traj_name, contact_file;
	int skip, c, nbin, ntraj, nsave, n_proj_save, max_save, tmpi, tmpj;
	proj_type *rc, *proj, *nu_proj, *pr1, *pr2, *pv_tmp, *rcut_tmp;
	int *ai_tmp;
	int *isTP;
	string  pvec_file = "NONE";
	string  istp_file = "NONE";
	vector<int> ii,jj;
	const int max_pl = 100000;
	int pl;
	double alpha = 0.2;
	double sigma_rc = 2.0;
	int nmc = 1000;
	int TRUE = 1;
	int FALSE = 0;
	double lambda;
	double mean_a,var_a,mean_b,var_b;
	double T, rc_i, rc_i_nu;
	double default_cutoff = 12.0;		// Angstroms
	int move_attempts[nmove], move_successes[nmove];
	int random_seed = 100;
	char decider = 'q';		// = 'q' (integrated pTP/sdev)
					// = 'i' (integrated pTP)
					// = 'm' (maximum)
					// = 'g' (maximum from gaussian fit)
	// active moves = 3 for all contacts and 4 if # active contacts < total contacts
	FILE *proj_out, *pvec_out, *pvec_traj; //, *ptp_out;
	string proj_name, pvec_name, ptp_name;
	*/

	if (argc < 3) {
		fprintf(stderr,Usage);
		exit(1);
	}
	while (1) {
		c=getopt(argc,argv,"ho:L:H:n:b:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'o':
				oname = optarg;
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
			case 'b':
				nblock = atoi(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	nproj = (argc-optind)/2;	// half dcd, half istp
	proj_names.resize(nproj);
	istp_names.resize(nproj);
	for (int i=0; i<nproj; i++) {
		proj_names[i] = argv[optind+2*i];
		istp_names[i] = argv[optind+2*i+1];
	}
	read_files(proj_names,proj);
	read_files(istp_names,isTP);
	if (proj.size() != isTP.size()) {
		fprintf(stderr,"projection and istp data of different sizes!\n");
		exit(1);
	}
	ndat = proj.size();
	int blocksize=ndat/nblock;
	p_r.resize(nblock);
	p_rTP.resize(nblock);
	p_TPr.resize(nblock);
	fprintf(stdout,"Number of data points = %i\n",ndat);
	fprintf(stdout,"Number of blocks = %i\n",nblock);
	fprintf(stdout,"Block size = %i\n",blocksize);
	double rate, rate_sum, rate_sum2;
	rate_sum = rate_sum2 = 0.0;
	pTP_sum = pTP_sum2 = 0.0;
	for (int b=0; b<nblock; b++) {
		fprintf(stdout,"Got here block %i\n",b);
		pTPr_calc(isTP, proj, bin_lo, bin_hi, nbin, p_r[b], p_rTP[b], 
				p_TPr[b],leq,ltp,ntp,b*blocksize,(b+1)*blocksize); // block
		ptp = double(ltp)/double(leq);
		pTP_sum += ptp;
		pTP_sum2 += ptp*ptp;
		rate = double(ntp)/double(leq);
		rate_sum += rate;
		rate_sum2 += rate*rate;
	}
	fprintf(stdout,"Got here A\n");
	double pTP_mean = pTP_sum/double(nblock);
	double pTP_sdev = sqrt(pTP_sum2/double(nblock)-pTP_mean*pTP_mean);
	rate_sum = rate/double(nblock);
	rate_sum2 = sqrt(rate_sum2/double(nblock)-rate_sum*rate_sum);
	pTPr_calc(isTP, proj, bin_lo, bin_hi, nbin, p_r_glob, p_rTP_glob, p_TPr_glob,
			leq,ltp,ntp,0,ndat); // global
	fprintf(stdout,"Got here B\n");
	write_ptp(oname.c_str(), bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr,
			p_r_glob, p_rTP_glob, p_TPr_glob,leq,ltp,ntp,
			pTP_mean, pTP_sdev, rate_sum,rate_sum2);
	fprintf(stdout,"Got here C\n");

	return 0;
}

