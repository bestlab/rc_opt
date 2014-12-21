/*
 * find optimal reaction coordinate r (maximize area of P(TP|r) distro)
 * using monte carlo simulation
 *
 * with Boltzmann sampling at constant temperature
 *
 * this version works only for equilibrium trajectories
 *
 * The MC moves in this version have been adapted to (try to) satisfy
 * detailed balance
 *
 * constant size of reaction coordinate (N) version
 *
 * including possibility of using all contacts
 *
 * and also variable contact distances
 *
 * 3/21/05 GCMC added 
 * 
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
//#include "TrajFile.h"
#include "traj.h"
#include "random_gen.h"
#include "config.h"
#include <sys/time.h>

using namespace std;

const char *Usage = "\n"
"	Usage: \n"
"		rc_opt [-r random_seed] [-d g|i|q|m|a|b|A] \n"
"			-L proj_lo -H proj_hi -n proj_nbin [-s skip] -p xyz.pdb\n"
"			[-c intial_contact_vector ]\n"
"			-q contact_list -m n_mc -T temperature -f outp_freq\n"
"			-o outp_root -N active_contacts\n"
"			-g mu \n"
"			file1.dcd file2.dcd ... fileN.dcd\n"
"	where:\n"
"		-d <opt> determines the optimization criterion as follows:\n"
"			<opt>=g: Optimize maximimum of Gaussian\n"
"			<opt>=i: Optimize integral of p(TP|r)\n"
"			<opt>=q: Optimize integral of p(TP|r) / variance of p(TP|r)\n"
"			<opt>=m: Directly optimize maximum of p(TP|r)\n"
"			<opt>=a: As for g but with mean constrained to that of the\n"
"				equilibrium distribution [RECOMMENDED]\n"
"			<opt>=b: As for q but with mean constrained to that of the\n"
"				equilibrium distribution\n"
"			<opt>=A: Optimize (<r>_f-<r>_u)^2/(var(r)_f*var(r)_u)\n"
"			<opt>=h: As for g but with mean constrained to 0.5\n"
"		-L proj_lo -H proj_hi determine the range for binning the \n"
"			projection coordinate\n"
"		-n proj_bin is the number of bins for the projection coordinate\n"
"		-g mu specifies chemical potential mu*active_contacts and GCMC\n"
"		-m nmc is the number of optimization iterations to do\n"
"		-o outp_root is a base name for output files\n"
"		-N active contacts is the number of contacts to use in rxn coord \n"
"			projection\n"
"	Trajectories do not have to be contiguous as they are assumed independent\n\n";

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

double save_pvec(const char *fname, vector<int> &active, vector<proj_type> &pvec, vector<proj_type> &rcut, int n) 
{
		FILE *pvec_out = fopen(fname,"w"); 
		double norm = 0.0;
		for (int i=0; i<n; i++) {
			fprintf(pvec_out,"%6i %12.8e %12.5f\n",active[i],
					pvec[i],rcut[i]);
			norm += pvec[i]*pvec[i];
		}
		fclose(pvec_out);
		return norm;
}

void extend_rc_vec(vector<proj_type> &pvec, vector<proj_type> &rcut_vec, 
		vector<int> &active_ind, int active_contacts, int dim,
		int existing, park_miller &gen, double defcut)
{
	double norm;
	bool used;
	pvec.resize(active_contacts);
	rcut_vec.resize(active_contacts);
	active_ind.resize(active_contacts);
	norm = 0.0;
	for (int i=existing; i<active_contacts; i++) {
		// check this not used already
		used = true;
		while (used) {
			active_ind[i] = int(floor(gen.random()*dim));
			used = false;
			for (int j=0; j<i; j++) {
				if (active_ind[i] == active_ind[j]) {
					used = true;
					break;
				}
			}
		}
		pvec[i] = 0.0;
		rcut_vec[i] = defcut;
	}
	return;
}

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

int read_trajectories(vector<string> &dcd_names, crd_type *X, crd_type *Y, crd_type *Z, 
		int natom, int skip,
		int max_save, int &nsave)
{
	crd_type *xtmp, *ytmp, *ztmp;
	int natom_trj;
	int ntraj = dcd_names.size();
	//DCDITrajFile inptraj;
	BaseITrajFile *inptraj;
	xtmp = new crd_type[natom];
	ytmp = new crd_type[natom];
	ztmp = new crd_type[natom];
	nsave = 0;
	for (int t=0; t<ntraj; t++) {
		string st;
		inptraj = open_trajectory(dcd_names[t].c_str());
		//open_trajectory(st,inptraj);
		natom_trj = inptraj->num_atoms();
		if (natom_trj != natom) {
			fprintf(stderr,"# atoms in trajectory and structure don't match!\n");
			fprintf(stderr,"	trajectory: %i\n",natom_trj);
			fprintf(stderr,"	structure: %i\n",natom);
			exit(1);
		}
		//int nframe = inptraj->total_frames();
		int i = 0;
		while (inptraj->frames_left()) {
		//for (int i=0; i<nframe; i++) {
			if (i%skip != 0) {
				inptraj->read_frame(xtmp,ytmp,ztmp,natom);
			} else {
				inptraj->read_frame(&X[nsave*natom],&Y[nsave*natom],&Z[nsave*natom],natom);
				nsave++;
			}
			if (nsave == max_save) {
				fprintf(stderr,"nsave > maxsave!!! (How does that happen?)\n");
				delete [] xtmp;
				delete [] ytmp;
				delete [] ztmp;
				return 1;
			}
		}
		inptraj->close();
	}
	delete [] xtmp;
	delete [] ytmp;
	delete [] ztmp;
	return 0;
}

proj_type norm(vector<proj_type> &vec)
{
	int s = vec.size();
	double norm = 0.0;
	for (int i=0; i<s; i++) {
		norm += vec[i]*vec[i];
	}
	return sqrt(norm);
}

void choose2(park_miller &gen, int nactive, int &elem_i, int &elem_j)
{
	elem_i = elem_j = int(floor(gen.random()*nactive));
	while (elem_j == elem_i) {
		elem_j = int(floor(gen.random()*nactive));
	}
	return;
}

void generate_random_vec(vector<proj_type> &pvec, vector<proj_type> &rcut_vec, 
		vector<int> &active_ind, int active_contacts, int dim,
		park_miller &gen, double defcut, bool posweights)
{
	double sum;
	bool used;
	pvec.resize(active_contacts);
	rcut_vec.resize(active_contacts);
	active_ind.resize(active_contacts);
	sum = 0.0;
	for (int i=0; i<active_contacts; i++) {
		// check this not used already
		used = true;
		while (used) {
			active_ind[i] = int(floor(gen.random()*dim));
			used = false;
			for (int j=0; j<i; j++) {
				if (active_ind[i] == active_ind[j]) {
					used = true;
					break;
				}
			}
		}
		pvec[i] = gen.random();
		if (posweights) 
			pvec[i] = abs(pvec[i]);
		sum += pvec[i];
		rcut_vec[i] = defcut;
	}
	for (int i=0; i<active_contacts; i++) {
		pvec[i]/=sum; 
	}
	return;
}

void curtime(long &s, long & us)
{
        struct timeval curtime;
        gettimeofday( &curtime, 0 );
        s = curtime.tv_sec;
        us = curtime.tv_usec;
}

void write_ptp(const char *fname, double bin_lo, double bin_hi, int nbin, 
		vector<proj_type> &p_r, vector<proj_type> &p_rTP, vector<proj_type> &p_TPr)
{
	double dr, pri; 
	FILE * ptp_out = fopen(fname,"w"); 
	if (ptp_out == NULL) {
		fprintf(stderr,"Could not open file %s for output\n",fname);
		exit(1);
	}
	dr = (bin_hi-bin_lo)/double(nbin);
	for (int i=0; i<nbin; i++) {
		pri = bin_lo+(double(i)+0.5)*dr;
		fprintf(ptp_out,"%12.5f %12.5f  %12.5f  %12.5f\n",pri,p_r[i],p_rTP[i],p_TPr[i]);
	}
	fclose(ptp_out);
}

proj_type pTP_max(vector<proj_type> &pTPr) 
{
	proj_type tmp = 0.0;
	for (int i=0; i<pTPr.size(); i++) 
		if (pTPr[i]>tmp) 
			tmp = pTPr[i];
	return tmp;
}

proj_type pTP_sum(vector<proj_type> &pTPr) 
{
	proj_type tmp = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp += pTPr[i];
	}
	return tmp;
}

proj_type pTP_quotient(vector<proj_type> &pTPr) 
{
	double tmp, sum, var, mean, sdev;
	// first calc mean, then var for stability
	mean = 0.0;
	sum = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		sum += pTPr[i];
		mean += double(i)*pTPr[i];
	}
	mean /= sum;
	var = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		var += pTPr[i]*tmp*tmp;
	}
	sdev = sqrt(var/sum);
	return sum/sdev;
}

proj_type pTP_quotient_eq(vector<proj_type> &pTPr, vector<proj_type> &peq) 
{
	// calc integral(p(TP|r))/sdev,
	// but calc variance using the EQUILIBRIUM
	// mean
	double tmp, sum, var, mean, sdev;
	mean = 0.0;
	sum = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		sum += peq[i];
		mean += double(i)*peq[i];
	}
	mean /= sum;
	sum = 0.0;
	var = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		sum += pTPr[i];
		var += pTPr[i]*tmp*tmp;
	}
	sdev = sqrt(var/sum);
	return sum/sdev;
}

proj_type pTP_gauss(vector<proj_type> &pTPr) 
{
	double tmp, sum, var, mean, sdev, num, den;
	// first calc mean, then var for stability
	mean = 0.0;
	sum = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		sum += pTPr[i];
		mean += double(i)*pTPr[i];
	}
	mean /= sum;
	var = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		var += pTPr[i]*tmp*tmp;
	}
	//sdev = sqrt(var/sum);
	var /= sum;
	sdev = sqrt(var);
	num = 0.0; den = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		num += pTPr[i]*exp(-tmp*tmp/(2.0*var));
		den += exp(-tmp*tmp/(var));
	}
	return num/den;
}

proj_type pTP_gauss_eq(vector<proj_type> &pTPr, vector<proj_type> &peq) 
{
	// fit gaussian to p(TP|r), but use <r> from the
	// EQUILIBRIUM distribution - hopefully this will
	// shut down any spurious side peaks from atypical
	// transition paths
	double tmp, sum, var, mean, sdev, num, den;
	// Eq. mean
	mean = 0.0;
	sum = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		sum += peq[i];
		mean += double(i)*peq[i];
	}
	mean /= sum;
	sum = 0.0;
	var = 0.0;
	// pTPR var
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		sum += pTPr[i];
		var += pTPr[i]*tmp*tmp;
	}
	var/=sum;
	sdev = sqrt(var);
	num = 0.0; den = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		num += pTPr[i]*exp(-tmp*tmp/(2.0*var));
		den += exp(-tmp*tmp/(var));
	}
	return num/den;
}

proj_type pTP_gauss_half(vector<proj_type> &pTPr, double bin_lo, double bin_hi, int nbin) 
{
	// fit gaussian to p(TP|r), but use <r> = 0.5
	double tmp, sum, var, mean, sdev, num, den;
	double binw;
	int binhalf;
	if (bin_lo >= 0.5 || bin_hi <= 0.5) {
		fprintf(stderr,"With 'h' decider, 0.5 must be within bin range\n");
		exit(1);
	}
	binw = (bin_hi-bin_lo)/double(nbin);
	mean = (0.5-bin_lo)/binw;
	// Eq. mean
	//mean = 0.5;
	//sum = 0.0;
	//for (int i=0; i<pTPr.size(); i++) {
	//	sum += peq[i];
	//	mean += double(i)*peq[i];
	//}
	//mean /= sum;
	sum = 0.0;
	var = 0.0;
	// pTPR var
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		sum += pTPr[i];
		var += pTPr[i]*tmp*tmp;
	}
	var/=sum;
	sdev = sqrt(var);
	num = 0.0; den = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		num += pTPr[i]*exp(-tmp*tmp/(2.0*var));
		den += exp(-tmp*tmp/(var));
	}
	return num/den;
}

proj_type attila_crit(int *isTP, proj_type *proj, int nsave, double &mean_a, double &var_a,
		double &mean_b, double &var_b)
{
	double tmpf, sum_a, sum_a2, sum_b, sum_b2, crit;
	int na, nb;
	sum_a = 0.0;
	sum_a2 = 0.0;
	sum_b = 0.0;
	sum_b2 = 0.0;
	na = 0;
	nb = 0;
	for (int i=0; i<nsave; i++) {
		if (isTP[i] == 0) {
			tmpf = proj[i];
			sum_a += tmpf;
			sum_a2 += tmpf*tmpf;
			na++;
		} else if (isTP[i] == 2) {
			tmpf = proj[i];
			sum_b += tmpf;
			sum_b2 += tmpf*tmpf;
			nb++;
		}
	}
	mean_a = sum_a/double(na);
	mean_b = sum_b/double(nb);
	var_a = sum_a2/double(na) - mean_a*mean_a;
	var_b = sum_b2/double(nb) - mean_b*mean_b;
	crit = mean_a-mean_b;
	crit = crit*crit/(var_a+var_b);
	return crit;
}

void pTP_stats(vector<proj_type> &pTPr, double &mean, double &sdev) 
{
	double tmp, sum, var;
	// first calc mean, then var for stability
	mean = 0.0;
	sum = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		sum += pTPr[i];
		mean += double(i)*pTPr[i];
	}
	mean /= sum;
	var = 0.0;
	for (int i=0; i<pTPr.size(); i++) {
		tmp = double(i) - mean;
		var += pTPr[i]*tmp*tmp;
	}
	sdev = sqrt(var/sum);
	return;
}

proj_type dotp(vector<proj_type> &a, vector<proj_type> &b)
{
	int dim = a.size();
	double p = 0.0;
	if (dim != b.size()) {
		cerr << "dotp: vectors of different sizes!" << endl;
	}
	for (int i=0; i<dim; i++) 
		p += a[i]*b[i];
	return p;
}

inline int index(int i, int j) 
{
	return j>=i ? -1: i*(i-1)/2+j;
}

void get_ij(int ind, int &i, int &j) 
{
	i = int(floor((1.0+sqrt(1.0+8.0*ind))/2.0));
	j = ind - i*(i-1)/2;
	return;
}

void pTPr_calc(int *isTP, proj_type *proj, int nsave, double bin_lo, double bin_hi, int nbin,
		vector<proj_type> &p_r, vector<proj_type> &p_rTP, vector<proj_type> &p_TPr)
{
	int bin, leq, ltp;;
	double r, dr;
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
	for (int s=0; s<nsave; s++) {
		r = proj[s];
		bin = int(floor((r-bin_lo)/dr));
		if (bin>=0 && bin<nbin) {
			p_r[bin] += 1.0;
			leq++;
			if (isTP[s]==1) {
				p_rTP[bin] += 1.0;
				ltp++;
			}
		} else {
			fprintf(stdout,"r = %12.5e, bin = %i\n", r, bin);
			exit(1);
		}
	}
	//fprintf(stdout,"leq = %i, ltp = %i\n", leq, ltp);
	for (int i=0; i<nbin; i++) {
		if (p_r[i] > 10.0) {
			p_TPr[i] = p_rTP[i]/p_r[i];
		} else {
			p_TPr[i] = 0.0;
		}
		p_r[i] /= double(leq);
		p_rTP[i] /= double(ltp);
		//fprintf(stdout,"p_r, p_rTP, p_TPr = %12.5f %12.5f %12.5f\n",
			//	p_r[i], p_rTP[i], p_TPr[i]) 
	}
	return;
}

void proj_update(crd_type *X, crd_type *Y, crd_type *Z, int natom,
		int nsave, proj_type *proj, proj_type *nu_proj, int p, 
		int q, double dp, double dq, double rc_p, double rc_q)
{
	/* only update projections based on exchange of weight lambda between components
	 * p and q */
	int pi, pj, qi, qj;
	double drp, drq, dx, dy, dz;
	double qp, qq;
	int ind;
	for (int i=1; i<natom; i++) {
		for (int j=0; j<i; j++) {
			ind = index(i,j);
			if (ind == p) {
				pi = i;
				pj = j;
			} else if (ind == q) {
				qi = i;
				qj = j;
			}
		}
	}
	if (p == q) {
		qi = pi;
		qj = pj;
	}
	for (int s=0; s<nsave; s++) {
		dx = X[s*natom+pi]-X[s*natom+pj];
		dy = Y[s*natom+pi]-Y[s*natom+pj];
		dz = Z[s*natom+pi]-Z[s*natom+pj];
		drp = sqrt(dx*dx+dy*dy+dz*dz);
		qp = contact(drp, rc_p);
		dx = X[s*natom+qi]-X[s*natom+qj];
		dy = Y[s*natom+qi]-Y[s*natom+qj];
		dz = Z[s*natom+qi]-Z[s*natom+qj];
		drq = sqrt(dx*dx+dy*dy+dz*dz);
		qq = contact(drq, rc_q);
		nu_proj[s] = proj[s] + qp*dp + qq*dq;
	}
}

int checkpos(proj_type *proj, int n)
{
	for (int i=0; i<n; i++) {
		if (proj[i] < 0.0)
			return 1;
	}
	return 0;
}

void proj_calc_active(vector<proj_type> &c, vector<proj_type> &rcut, vector<int> &a, 
		crd_type *X, crd_type *Y, crd_type *Z, int natom,
		int nsave, proj_type *proj)
{
	vector<proj_type> p;
	proj_type dx,dy,dz,dr,r,eqrc;
	int ind, i,j;
	int dim = c.size();
	p.resize(dim);
	for (int s=0; s<nsave; s++) {
		for (int ind=0; ind<dim; ind++) {
			get_ij(a[ind],i,j);
			dx = X[s*natom+i]-X[s*natom+j];
			dy = Y[s*natom+i]-Y[s*natom+j];
			dz = Z[s*natom+i]-Z[s*natom+j];
			dr = sqrt(dx*dx+dy*dy+dz*dz);
			p[ind] = contact(dr,rcut[ind]);
		}
		r = dotp(c,p);
		proj[s] = r;
	}
	return;
}

void proj_calc(vector<proj_type> &c, crd_type *X, crd_type *Y, crd_type *Z, int natom,
		int nsave, proj_type *proj)
{
	vector<proj_type> p;
	double dx,dy,dz,dr,r,eqrc;
	int ind;
	int dim = c.size();
	p.resize(dim);
	for (int s=0; s<nsave; s++) {
		for (int i=1; i<natom; i++) {
			for (int j=0; j<i; j++) {
				ind = index(i,j);
				dx = X[s*natom+i]-X[s*natom+j];
				dy = Y[s*natom+i]-Y[s*natom+j];
				dz = Z[s*natom+i]-Z[s*natom+j];
				dr = sqrt(dx*dx+dy*dy+dz*dz);
				p[ind] = contact(dr,12.0);
			}
		}
		r = dotp(c,p);
		proj[s] = r;
	}
	return;
}

void read_active_vec(vector<proj_type> &v, vector<proj_type> &rcut,
		vector<int> &active_ind, const char *file) 
{
	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
	int tmpi;
	double tmpf;
	FILE *f = fopen(file,"r");
	if (f == NULL) {
		fprintf(stderr,"Could not open file %s to read vector\n",file);
		exit(1);
	}
	fgets(buf,buf_len,f);
	int dim = 0;
	while (feof(f) == 0) {
		dim++;
		fgets(buf,buf_len,f);
	}
	fclose(f);
	// allocate space
	v.resize(dim);
	active_ind.resize(dim);
	rcut.resize(dim);
	// read conts
	f = fopen(file,"r");
	double norm = 0.0;
	// file format has three columns: col1 = contact index, col2 = weight, col3 = cutoff
	// ** only active weights are saved
	for (int i=0; i<dim; i++) {
		fgets(buf,buf_len,f);
		tmpi = atoi(strtok(buf," \t"));
		tmpf = atof(strtok(NULL," \t"));
		active_ind[i] = tmpi;
		v[i] = tmpf;
		norm += tmpf*tmpf;
		rcut[i] = atof(strtok(NULL," \t"));
		fprintf(stdout,"%6i %12.8e %12.5f\n",active_ind[i],v[i],rcut[i]);
	}
	// ensure norm = 1 (since reading from text file could be small errors)
	/*
	norm = sqrt(norm);
	for (int i=0; i<dim; i++) {
		v[i] /= norm;
	}
	fclose(f);
	*/
	return;
}

int main(int argc, char **argv)
{
	park_miller uniform;
	box_muller gaussian;
	//DCDITrajFile inptraj;
	BaseITrajFile *inptraj;
	bool change_cutoffs = false;
	bool noswap = false;
	bool posweights = false;
	int nmove = 7;
	double weightsum;
	int natom, dim, move_type;		
	vector<proj_type> pvec, rcut_vec, best_pvec, best_rcut_vec;
	vector<int> active_ind, best_active_ind;
	double ptot, ptot_trial, psdev, psdev_trial, pmean, pmean_trial, crit, crit_trial, best_crit;
	vector<string> dcd_names, istp_names;
	string pdb_name, oroot, qfile,pvec_traj_name, contact_file, pvec_tmpf;
	int skip, c, nbin, ntraj, nsave, n_proj_save, max_save, tmpi, tmpj;
	double bin_lo, bin_hi;
	crd_type *X, *Y, *Z;
	proj_type *rc, *proj, *nu_proj, *pr1, *pr2, *pv_tmp, *rcut_tmp;
	int *ai_tmp;
	int *isTP;
	string  pvec_file = "NONE";
	string  istp_file = "NONE";
	vector<int> ii,jj;
	const int max_pl = 100000;
	double mu = -1.0;	// no gcmc
	int pl;
	vector<proj_type> p_r, p_rTP, p_TPr;
	double alpha = 0.2;
	int nmc = 1000;
	int SOTRUE = 1;
	int SOFALSE = 0;
	double lambda;
	double mean_a,var_a,mean_b,var_b;
	double T, rc_i, rc_i_nu, w_i_nu;
	double default_cutoff = 12.0;		// Angstroms
	double default_cutoff_sdev = -1.0;		// Angstroms
	double sigma_rc = 2.0;
	int move_attempts[nmove], move_successes[nmove];
	double avew, sdevw;
	int random_seed = 100;
	char decider = 'q';		// = 'q' (integrated pTP/sdev)
					// = 'i' (integrated pTP)
					// = 'm' (maximum)
					// = 'g' (maximum from gaussian fit)
	long clock_start, clock_stop, us;
	int outp_freq, active_contacts, active_contacts_trial;
	// active moves = 3 for all contacts and 4 if # active contacts < total contacts
	FILE *proj_out, *pvec_out, *pvec_traj; //, *ptp_out;
	string proj_name, pvec_name, ptp_name;

	curtime(clock_start,us);
	if (argc < 3) {
		fprintf(stderr,Usage);
		exit(1);
	}
	active_contacts = 0;
	skip = 10;
	T = 0.0;	// default reduces to original algorithm
	outp_freq = 100;
	while (1) {
		c=getopt(argc,argv,"hL:H:n:m:s:p:c:o:q:d:r:T:f:N:x:FPg:K:k:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'L':
				bin_lo = atof(optarg);
				break;
			case 'H':
				bin_hi = atof(optarg);
				break;
			case 'K':
				default_cutoff = atof(optarg);
				break;
			case 'k':
				default_cutoff_sdev = atof(optarg);
				break;
			case 'g':
				mu = atof(optarg);
				break;
			case 'n':
				nbin = atoi(optarg);
				break;
			case 'm':
				nmc = atoi(optarg);
				break;
			case 's':
				skip = atoi(optarg);
				break;
			case 'p':
				pdb_name = optarg;
				break;
			case 'c':
				pvec_file = optarg;
				break;
			case 'q':
				contact_file = optarg;
				break;
			case 'o':
				oroot = optarg;
				break;
			case 'd':
				decider = optarg[0];
				break;
			case 'T':
				T = atof(optarg);
				break;
			case 'r':
				random_seed = atoi(optarg);
				break;
			case 'f':
				outp_freq = atoi(optarg);
				break;
			case 'N':
				active_contacts = atoi(optarg);
				break;
			case 'x':
				change_cutoffs = true;
				sigma_rc = atof(optarg);
				break;
			case 'F':
				noswap = true;
				break;
			case 'P':
				posweights = true;
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	uniform.initialize(random_seed);
	ntraj = (argc-optind)/2;	// half dcd, half istp
	dcd_names.resize(ntraj);
	istp_names.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		dcd_names[i] = argv[optind+2*i];
		istp_names[i] = argv[optind+2*i+1];
	}
	inptraj = open_trajectory(dcd_names[0].c_str());
	natom = inptraj->num_atoms();
	inptraj->close();
	dim = natom*(natom-1)/2;
	if (active_contacts > dim) {
		fprintf(stderr,"Number of active contacts cannot be larger than contact vector dimension!\n");
		exit(1);
	}
	//pvec.resize(active_contacts);
	//rcut_vec.resize(active_contacts);
	//best_pvec.resize(active_contacts);
	//best_rcut_vec.resize(active_contacts);
	//active_ind.resize(active_contacts);
	//best_active_ind.resize(active_contacts);
	if (string("NONE") == pvec_file) { 		// starting from random initial contacts
		if (active_contacts == 0 || active_contacts == dim) {
			active_contacts = dim;
			fprintf(stdout,">> Setting active contacts to all contacts (%i)\n",dim);
		} else if ( active_contacts > dim ) {
			fprintf(stderr,"?! Number of active contacts (%i) cannot be larger than contact vector (%i)!\n",
					active_contacts,dim);
			exit(1);
		}
		fprintf(stdout,">> Generating random initial weights\n");
		generate_random_vec(pvec,rcut_vec,active_ind,active_contacts,
				dim,uniform,default_cutoff,posweights);
		pvec_name = oroot + string("_pvec_random.dat");
		pvec_out = fopen(pvec_name.c_str(),"w"); 
		for (int i=0; i<active_contacts; i++) {
			fprintf(pvec_out,"%6i %12.8e %12.5f\n",active_ind[i],pvec[i],rcut_vec[i]);
		}
		fclose(pvec_out);
	} else {					// read contacts from file
		fprintf(stdout,">> Reading initial weights from file %s\n",pvec_file.c_str());
		read_active_vec(pvec,rcut_vec,active_ind,pvec_file.c_str());
		int nc_read = active_ind.size();
		fprintf(stdout,">> Number of contacts in file = %i\n",nc_read);
		if (active_contacts == 0) {				// default - take number from file
			fprintf(stdout,">> Taking contact vector size from file\n");
			active_contacts = active_ind.size();
		} else if (active_contacts< active_ind.size()) { 	// must have at least number from file
			fprintf(stderr,"?! number of active contacts smaller than vector in file!\n");
			exit(1);
		} else if (active_contacts > active_ind.size()){	// extend rc size
			fprintf(stdout,">> Extending contact vector in file %s (size %i)\n",
					pvec_file.c_str(),nc_read);
			fprintf(stdout,">> to size %i by adding random contacts with zero weight\n",
					active_contacts);
			extend_rc_vec(pvec, rcut_vec, active_ind, 
					active_contacts, dim, nc_read, 
					uniform, default_cutoff);
		}
	}
	fflush(stdout);
	fprintf(stdout,"GOT HERE 0\n");
	best_rcut_vec.resize(active_contacts);
	best_pvec.resize(active_contacts);
	best_active_ind.resize(active_contacts);
	if (posweights) {
		for (int i=0; i<active_contacts; i++) {
			pvec[i] = abs(pvec[i]);
		}
	}

	if (active_contacts == dim) 
		noswap = true;
	pv_tmp = new proj_type[dim];
	rcut_tmp = new proj_type[dim];
	ai_tmp = new int[dim];
	fprintf(stdout,"GOT HERE 0.5\n");
	for (int i=0; i<active_contacts; i++) {			// worst case scenario!
		best_active_ind[i] = active_ind[i];		
		best_pvec[i] = pvec[i];				
		best_rcut_vec[i] = rcut_vec[i];				
	}
	fprintf(stdout,"GOT HERE 1\n");
	max_save = 0;
	fprintf(stdout,"=========================================================\n");
	fprintf(stdout,"                PDB file = %s\n", pdb_name.c_str());
	fprintf(stdout,"         Number of atoms = %i\n", natom);
	fprintf(stdout,"Contact vector dimension = %i\n", dim);
	fprintf(stdout,"Number of Active Weights = %i\n", active_contacts);
	fprintf(stdout,"Trajectory files:\n");
	for (int i=0; i<ntraj; i++) {
		fprintf(stdout,"%s\n",dcd_names[i].c_str());
		inptraj = open_trajectory(dcd_names[i].c_str());
		max_save += inptraj->total_frames();
		inptraj->close();
	}
	fprintf(stdout,"isTP files:\n");
	for (int i=0; i<ntraj; i++) {
		fprintf(stdout,"%s\n",istp_names[i].c_str());
	}
	fprintf(stdout,"Total number of frames = %i\n",max_save);
	fprintf(stdout,"           Random seed = %i\n",random_seed);
	max_save = max_save/skip + 100;
	fprintf(stdout,"# Frames to save = %i (skipping every %i)\n",max_save,skip);
	if (noswap) {
		fprintf(stdout,"Not swapping active contacts!\n");
	}
	fprintf(stdout,"=========================================================\n");
	X = new crd_type[max_save*natom];
	Y = new crd_type[max_save*natom];
	Z = new crd_type[max_save*natom];
	rc = new proj_type[max_save];
	isTP = new int[max_save];
	pr1 = new proj_type[max_save];
	pr2 = new proj_type[max_save];
	proj = pr1;
	nu_proj = pr2;
	int bpath;
	double eqrc;
	if (read_trajectories(dcd_names,X,Y,Z,natom,skip,max_save,nsave)) {
		fprintf(stderr,"nsave > maxsave!!! (How does that happen?)\n");
		delete [] X;
		delete [] Y;
		delete [] Z;
		exit(1);
	}
	read_istp(istp_names,isTP,skip,nsave);

	// calculate initial pTPr
	proj_calc_active(pvec, rcut_vec, active_ind, X, Y, Z, natom, nsave, proj);
	pTPr_calc(isTP,proj, nsave, bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr);
	switch (decider) {
		case 'g':
			crit = pTP_gauss(p_TPr);
			break;
		case 'i':
			crit = pTP_sum(p_TPr);
			break;
		case 'q':
			crit = pTP_quotient(p_TPr);
			break;
		case 'm':
			crit = pTP_max(p_TPr);
			break;
		case 'a':
			crit = pTP_gauss_eq(p_TPr, p_r);
			break;
		case 'b':
			crit = pTP_quotient_eq(p_TPr, p_r);
			break;
		case 'h':
			crit = pTP_gauss_half(p_TPr,bin_lo,bin_hi,nbin);
			break;
		case 'A':
			crit = attila_crit(isTP,proj,nsave,mean_a,var_a,mean_b,var_b);
			break;
		default:
			fprintf(stderr,"Unknown optimization criterion %1c\n",decider);
			exit(1);
			break;
	}
	if (mu>=0) 
		crit *= exp(-mu*active_contacts);
	// write projection onto initial coordinate
	proj_name = oroot + string("_proj_ini.dat");
	proj_out = fopen(proj_name.c_str(),"w"); 
	for (int i=0; i<nsave; i++) {
		fprintf(proj_out,"%12.8e\n",proj[i]);
	}
	fclose(proj_out);
	ptp_name = oroot + string("_pTP_ini.dat");
	write_ptp(ptp_name.c_str(), bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr);
	pvec_traj_name = oroot + string("_pvec_trj.dat");
	pvec_traj = fopen(pvec_traj_name.c_str(),"w");
	n_proj_save = nmc/outp_freq;	// number of frames to save
	fwrite(&n_proj_save,sizeof(int),1,pvec_traj);
	fwrite(&dim,sizeof(int),1,pvec_traj);
	fwrite(&active_contacts,sizeof(int),1,pvec_traj);
	if (change_cutoffs) {
		fwrite(&SOTRUE,sizeof(int),1,pvec_traj);
	} else {
		fwrite(&SOFALSE,sizeof(int),1,pvec_traj);
		proj_type *tmpvr = new proj_type[active_contacts];
		for (int i=0; i<active_contacts; i++) 
			tmpvr[i] = rcut_vec[i];
		fwrite(tmpvr,sizeof(proj_type),active_contacts,pvec_traj);
		delete [] tmpvr;
	}
	if (noswap) {
		fwrite(&SOTRUE,sizeof(int),1,pvec_traj);
		int *tmpvi = new int[active_contacts];
		for (int i=0; i<active_contacts; i++) 
			tmpvi[i] = active_ind[i];
		fwrite(tmpvi,sizeof(int),active_contacts,pvec_traj);
		delete [] tmpvi;
	} else {
		fwrite(&SOFALSE,sizeof(int),1,pvec_traj);
	}
	if (mu > 0.0) { //prepare for gcmc
		avew = 0.0;
		sdevw = 0.0;
		for (int i=0; i<active_contacts; i++) {
			avew += pvec[i];
			sdevw += pvec[i]*pvec[i];
		}
		avew /= double(active_contacts);
		sdevw = sqrt(sdevw/double(active_contacts)-avew*avew);
	}
	best_crit = 0.0;
	int elem_i, elem_j;
	double pi, pj, pi2, pj2, pi_nu, pj_nu;
	double mag;
	for (int i=0; i<nmove; i++) {
		move_attempts[i] = 0;
		move_successes[i] = 0;
	}
	for (int s = 0; s<nmc; s++) {
		move_type = int(floor(uniform.random()*nmove));
		// keep choosing till we don't have a conflict
		while ((move_type==3 && noswap) || (move_type==4 && (!change_cutoffs))
				|| (move_type==2 && posweights) 
				|| ((move_type==5 || move_type==6) && mu<=0.0)) { 
			move_type = int(floor(uniform.random()*nmove));
		}

		//move_type = 3;
		move_attempts[move_type]++;
		if (move_type == 3) {
			// change active elements
			// choose an active element elem_i
			elem_i = int(floor(uniform.random()*active_contacts));
			tmpi = active_ind[elem_i];
			// check that tmpj is not an active element
			bool active;
			active = true;
			while (active) {
				tmpj = int(floor(uniform.random()*dim));
				active = false;
				for (int i=0; i<active_contacts; i++) {
					if (active_ind[i] == tmpj) {
						active = true;
						break;
					}
				}
			}
			proj_update(X, Y, Z, natom, nsave, proj, nu_proj, 
					tmpi,tmpj,-pvec[elem_i],pvec[elem_i],
					rcut_vec[elem_i],rcut_vec[elem_i]); 
		} else if (move_type==4) {
			// change cut-off
			elem_i = int(floor(uniform.random()*active_contacts));
			tmpi = active_ind[elem_i];
			rc_i = rcut_vec[elem_i];
			rc_i_nu = -1.0;
			while (rc_i_nu <= 0.0)
				rc_i_nu = rc_i + gaussian.random()*sigma_rc;
			proj_update(X, Y, Z, natom, nsave, proj, nu_proj, 
					tmpi,tmpi,-pvec[elem_i],pvec[elem_i],
					rc_i,rc_i_nu); 
		} else if (move_type==5) {
			// grand-canonical insertion
			// choose an 'inactive' contact
			bool active = true;
			while (active) {
				tmpi = int(floor(uniform.random()*dim));
				active = false;
				for (int i=0; i<active_contacts; i++) {
					if (active_ind[i] == tmpi) {
						active = true;
						break;
					}
				}
			}
			if (default_cutoff_sdev > 0) {
				rc_i_nu = -1.0;
				while (rc_i_nu <= 0.0)
					rc_i_nu = default_cutoff + gaussian.random()*default_cutoff_sdev;
			} else {
				rc_i_nu = default_cutoff;
			}
			w_i_nu = -1.0;
			while (w_i_nu <= 0.0)
				w_i_nu = avew + gaussian.random()*sdevw;
			proj_update(X, Y, Z, natom, nsave, proj, nu_proj, 
					tmpi,0,w_i_nu,0.0,
					rc_i_nu,10.0); 
			active_contacts_trial = active_contacts+1;
		} else if (move_type==6) {
			// grand-canonical deletion
			// choose active contact
			elem_i = int(floor(uniform.random()*active_contacts));
			tmpi = active_ind[elem_i];
			proj_update(X, Y, Z, natom, nsave, proj, nu_proj, 
					tmpi,0,-pvec[elem_i],0.0,
					rcut_vec[elem_i],10.0); 
			active_contacts_trial = active_contacts-1;
		} else {
			choose2(uniform,active_contacts,elem_i,elem_j);
			pi = pvec[elem_i]; pi2 = pi*pi;
			pj = pvec[elem_j]; pj2 = pj*pj;
			if (move_type == 0) {
				// reassign relative weights
				lambda = uniform.random()*(pi2+pj2);
				pi_nu = sqrt(lambda)*sign(pi);
				pj_nu = sqrt(pi2+pj2-lambda)*sign(pj);
			} else if (move_type == 1) {
				// swap weights
				pi_nu = pj;
				pj_nu = pi;
			} else if (move_type == 2) {
				// negate weights
				pi_nu = -pi;
				pj_nu = pj;
			}
			proj_update(X, Y, Z, natom, nsave, proj, nu_proj, 
					active_ind[elem_i],active_ind[elem_j],
					pi_nu-pi,pj_nu-pj,
					rcut_vec[elem_i],rcut_vec[elem_j]); 
		}
		// DEBUG -------------------------------------
		pvec_tmpf = oroot + string("_pvec_sofar.dat");
		weightsum = save_pvec(pvec_tmpf.c_str(), active_ind, pvec, rcut_vec, 
				active_contacts); 
		fprintf(stdout,"move %i; wsum = %12.6e; move type = %i\n",s,weightsum,move_type);
		// DEBUG -------------------------------------
		if (decider != 'A') 
			pTPr_calc(isTP,nu_proj, nsave, bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr);
		switch (decider) {
			case 'g':
				crit_trial = pTP_gauss(p_TPr);
				break;
			case 'i':
				crit_trial = pTP_sum(p_TPr);
				break;
			case 'q':
				crit_trial = pTP_quotient(p_TPr);
				break;
			case 'm':
				crit_trial = pTP_max(p_TPr);
				break;
			case 'a':
				crit_trial = pTP_gauss_eq(p_TPr, p_r);
				break;
			case 'b':
				crit_trial = pTP_quotient_eq(p_TPr, p_r);
				break;
			case 'h':
				crit_trial = pTP_gauss_half(p_TPr,bin_lo,bin_hi,nbin);
				break;
			case 'A':
				crit_trial = attila_crit(isTP,nu_proj,nsave,mean_a,var_a,mean_b,var_b);
				break;
			default:
				fprintf(stderr,"Unknown optimization criterion %1c\n",decider);
				exit(1);
				break;
		}
		// shouldn't this be crit_trial/pow(active_contacts_trial,mu)
		if (mu > 0) 	// gcmc
			crit_trial /= pow(active_contacts_trial,mu);
		// old version!
		//	crit_trial *= exp(-mu*active_contacts_trial);
			
		if (crit_trial > crit) {
			proj_type *tmp = proj;
			proj = nu_proj;
			nu_proj = tmp;
			crit = crit_trial;
			if (move_type == 3) {
				active_ind[elem_i] = tmpj;
			} else if (move_type == 4) {
				rcut_vec[elem_i] = rc_i_nu;
			} else if (move_type == 5) {
				active_contacts = active_contacts_trial;
				// expand vectors
				pvec.resize(active_contacts);
				rcut_vec.resize(active_contacts);
				active_ind.resize(active_contacts);
				// insert
				pvec[active_contacts-1] = w_i_nu;
				rcut_vec[active_contacts-1] = rc_i_nu;
				active_ind[active_contacts-1] = tmpi;
			} else if (move_type == 6) {
				active_contacts = active_contacts_trial;
				// shift vectors
				for (int i=elem_i; i<active_contacts-1; i++) {
					pvec[i] = pvec[i+1];
					rcut_vec[i] = rcut_vec[i+1];
					active_ind[i] = active_ind[i+1];
				}
				// shrink
				pvec.resize(active_contacts);
				rcut_vec.resize(active_contacts);
				active_ind.resize(active_contacts);
			} else {
				pvec[elem_i] = pi_nu;
				pvec[elem_j] = pj_nu;
			}
			move_successes[move_type]++;
			if (crit>best_crit) {
				// save new best weight matrix - write out at end
				if (best_active_ind.size() != active_contacts) {
					best_active_ind.resize(active_contacts);
					best_pvec.resize(active_contacts);
					best_rcut_vec.resize(active_contacts);
				}
				for (int i=0; i<active_contacts; i++) {
					best_active_ind[i] = active_ind[i];
					best_pvec[i] = pvec[i];
					best_rcut_vec[i] = rcut_vec[i];
				}
				pvec_tmpf = oroot + string("_pvec_sofar.dat");
				weightsum = save_pvec(pvec_tmpf.c_str(), best_active_ind, best_pvec, best_rcut_vec, 
						active_contacts); 
				ptp_name = oroot + string("_pTP_sofar.dat");
	        		write_ptp(ptp_name.c_str(), bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr);
				best_crit = crit;
				fprintf(stdout,"!! New best weight matrix: move # %i  crit = %12.5f  contacts = %i  wsum=%12.5f\n", 
						s, crit, active_contacts,weightsum);
				/*
				if (posweights) {
					int pos_proj = checkpos(proj,nsave);
					if (pos_proj) {
						fprintf(stdout,"Some negative projections found, recalculating...\n");
						proj_calc_active(pvec, rcut_vec, active_ind, X, Y, Z, natom, nsave, proj);
						pos_proj = checkpos(proj,nsave);
						if (pos_proj) {
							fprintf(stdout,"Still negative pro projections found, quitting\n");
							exit(1);
						}
					}
				}
				*/
			}
		} else if (T>1e-20) {
			double rand = uniform.random();
			if (rand < pow(crit_trial/crit,1.0/T)) {
				proj_type *tmp = proj;
				proj = nu_proj;
				nu_proj = tmp;
				crit = crit_trial;
				if (move_type == 3) {
					active_ind[elem_i] = tmpj;
				} else if (move_type == 4) {
					rcut_vec[elem_i] = rc_i_nu;
				} else if (move_type == 5) {
					active_contacts = active_contacts_trial;
					// expand vectors
					pvec.resize(active_contacts);
					rcut_vec.resize(active_contacts);
					active_ind.resize(active_contacts);
					// insert
					pvec[active_contacts-1] = w_i_nu;
					rcut_vec[active_contacts-1] = rc_i_nu;
					active_ind[active_contacts-1] = tmpi;
				} else if (move_type == 6) {
					active_contacts = active_contacts_trial;
					// shift vectors
					for (int i=elem_i; i<active_contacts-1; i++) {
						pvec[i] = pvec[i+1];
						rcut_vec[i] = rcut_vec[i+1];
						active_ind[i] = active_ind[i+1];
					}
					// shrink
					pvec.resize(active_contacts);
					rcut_vec.resize(active_contacts);
					active_ind.resize(active_contacts);
				} else {
					pvec[elem_i] = pi_nu;
					pvec[elem_j] = pj_nu;
				}
				move_successes[move_type]++;
			}
		} 
		// write this once we know whether it's accepted!
		if (s%outp_freq == 0) {
			fprintf(stdout,"MC move # %i  crit = %12.8f  contacts = %i\n", s, crit,active_contacts);
			for (int i=0; i<active_contacts; i++) {
				pv_tmp[i] = pvec[i];
				ai_tmp[i] = active_ind[i];
				rcut_tmp[i] = rcut_vec[i];
			}
			if (mu > 0.0) 
				fwrite(&active_contacts,sizeof(int),1,pvec_traj);
			fwrite(pv_tmp,sizeof(proj_type),active_contacts,pvec_traj);
			if (!noswap)
				fwrite(ai_tmp,sizeof(int),active_contacts,pvec_traj);
			if (change_cutoffs) 
				fwrite(rcut_tmp,sizeof(proj_type),active_contacts,pvec_traj);
		}
	}
	fprintf(stdout,"MC run finished; crit = %12.8f\n", crit);
	fprintf(stdout,"%4s  %8s  %8s  %12s\%\n", "move", 
			"attempts", "successes", "% success");
	for (int i=0; i<nmove; i++) {
		double psuc=0.0;
		if (move_attempts[i]>0) 
			psuc = double(move_successes[i])/double(move_attempts[i])*100.0;
		fprintf(stdout,"%4i  %8i  %8i  %12.5f\%\n", i, 
				move_attempts[i], move_successes[i], psuc);
	}
	proj_calc_active(best_pvec, best_rcut_vec, best_active_ind, X, Y, Z, natom, nsave, proj);
	pTPr_calc(isTP,proj, nsave, bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr);
	ptp_name = oroot + string("_pTP_fin.dat");
	write_ptp(ptp_name.c_str(), bin_lo, bin_hi, nbin, p_r, p_rTP, p_TPr);
	pvec_name = oroot + string("_pvec_fin.dat");
	pvec_out = fopen(pvec_name.c_str(),"w"); 
	double norm = 0.0;
	for (int i=0; i<active_contacts; i++) {
		fprintf(pvec_out,"%6i %12.8e %12.5f\n",best_active_ind[i],best_pvec[i],best_rcut_vec[i]);
		norm += best_pvec[i]*best_pvec[i];
	}
	fclose(pvec_out);
	fprintf(stdout,"norm(pvec) = %12.8f\n", sqrt(norm));
	proj_name = oroot + string("_proj_fin.dat");
	proj_out = fopen(proj_name.c_str(),"w"); 
	for (int i=0; i<nsave; i++) {
		fprintf(proj_out,"%12.8e\n",proj[i]);
	}
	fclose(proj_out);
	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] rc;
	delete [] isTP;
	delete [] pr1;
	delete [] pr2;
	delete [] pv_tmp;
	delete [] rcut_tmp;

	curtime(clock_stop,us);
	long time = clock_stop-clock_start;
	int hours = time/3600;
	int minutes = (time%3600)/60;
	int seconds = time%60;
	fprintf(stdout,"Run time = %i hours %i mins %i secs\n", hours, minutes, seconds);

	return 0;
}
