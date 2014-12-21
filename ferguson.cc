/* 
 * design a protein to be two-state, no matter what!
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
#include "random_gen.h"
#include "config.h"
#include <sys/time.h>

void choose2(park_miller &gen, int nc, int &elem_i, int &elem_j)
{
	elem_i = elem_j = int(floor(gen.random()*nc));
	while (elem_j == elem_i) {
		elem_j = int(floor(gen.random()*nc));
	}
	return;
}

void read_golist(vector<int> &i, vector<int> &j, vector<proj_type> &rij, 
		vector<proj_type> &eij, const char *file) 
{
	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
	int tmpi;
	double tmpf;
	FILE *f = fopen(file,"r");
	if (f == NULL) {
		fprintf(stderr,"Could not open file %s to read go contacts\n",file);
		exit(1);
	}
	fgets(buf,buf_len,f);
	int dim = 0;
	// count contacts
	while (feof(f) == 0) {
		dim++;
		fgets(buf,buf_len,f);
	}
	fclose(f);
	// allocate space
	i.resize(dim);
	j.resize(dim);
	rij.resize(dim);
	eij.resize(dim);
	// read contacts
	f = fopen(file,"r");
	double norm = 0.0;
	// file format has three columns: col1 = contact index, col2 = weight, col3 = cutoff
	// ** only active weights are saved
	for (int p=0; p<dim; p++) {
		fgets(buf,buf_len,f);
		i[p] = atoi(strtok(buf," \t"))-1;
		j[p] = atoi(strtok(NULL," \t"))-1;
		rij[p] = atof(strtok(NULL," \t"));
		eij[p] = atof(strtok(NULL," \t"));
		//fprintf(stdout,"%4i %4i %12.6f %12.6f\n",i[p],j[p],rij[p],eij[p]);
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

void read_rc(vector<string> &q_names, proj_type *rc, int nframe)
{
	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
	int tmpi;
	FILE *f;
	double tmpf;
	int nf =  q_names.size();
	int ind = 0;
	for (int ff = 0; ff<nf; ff++) {
		f = fopen(q_names[ff].c_str(),"r");
		if (f == NULL) {
			fprintf(stderr,"Could not open file %s to read go contacts\n",
					q_names[ff].c_str());
			exit(1);
		}
		double norm = 0.0;
		// rc in col 2
		while (feof(f) == 0) {
			fgets(buf,buf_len,f);
			if (feof(f) != 0) break;
			tmpi = atoi(strtok(buf," \t"));
			rc[ind++] = atof(strtok(NULL," \t"));
			//fprintf(stdout,"rc[%i] = %8.3f\n",ind,rc[ind-1]);
		}
		fclose(f);
	}
	return;
}
void calc_contact_energy(vector<string> &dcdfiles, vector<int> &ii, vector<int> &jj,
		vector<proj_type> &rij, proj_type *EIJ, char ptype)
{
	int cumframe, nc, cumi, natom,i,j;;
	nc = ii.size();
	DCDITrajFile inptraj;
	double ee,sr,dx,dy,dz,dr,sig;
	double *X,*Y,*Z;

	inptraj.open(dcdfiles[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	X = new double[natom];
	Y = new double[natom];
	Z = new double[natom];

	cumframe = 0;
	for (int f=0; f<dcdfiles.size(); f++) {
		inptraj.open(dcdfiles[f].c_str());
		int nframe = inptraj.total_frames();
		for (int frame = 0; frame < nframe; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			for (int p=0; p<nc; p++) {
				i = ii[p];
				j = jj[p];
				sig = rij[p];
				dx = X[i]-X[j];
				dy = Y[i]-Y[j];
				dz = Z[i]-Z[j];
				dr = sqrt(dx*dx+dy*dy+dz*dz);
				if (ptype == 'm') {
					ee = (1.0-exp(-1.7*(dr-sig)));
					ee = ee*ee;
				} else if (ptype == 'k') {
					sr = sig/dr;
					ee = 13.*pow(sr,12)-18.*pow(sr,10)
						+4.*pow(sr,6);
				}
				EIJ[cumframe*nc+p] = ee;
				if (cumframe == 0) {
					fprintf(stdout,"%4i %4i %8.3f\n",i,j,ee);
				}
			}
			cumframe++;
		}
		inptraj.close();
	}

	delete [] X;
	delete [] Y;
	delete [] Z;
	return;

}

void update_weights(proj_type *EIJ, proj_type *w, proj_type *trial_w,
		int ind, proj_type dE, int nframe, int nc, double beta)
{
	for (int f=0; f<nframe; f++) {
		trial_w[f] = w[f]*exp(-beta*dE*EIJ[f*nc+ind]);
	}
}

double opt_crit(proj_type *w, proj_type *q, const double &x, const int &n)
{
	double Nf, Nu, mf, mu, sf, su, ww, qq, dm, Nts, NN;

	Nf = 0.;
	Nu = 0.;
	Nts = 0.;
	NN = 0.;
	mu = 0.;
	mf = 0.;
	su = 0.;
	sf = 0.;
	for (int i=0; i<n; i++) {
		qq = q[i];
		ww = w[i];
		//fprintf(stdout,"qq=%12.6f, x=%12.6f, ww=%12.6f\n",qq, x, ww);
		if (qq<x) {
			//fprintf(stdout,"<<<<<<<<<<<<\n");
			Nu += ww;
			mu += qq*ww;
			su += qq*qq*ww;
		} else {
			//fprintf(stdout,">>>>>>>>>>>>\n");
			Nf += ww;
			mf += qq*ww;
			sf += qq*qq*ww;
		}
		if (qq<0.60&&qq>0.4) {
			Nts += ww;
		}
		NN += ww;
	}
	//fprintf(stdout,"mu=%12.6f; mf=%12.6f; Nu=%12.6f; Nf=%12.6f\n",mu,mf,Nu,Nf);
	mu /= Nu;
	su = su/Nu - mu*mu;
	mf /= Nf;
	sf = sf/Nf - mf*mf;
	dm = mu-mf;
	fprintf(stdout,"mu=%12.6f; mf=%12.6f; su=%12.6f; sf=%12.6f\n",mu,mf,su,sf);
	return dm*dm/(sf+su);
	//return dm*dm; //*NN/Nts;
}

using namespace std;

const char *Usage = "\n"
"	Usage: \n"
"		ferguson [-r random_seed] [-d A] \n"
"			[-s skip] [-T temp] [-B MC temp]\n"
"			-f folded_lo -F folded_hi -u unfolded_lo -U unfolded_hi \n"
"			-g golist -m n_mc -n outp_freq\n"
"			-o outp_root\n"
"			file1.dcd q1.dat file2.dcd q2.dat... fileN.dcd qN.dat\n"
"	where:\n"
"		-d <opt> determines the optimization criterion as follows:\n"
"			<opt>=A: Optimize (<r>_f-<r>_u)^2/(var(r)_f+var(r)_u) [default]\n"
"               -x dividing value between U and F\n"
"               -M maximum fractional change in parameters\n"
"               -S standard deviation of trial changes in parameters\n"
"		-T temp is the temperature at which the simulation was run\n"
"		-B MC temp is the temperature for MC sampling\n"
"		-n outp_freq: frequency for writing MC output\n"
"		-g golist: file with go interactions\n"
"		-m nmc is the number of optimization iterations to do\n"
"		-o outp_root is a base name for output files\n"
"		-s skip: only use every skipth frame\n"
"	Trajectories do not have to be contiguous as they are assumed independent\n\n";

int main(int argc, char **argv)
{
	park_miller uniform;
	box_muller gaussian;
	int c,skip,nprint,nmc,nframe,random_seed,ntraj,natom,nc;
	double simT, mcT, flo, fhi, ulo, uhi, divide, maxchange, sigchange;
	double tmpe,crit,trial_crit,simbeta;
	char decider,ptype;
	proj_type *EIJ, *rc, *w1, *w2, *w, *trial_w;
	DCDITrajFile inptraj;
	vector<int> ii,jj;
	vector<proj_type> rij,eij,eij_orig,eij_best;
	vector<string> dcd_names, q_names;
	double binlo, binhi;
	int nbin;
	string go_file,oroot;
	if (argc < 3) {
		fprintf(stderr,Usage);
		exit(1);
	}
	// defaults
	ptype = 'm';
	decider = 'A';
	skip = 1;
	simT = 300.;
	mcT = 1.0;
	flo = 0.8;
	fhi = 1.0;
	ulo = 0.0;
	uhi = 0.4;
	nmc = 1000;
	nprint = 100;
	oroot = "junk";
	random_seed = 27041994;
	sigchange = 0.05;
	maxchange = 0.2;
	nbin = 40;
	binlo = binhi = 0.;

	while (1) {
		c=getopt(argc,argv,"hd:s:T:B:x:g:m:n:o:r:F:p:L:H:N:M:S:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'd':
				decider = optarg[0];
				break;
			case 'p':
				ptype = optarg[0];
				break;
			case 's':
				skip = atoi(optarg);
				break;
			case 'T':
				simT = atof(optarg);
				break;
			case 'B':
				mcT = atof(optarg);
				break;
			case 'N':
				nbin = atoi(optarg);
				break;
			case 'L':
				binlo = atof(optarg);
				break;
			case 'H':
				binhi = atof(optarg);
				break;
			case 'M':
				maxchange = atof(optarg);
				break;
			case 'S':
				sigchange = atof(optarg);
				break;
				/*
			case 'f':
				flo = atof(optarg);
				break;
			case 'F':
				fhi = atof(optarg);
				break;
			case 'u':
				ulo = atof(optarg);
				break;
			case 'U':
				uhi = atof(optarg);
				break;
				*/
			case 'x':
				divide = atof(optarg);
				break;
			case 'g':
				go_file = optarg;
				break;
			case 'm':
				nmc = atoi(optarg);
				break;
			case 'n':
				nprint = atoi(optarg);
				break;
			case 'o':
				oroot = optarg;
				break;
			case 'r':
				random_seed = atoi(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	ntraj = (argc-optind)/2;	// half dcd, half q
	dcd_names.resize(ntraj);
	q_names.resize(ntraj);
	for (int i=0; i<ntraj; i++) {
		dcd_names[i] = argv[optind+2*i];
		q_names[i] = argv[optind+2*i+1];
	}
	fprintf(stdout,"=========================================================\n");
	fprintf(stdout,"         Go contact list = %s\n", go_file.c_str());
	fprintf(stdout,"Trajectory files:\n");
	nframe=0;
	for (int i=0; i<ntraj; i++) {
		fprintf(stdout,"%s\n",dcd_names[i].c_str());
		inptraj.open(dcd_names[i].c_str());
		nframe += inptraj.total_frames();
		natom = inptraj.num_atoms();
		inptraj.close();
	}
	fprintf(stdout,"Q files:\n");
	for (int i=0; i<ntraj; i++) {
		fprintf(stdout,"%s\n",q_names[i].c_str());
	}
	fprintf(stdout,"Total number of frames = %i\n",nframe);
	fprintf(stdout,"           Random seed = %i\n",random_seed);
	nframe = nframe/skip;
	fprintf(stdout,"# Frames to save = %i (skipping every %i)\n",nframe,skip);
	if  ( ptype == 'm' ) {
		fprintf(stdout,"Using morse potential\n");
	} else if ( ptype == 'k' ) {
		fprintf(stdout,"Using karanicolas&brooks potential\n");
	} 
	fprintf(stdout,"=========================================================\n");
	uniform.initialize(random_seed);
	read_golist(ii, jj, rij, eij, go_file.c_str());
	nc = ii.size();
	simbeta = 1./(8.31*simT/4184.);
	// allocate space
	// ================================
	EIJ = new proj_type[nframe*nc];
	rc = new proj_type[nframe];
	w1 = new proj_type[nframe];
	w2 = new proj_type[nframe];
	// ================================
	w = w1;
	trial_w = w2;
	eij_orig.resize(nc);
	eij_best.resize(nc);
	fprintf(stdout,"=============================================\n");
	fprintf(stdout,"ORIGINAL PARAMETERS:\n");
	fprintf(stdout,"=============================================\n");
	fprintf(stdout,"%4s %4s %12s %12s\n","I","J","RIJ [A]","EIJ kcal/mol");
	// copy original energies and initialize weights
	for (int p=0; p<nframe; p++) {
		w[p] = 1.0;
	}
	for (int p=0; p<nc; p++) {
		eij_orig[p] = eij[p];
		fprintf(stdout,"%4i %4i %12.6f %12.6f\n",ii[p],jj[p],rij[p],eij_orig[p]);
	}
	calc_contact_energy(dcd_names, ii, jj, rij, EIJ, ptype);
	//exit(0);
	read_rc(q_names,rc,nframe);
	crit = opt_crit(w, rc, divide, nframe);
	fprintf(stdout,"divide = %12.6f\n",divide);
	fprintf(stdout,"initial value of opt. crit = %12.6f\n",crit);
	//fprintf(stdout,"w[0] = %12.6f; w[10] = %12.6f;\n",w[0],w[10]);
	int move_successes=0;
	double best_crit = 0.;
	int elem_i, elem_j;
	double ei,ej,emin;
	for (int s = 0; s<nmc; s++) {
		// pick a contact
		choose2(uniform,nc,elem_i,elem_j);
		ei = eij[elem_i];
		ej = eij[elem_j];
		emin = ei<ej ? ei : ej;
		double dE = 1.e6;
		while (abs(eij_orig[elem_i]-(eij[elem_i]+dE))/eij_orig[elem_i] > maxchange 
				|| abs(eij_orig[elem_j]-(eij[elem_j]-dE))/eij_orig[elem_j] > maxchange) {
			dE = emin * gaussian.random()*sigchange;
			//fprintf(stdout,"dE = %12.6e\n", dE);
                        //fprintf(stdout,"%12.6e %12.6e\n",(eij[elem_i]+dE)/eij_orig[elem_i], maxchange);
			//fprintf(stdout,"dE: %12.6f %12.6f %12.6f\n", dE, tmpe, eij_orig[ind]);
		}
		fprintf(stdout,"GOT HERE 0\n");
		//update_weights(EIJ,w,trial_w,ind,dE);
		// update weights and evaluate averages
		for (int f=0; f<nframe; f++) {
			trial_w[f] = w[f]*exp(-simbeta*dE*(EIJ[f*nc+elem_i]-EIJ[f*nc+elem_j]));
		}
		//fprintf(stdout,"w[0] = %12.6f; w[10] = %12.6f;\n",w[0],w[10]);
		//fprintf(stdout,"trial_w[0] = %12.6f; trial_w[10] = %12.6f;\n",trial_w[0],trial_w[10]);
		trial_crit = opt_crit(trial_w, rc, divide, nframe);
		fprintf(stdout,"trial_crit  = %12.6f\n", trial_crit);
		if (trial_crit > crit) {
			proj_type *tmp = w;
			w = trial_w;
			trial_w = tmp;
			crit = trial_crit;
			eij[elem_i] += dE;
			eij[elem_j] -= dE;
			move_successes++;
			fprintf(stdout,"accepting, higher crit\n");
			if (crit>best_crit) {
				// save new best energies - write out at end
				for (int t=0; t<nc; t++) {
					eij_best[t] = eij[t];
				}
				fprintf(stdout,"!! New best energies: move # %i  crit = %12.5f\n", 
						s, crit);
			}
		} else if (mcT>1e-20) {
			double rand = uniform.random();
			if (rand < pow(trial_crit/crit,1.0/mcT)) {
				fprintf(stdout,"accepting, lower crit\n");
				proj_type *tmp = w;
				w = trial_w;
				trial_w = tmp;
				crit = trial_crit;
			//	eij[ind] += dE;
				eij[elem_i] += dE;
				eij[elem_j] -= dE;
				move_successes++;
			}
		} 
		/*
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
	fprintf(stdout,"MC run finished; crit = %12.8f\n", crit);
	fprintf(stdout,"%4s  %8s  %8s  %12s\%\n", "move", 
			"attempts", "successes", "% success");
*/

	}

	// form histogram of unweighted and reweighted p(Q) & save
	if (binhi>binlo) {
		vector<long int> hist;
		vector<double> hist_rw;
		hist.resize(nbin);
		hist_rw.resize(nbin);
		for (int b=0; b<nbin; b++) {
			hist[b] = 0L;
			hist_rw[b] = 0.0;
		}
		double dq = (binhi-binlo)/double(nbin);
		double wsum = 0.0;
		for (int f=0;f<nframe;f++) {
			double qq = rc[f];
			int bin = int(floor((qq-binlo)/dq));
			double dE = 0.;
			for (int t=0;t<nc;t++) {
				dE += (eij_best[t]-eij_orig[t])*EIJ[f*nc+t];
			}
			double ww = exp(-simbeta*dE);
			hist[bin]+=1;
			hist_rw[bin]+= ww;
			wsum += ww;
		}
		string histfile = oroot+string("_hist.dat");
		FILE * histout = fopen(histfile.c_str(),"w");
		for (int b=0; b<nbin; b++) {
			fprintf(histout,"%12.6f %12.6f %12.6f\n",binlo+b*dq,
					float(hist[b])/(float(nframe)*dq),
					hist_rw[b]/(wsum*dq));
		}
		fclose(histout);
	}
	
	fprintf(stdout,"nframe = %i\n", nframe);
	fprintf(stdout,"nc = %i\n", nc);
	// write out eij_nu

	string eijfile = oroot+string("_eij.dat");
	FILE * eijout = fopen(eijfile.c_str(),"w");
	for (int t=0;t<nc;t++) {
		fprintf(eijout,"%5i %5i %12.6f %12.6f\n",
				ii[t]+1,jj[t]+1,eij_best[t],eij_best[t]-eij_orig[t]);
	}
	fclose(eijout);

	return 0;
}

