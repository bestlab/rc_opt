/* 
 *
 * Find largest contiguous sequence of "native" structure
 *
 * either (i) just longest sequence
 *     or (ii) sequence with largest associated native contact
 *             energy/weight
 *
 * output (i) total length or energy/weight of native sequence
 *        (ii) first residue in sequence
 *        (iii) last residue in sequence
 *
 * Compute a general contact based reaction coordinate
 * based on the following definition of Q:
 *
 * Q = sum_(i,j) w_ij H(r_ij,rcut_ij)
 *
 * where w_ij is a weight factor
 * H(r,rc) is a step function of form H(r,rc) = 1.0/(1.0+exp(beta*(r-rc)))
 * r_ij is the distance between atoms i and j
 * rcut_ij is a cutoff distance specific to the pair (i,j)
 *
 * sum is over contacts in longest nlike seq.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TrajFile.h"

#define DEBUG

string usage = "\n\n	Usage:\n		qcontig go.qlist ncout.dat go1.dcd ... goN.dcd\n	where\n		* go.qlist is the input native contact list (pairs of\n		  atom numbers and 'native' distances)\n		* ncout.dat is the file to write the fraction of native\n		  contacts during simulations\n		  including 'native' distances\n		* goi.dcd are the simulation trajectories\n\n";

void ReadContacts(const string &file, vector<int> &i, vector<int> &j,
		vector<double> &rmin_ij, vector<double> &w_ij)
{
	ifstream inp(file.c_str());
	const int bufsz = 1024;
	char buf[bufsz];
	int p,q,nc;
	nc = 0;
	while (inp.good()) {
		inp.getline(buf,bufsz,'\n');
		nc++;
	}
	nc--;
	inp.close();
	i.resize(nc);
	j.resize(nc);
	rmin_ij.resize(nc);
	w_ij.resize(nc);
	ifstream inp2(file.c_str());
	for (int t=0; t<nc; t++) {
		inp2 >> i[t] >> j[t] >> rmin_ij[t] >> w_ij[t];
	}
	inp2.close();
	return;
}

void purge_local(vector<int> &i, vector<int> &j, vector<double> &rij, 
		vector<double> &wij, int min) 
{
	int nlocal = 0;
	for (int t=0; t<i.size();t++) {
		if ((j[t]-i[t]) >= min) {
			i[nlocal] = i[t];
			j[nlocal] = j[t];
			rij[nlocal] = rij[t];
			wij[nlocal] = wij[t];
			nlocal++;
		}
	}
	i.resize(nlocal);
	j.resize(nlocal);
	rij.resize(nlocal);
	wij.resize(nlocal);
	return;
}

void qene(vector<int> &i, vector<int> &j,vector<double> &rmin_ij,
		vector<double> &w_ij,vector<double> &DR,
		int lo, int hi, double &qtot, double &qss)
{
	int ii, jj;
	double eij,s,s2,s4;
	qtot = 0.0;
	qss = 0.0;
	for (int t=0; t<i.size(); t++) {
		ii = i[t];
		jj = j[t];
		s=rmin_ij[t]/DR[t];
		s2=s*s;
		s4=s2*s2;
		eij = w_ij[t]*(13.*s4*s4*s4-18.*s4*s4*s2+4.*s4*s2);
		qtot += eij;
		if ((ii>=lo && ii<hi)&&(jj>=lo&&jj<hi)) 
			qss += eij;
	}
	return;
}

int find_contig(double *QQ, double *QTOT, int natom, int nexti,
		int &first, int &last, double qcrit=0.9, double mintot = 1.1)
{
	//fprintf(stdout,"GOT HERE 3\n");
	double qq, qtot;
	int i=nexti;
	while (i<natom) { 
		qq = QQ[i];
		qtot = QTOT[i];
		i++;
		if (qq < qcrit && qtot > mintot)
			break;
		//fprintf(stdout,"%i\n",i);
	}
	first = nexti;
	last = i; 	// actually last is last+1 !
	return i;
}

int main( int argc, char ** argv )
{
	string qlist, output;
	vector<string> dcd_names;
	int nselect,ncontact,natom,ii,jj;
	int * selection;
	vector<int> i, j;
	vector<double> rcut_ij, w_ij, Eij, Qij, DR, rmin_ij;
	double com[3], mat[3][3], dr,dx,dy,dz,beta;
	int minseq;
	double qcrit = 0.5;
	double fudge_lambda = 1.2;
	double s,s2,s4;
	int lfirst, llast, local;

	beta = 5.0;
	if (argc < 4) {
		cout << usage << endl;
		exit(0);
	}
	local = atoi(argv[1]);
	qcrit = atof(argv[2]);
	minseq = atoi(argv[3]);
	qlist = string(argv[4]);
	output = string(argv[5]);
	dcd_names.resize(argc-6);
	for (int t=0; t<dcd_names.size(); t++) {
		dcd_names[t] = string(argv[6+t]);
	}
	ReadContacts(qlist,i,j,rmin_ij,w_ij);
	fprintf(stdout,"ncontact = %i\n", i.size());
	purge_local(i,j,rmin_ij, w_ij,local);
	fprintf(stdout,"nonlocal = %i\n", i.size());
	DCDITrajFile inptraj;
	inptraj.open(dcd_names[0].c_str());
	natom = inptraj.num_atoms();
	ncontact = i.size();
	inptraj.close();
	//double *QQ = new double[ncontact];
	double *X = new double[natom];
	double *Y = new double[natom];
	double *Z = new double[natom];
	Eij.resize(ncontact);
	Qij.resize(ncontact);
	int frames = 0;
	FILE *outp = fopen(output.c_str(),"w");
	if ( outp == NULL ) {
		fprintf(stderr,"Could not open file %s\n",output.c_str());
		exit(1);
	}
	rcut_ij.resize(rmin_ij.size());
	for (int t=0; t<rmin_ij.size(); t++) 
		rcut_ij[t] = rmin_ij[t]*fudge_lambda;
	//ofstream outp(output.c_str());
	double ess, etot, qss;
	for (int trajfile=0; trajfile<dcd_names.size(); trajfile++) {
		// iterate over frames ...
		inptraj.open( dcd_names[trajfile].c_str() );
		int no_frames = inptraj.total_frames();
		for (int frame = 0; frame < no_frames; frame++) {
			inptraj.read_frame(X,Y,Z,natom);
			etot = 0.0;
			for (int t=0; t<i.size(); t++) {
				ii = i[t]-1;
				jj = j[t]-1;
				dx = X[ii]-X[jj];
				dy = Y[ii]-Y[jj];
				dz = Z[ii]-Z[jj];
				dr = sqrt(dx*dx+dy*dy+dz*dz);
				s=rmin_ij[t]/dr;
				s2=s*s;
				s4=s2*s2;
				double q = w_ij[t]/(1.0+exp(beta*(dr-fudge_lambda*rcut_ij[t])));
				double eij = w_ij[t]*(13.*s4*s4*s4-18.*s4*s4*s2+4.*s4*s2);
				Eij[t] = eij;
				Qij[t] = q;
				etot += eij;
				//fprintf(stdout,"---- %8.3f %8.3f %8.3f\n",dr,rmin_ij[t], eij);
				//DR[t] = dr;
				//QQ[t] = q;
			}
			double esingle, qsingle;
			int max_sslen = 0;
			lfirst = llast=0;
			for (int first=0; first<(natom-minseq); first++) {
				for (int last=(first+minseq); last<natom; last++) {
					double wsum = 0.0;
					int ncon = 0;
					qss=0.0;
					ess=0.0;
					for (int t=0; t<ncontact; t++) {
						ii=i[t]-1;
						jj=j[t]-1;
						if (ii<first || ii>=last)
							continue;
						if (jj<first || jj>=last)
							continue;
						ess += Eij[t];
						qss += Qij[t];
						wsum += w_ij[t];
						ncon += 1;
					}
					if (ncon>0) {
						qss/=wsum;
					} else{ 
						qss=0.0;
					}
					//fprintf(stdout,"%8.3f %8.3f %8.3f\n",qss,ess,etot);
					if (qss > qcrit && (last-first)>max_sslen) {
						max_sslen = last-first;
						lfirst = first;
						llast = last;
						esingle = ess;
						qsingle = qss;
					}
				}
			}

			/* OLD ALGORITHM (AS PRESENTED AT CECAM)
			while (nexti<natom) {
				if (QQ[nexti] > qcrit || QTOT[nexti]<mintot) {
					//fprintf(stdout,"GOT HERE 2\n");
					nexti = find_contig(QQ,QTOT,natom,\
							nexti,first,last,qcrit,mintot);
					//fprintf(stdout,"%i %i %i\n",
					//		nexti,first,last);
					if ((last-first)>longest) {
						longest = double(last-first);
						lfirst = first;
						llast = last;
					}
				} else { 
					nexti++;
				}
			}
			*/


			frames++;
			//fprintf(stdout,"GOT HERE 3.5\n");
			fprintf(outp,"%5i %5i %5i %8.3f %8.3f %8.3f\n",max_sslen,lfirst,llast,esingle,etot,qsingle);
			//fprintf(outp,"%5i %5i\n",first,last);
			//fprintf(stdout,"GOT HERE 4\n");
		}
		inptraj.close();
	}

	//delete [] QQ;
	//delete [] QTOT;
	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}

