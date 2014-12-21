/*
 * project contact lists onto arbitrary contact vectors
 */

#include <cstdlib>
#include <cstdio>
#include <string>
#include <list>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <iostream>

#include "TrajFile.h"

double dotp(vector<double> &a, vector<double> &b)
{
	double dim = a.size();
	double p = 0.0;
	if (dim != b.size()) {
		cerr << "dotp: vectors of different sizes!" << endl;
	}
	for (int i=0; i<dim; i++) 
		p += a[i]*b[i];
	return p;
}

inline int index(int i, int j) {
	return j>=i ? j*(j-1)/2+i : i*(i-1)/2+j;
}

inline double contact(double dr, double cut) {
	return 1.0/(1.0+exp(5.0*(dr-cut)));
}

void read_active_vec(vector<double> &v,	vector<int> &active_ind, vector<double> &rcut, const char *file) 
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
		tmpf = atof(strtok(NULL," \t"));
		rcut[i] = tmpf;
		//norm += tmpf*tmpf;
		fprintf(stdout,"%6i %12.8e %12.8f\n",active_ind[i],rcut[i]);
	}
	// ensure norm = 1 (since reading from text file could be small errors)
	//norm = sqrt(norm);
	//for (int i=0; i<dim; i++) {
	//	v[i] /= norm;
	//}
	fclose(f);
	return;
}

void get_ij(int ind, int &i, int &j) 
{
	i = int(floor((1.0+sqrt(1.0+8.0*ind))/2.0));
	j = ind - i*(i-1)/2;
	return;
}

double proj_calc_active(vector<double> &c, vector<int> &a, vector<double> &rcut,
		double *X, double *Y, double *Z)
{
	vector<double> p;
	double dx,dy,dz,dr,r,eqrc;
	int ind, i,j;
	int dim = c.size();
	p.resize(dim);
	for (int ind=0; ind<dim; ind++) {
		get_ij(a[ind],i,j);
		dx = X[i]-X[j];
		dy = Y[i]-Y[j];
		dz = Z[i]-Z[j];
		dr = sqrt(dx*dx+dy*dy+dz*dz);
		p[ind] = contact(dr,rcut[ind]);
		//p[ind] = contact(dr,12.0);
	}
	r = dotp(c,p);
	return r;
}

void read_vec(vector<int> &a,vector<double> &v, const char *file) 
{
	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
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
	a.resize(dim);
	// read conts
	f = fopen(file,"r");
	double norm = 0.0;
	for (int i=0; i<dim; i++) {
		fgets(buf,buf_len,f);
		a[i] = atoi(strtok(buf," \t"));
		fgets(buf,buf_len,f);
		v[i] = atof(strtok(NULL," \t"));
		norm += v[i]*v[i];
	}
	// ensure norm = 1 (since reading from text file could be small errors)
	norm = sqrt(norm);
	for (int i=0; i<dim; i++) {
		v[i] /= norm;
	}
	fclose(f);
	return;
}

/*
void read_contacts(const char *cont_file, vector<int> &ii,
		vector<int> &jj, vector<double> &rij)
{
	FILE *f;
	int nc;
	const int buf_len = 1024;
	char buf[buf_len];
	char *tok;

	// count contacts
	f = fopen(cont_file,"r");
	if (f == NULL) {
		fprintf(stderr,"Could not open file %s\n",cont_file);
		exit(1);
	}
	fgets(buf,buf_len,f);
	nc = 0;
	while (feof(f) == 0) {
		nc++;
		fgets(buf,buf_len,f);
	}
	fclose(f);
	// allocate space
	ii.resize(nc); jj.resize(nc); rij.resize(nc);
	// read conts
	f = fopen(cont_file,"r");
	for (int i=0; i<nc; i++) {
		fgets(buf,buf_len,f);
		tok = strtok(buf," \t");
		ii[i] = atoi(tok) - 1;
		tok = strtok(NULL," \t");
		jj[i] = atoi(tok) - 1;
		tok = strtok(NULL," \t");
		rij[i] = atof(tok);
	}
	fclose(f);
}
*/

int main(int argc, char** argv)
{
	const char *Usage = "Usage: constN_proj -v vec.dat -o proj.dat f1.dcd ... fN.dcd\n";
	int c, ind;
	const int BUF_LEN = 1024;
	//char cont_file[BUF_LEN];
	char vec_file[BUF_LEN];
	char outp_file[BUF_LEN];
	double proj, q;
	vector<string> dcd_files;
	vector<double> v, rcut;
	vector<int> a;
	double *X, *Y, *Z;
	DCDITrajFile inptraj;
	int dim, natom;
	FILE *proj_out;
	double cutoff = 12.0;
	// cmd line opts
	while (1) {
		//c=getopt(argc,argv,"hp:c:v:e:o:");
		c=getopt(argc,argv,"hv:o:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage);
				exit(0);
				break;
			case 'v':
				strcpy(vec_file,optarg);
				break;
			case 'o':
				strcpy(outp_file,optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);
		}
	}
	dcd_files.resize(argc-optind);
	for (int i=0; i<argc-optind; i++) {
		dcd_files[i] = argv[optind+i];
	}
	read_active_vec(v,a,rcut,vec_file);
	inptraj.open(dcd_files[0].c_str());
	natom = inptraj.num_atoms();
	inptraj.close();
	X = new double[natom];
	Y = new double[natom];
	Z = new double[natom];
	dim = v.size();
	proj_out = fopen(outp_file,"w");
	if (proj_out == NULL) {
		fprintf(stderr,"Could not open file \"%s\" for output\n",outp_file);
		exit(1);
	}
	int cum_frame = 0;
	for (int f=0; f<dcd_files.size(); f++) {
		inptraj.open(dcd_files[f].c_str());
		int nframe = inptraj.total_frames();
		for (int ff=0; ff<nframe; ff++) {
			cum_frame++;
			if (cum_frame%1000 == 0) { 
				fprintf(stdout,"processing frame # %i\n",cum_frame);
			}
			inptraj.read_frame(X,Y,Z,natom);
			proj = proj_calc_active(v, a, rcut, X, Y, Z);
			fprintf(proj_out, "%12.5e\n", proj);
		}
		inptraj.close();
	}
	fclose(proj_out);
	//
	//
	//
	delete [] X;
	delete [] Y;
	delete [] Z;
	return 0;
}



