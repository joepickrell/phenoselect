/*
 * CountData.h
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */

#ifndef COUNTDATA_H_
#define COUNTDATA_H_
#include "Settings.hpp"
#include "PhyloPop_params.h"


class CountData{
public:
	CountData( CountData*, vector<string>, gsl_matrix*, PhyloPop_params *, gsl_rng *);
	void set_cov_ran(gsl_matrix * , gsl_rng *);
	CountData(string, PhyloPop_params*);
	CountData(string); //no rescaling (ie. just the allele frequencies)
	void read_counts(string);
	void read_micro_data(string);
	map<string, int> pop2id;
	map<int, string> id2pop;
	map<int, double> mean_hzy;
	map<int, double> mean_ninds;
	map<int, double> mean_var; //for microsats
	map<int, int> id2nsnp;
	double get_freq(int, int);
	vector<vector<pair<int, int> > > allele_counts;
	vector<vector<vector<float> > > micro_lens;
	int npop, nsnp, ef_nsnp;
	string get_pops(); //in Newick format
	vector<string> list_pops(); //simple list
	string get_pop_in_index(int); //return the name of the population at a given index
	gsl_matrix *alfreqs, *scatter, *cov, *cov_var, *cov_var2;
	//gsl_matrix *cov_samp;
	map<string, map<string, vector<double> > > cov_samp;


	void correct_f2s(string, map<string, float>, float);	//correct f2s for mixture

	gsl_matrix *cov_cov;
	void set_alfreqs();
	void set_alfreqs_micro();
	void scale_alfreqs();

	void read_alfreqs(string); //for debugging, can input simulated allele frequencies

	void print_alfreqs(string);

	pair< vector<string>, vector<double> > get_freqs(int);
	pair< vector<string>, vector<double> > get_centered_freqs(int);
	PhyloPop_params* params;
	gsl_matrix *U, *scatter_prime;
	int ne, ne2, ncomp, nblock;
	int ncomp_ef;
	vector<string> rss;
	vector<string> pos;
	vector<string> chr;
	vector<string> a1;
	vector<string> a2;
	void set_ne();
	void set_ne2();
};

extern int rwishart(gsl_rng *, const int, const int, const gsl_matrix *, gsl_matrix *);
#endif /* COUNTDATA_H_ */
