/*
 * PhenoBF.h
 *
 *  Created on: Dec 13, 2012
 *      Author: pickrell
 */

#ifndef PHENOBF_H_
#define PHENOBF_H_
#include "CountData.h"
using namespace std;

class PhenoBF{
public:
	PhenoBF();
	PhenoBF(string, string, PhyloPop_params *);
	gsl_rng* r;
	PhyloPop_params * params;
	CountData* counts_control;
	CountData* counts_pheno;
	void initialize();
	double c1, c2, c3;
	vector<double> xa_control;
	vector<double> xa_pheno;
	int ncontrol;
	int npheno;
	double s;
	bool models;
	double single_llk0_control(int);
	double single_llk0_pheno(int);
	double single_llk1(int);
	double llk0();
	double llk1();
	double normal_ldens(double, double, double);

	double beta_prior, branch_prior, s_prior;
	double xa_update_sd, branch_update_sd, s_update_sd;
	int nit;
	int whichpop;
	bool accept(double, double, double, double);
	void single_iteration();
	void update_xa_control(int);
	void update_xa_pheno(int);
	void update_branches();
	void update_s();
	void run_MCMC(int);
	vector<double> stored_c1;
	vector<double> stored_c2;
	vector<double> stored_c3;
	vector<double> stored_s;
	vector<double> stored_llk;

	void print_stored(string);

};


#endif /* PHENOBF_H_ */
