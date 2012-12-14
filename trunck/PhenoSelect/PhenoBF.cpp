/*
 * PhenoBF.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: pickrell
 */


#include "PhenoBF.h"
using namespace std;

PhenoBF::PhenoBF(string infilecontrol, string infilepheno){
	counts_control = new CountData(infilecontrol);
	counts_pheno = new CountData(infilepheno);
	ncontrol = counts_control->nsnp;
	npheno = counts_pheno->nsnp;
	initialize();
    const gsl_rng_type * T;

    gsl_rng_env_setup();
    T = gsl_rng_ranlxs2;
    r = gsl_rng_alloc(T);
    int seed = (int) time(0);
    gsl_rng_set(r, seed);
    cout << "SEED: "<< seed << "\n";
}

double PhenoBF::single_llk0_control(int which){
	double toreturn = 0;
	double xa = xa_control[which];
	double f1 = counts_control->get_freq(0, which);
	double f2 = counts_control->get_freq(1, which);
	double f3 = counts_control->get_freq(2, which);
	//cout << xa << " "<< f1 << " "<< f2 << " "<< f3 << "\n";
	double l1 = normal_ldens(f1, xa, xa*(1-xa)*c1);
	double l2 =  normal_ldens(f2, xa, xa*(1-xa)*c2);
	double l3 =  normal_ldens(f3, xa, xa*(1-xa)*c3);
	toreturn += l1+l2+l3;
	//cout << l1 << " "<< l2 << " "<< l3 << "\n";
	return toreturn;
}

double PhenoBF::single_llk0_pheno(int which){
	double toreturn = 0;
	double xa = xa_pheno[which];
	double f1 = counts_pheno->get_freq(0, which);
	double f2 = counts_pheno->get_freq(1, which);
	double f3 = counts_pheno->get_freq(2, which);
	double l1 =  normal_ldens(f1, xa, xa*(1-xa)*c1);
	double l2 = normal_ldens(f2, xa, xa*(1-xa)*c2);
	double l3 = normal_ldens(f3, xa, xa*(1-xa)*c3);

	return toreturn;
}

double PhenoBF::single_llk1(int whichsnp){
	double toreturn = 0;
	double xa = xa_pheno[whichsnp];
	double f1 = counts_pheno->get_freq(0, whichsnp);
	double f2 = counts_pheno->get_freq(1, whichsnp);
	double f3 = counts_pheno->get_freq(2, whichsnp);

	double m1 = xa;
	double m2 = xa;
	double m3 = xa;
	if (whichpop == 1) m1 += xa*(1-xa)*s;
	else if (whichpop ==2) m2 += xa*(1-xa)*s;
	else if (whichpop ==3) m3 += xa*(1-xa)*s;
	double l1 =  normal_ldens(f1, m1, xa*(1-xa)*c1);
	double l2 = normal_ldens(f2, m2, xa*(1-xa)*c2);
	double l3 = normal_ldens(f3, m3, xa*(1-xa)*c3);

	return toreturn;
}


double PhenoBF::llk0(){
	double toreturn = 0;
	for (int i = 0; i < ncontrol; i ++) {
		//cout << i << " "<< single_llk0_control(i) << "\n";
		toreturn += single_llk0_control(i);
	}
	for (int i = 0; i < npheno; i++) toreturn += single_llk0_pheno(i);
	return toreturn;
}

double PhenoBF::llk1(){
	double toreturn = 0;
	for (int i = 0; i < ncontrol; i ++) {
		//cout << i << " "<< single_llk0_control(i) << "\n";
		toreturn += single_llk0_control(i);
	}
	for (int i = 0; i < npheno; i++) toreturn += single_llk1(i);
	return toreturn;
}


void PhenoBF::initialize(){
	xa_control.clear();
	xa_pheno.clear();
	s = 0;
	c1 = 0.01;
	c2 = 0.01;
	c3 = 0.01;
	beta_prior = 0.8;
	branch_prior = 10;
	s_prior = 0.1;
	nit = 0;

	for (int i = 0; i <ncontrol; i++) {
		double f1 = counts_control->get_freq(0, i);
		double f2 = counts_control->get_freq(1, i);
		double f3 = counts_control->get_freq(2, i);
		double m = (f1+f2+f3)/3.0;
		if (m < 1e-16) m = 1e-5;
		xa_control.push_back(m);
	}
	for (int i = 0; i < npheno; i++) {
		double f1 = counts_pheno->get_freq(0, i);
		double f2 = counts_pheno->get_freq(1, i);
		double f3 = counts_pheno->get_freq(2, i);
		double m = (f1+f2+f3)/3.0;
		if (m < 1e-16) m = 1e-5;
		xa_pheno.push_back( m );
	}
}

double PhenoBF::normal_ldens(double x, double mu, double sigma2){
	double toreturn =0;
	toreturn += -log(sqrt(sigma2));
	toreturn += -log(sqrt( 2* M_PI));
	toreturn += - ( (x-mu)*(x-mu) )/(2*sigma2);
	return toreturn;
}

void PhenoBF::single_iteration(){
	for (int i = 0; i < ncontrol; i++)	update_xa_control(i);
	for (int i = 0; i < npheno; i++)	update_xa_pheno(i);
	update_branches();
	update_s();
	nit++;
}

void PhenoBF::update_xa_control(int whichsnp){
	double previous = xa_control[whichsnp];
	double lk_previous = single_llk0_control(whichsnp);
	double toadd = gsl_ran_gaussian(r, xa_update_sd);
	double n = previous+toadd;
	xa_control[whichsnp] = n;
	double lk_new = single_llk0_control(whichsnp);
	double logratio = log(gsl_ran_beta_pdf(n, beta_prior, beta_prior))+ lk_new -  log(gsl_ran_beta_pdf(previous, beta_prior, beta_prior))-lk_previous;
}

void PhenoBF::update_branches(){

}

void PhenoBF::update_xa_pheno(int whichsnp){

}
void PhenoBF::update_s(){

}
