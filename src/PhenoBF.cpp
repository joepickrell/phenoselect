/*
 * PhenoBF.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: pickrell
 */


#include "PhenoBF.h"
using namespace std;

PhenoBF::PhenoBF(string infilecontrol, string infilepheno, PhyloPop_params* p){
	counts_control = new CountData(infilecontrol, p);
	counts_pheno = new CountData(infilepheno, p);
	ncontrol = counts_control->nsnp;
	npheno = counts_pheno->nsnp;
	params = p;
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
	double toreturn;
	double xa = xa_control[which];
	double f1 = counts_control->get_freq(0, which);
	double f2 = counts_control->get_freq(1, which);
	double f3 = counts_control->get_freq(2, which);
	//cout << xa << " "<< f1 << " "<< f2 << " "<< f3 << "\n";
	double l1 = normal_ldens(f1, xa, xa*(1-xa)*c1);
	double l2 =  normal_ldens(f2, xa, xa*(1-xa)*c2);
	double l3 =  normal_ldens(f3, xa, xa*(1-xa)*c3);
	toreturn = l1+l2+l3;
	if (isnan(toreturn)){
		cout << which << " "<< f1<< " "<< f2 << " "<< f3 << " "<< xa << " "<< l1 << " "<< l2 << " "<< l3 << "\n";
	}
	return toreturn;
}

double PhenoBF::single_llk0_pheno(int which){
	double toreturn;
	double xa = xa_pheno[which];
	double f1 = counts_pheno->get_freq(0, which);
	double f2 = counts_pheno->get_freq(1, which);
	double f3 = counts_pheno->get_freq(2, which);

	double l1 =  normal_ldens(f1, xa, xa*(1-xa)*c1);
	double l2 = normal_ldens(f2, xa, xa*(1-xa)*c2);
	double l3 = normal_ldens(f3, xa, xa*(1-xa)*c3);
	toreturn = l1+l2+l3;
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

	//if (whichsnp==0) cout << xa << " " << m1 << " "<< m2 << " "<< m3 << "\n";
	double l1 =  normal_ldens(f1, m1, xa*(1-xa)*c1);
	double l2 = normal_ldens(f2, m2, xa*(1-xa)*c2);
	double l3 = normal_ldens(f3, m3, xa*(1-xa)*c3);
	toreturn += l1+l2+l3;
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
	xa_update_sd = 0.01;
	s_update_sd = 0.01;
	branch_update_sd = 0.1;
	whichpop = 3;
	nit = 0;
	models = true;
	for (int i = 0; i <ncontrol; i++) {
		double f1 = counts_control->get_freq(0, i);
		double f2 = counts_control->get_freq(1, i);
		double f3 = counts_control->get_freq(2, i);
		double m = (f1+f2+f3)/3.0;
		if (m < 1e-16) m = 1e-5;
		if (m > 1-1e-16)m = 1-1e-5;
		xa_control.push_back(m);
	}
	for (int i = 0; i < npheno; i++) {
		double f1 = counts_pheno->get_freq(0, i);
		double f2 = counts_pheno->get_freq(1, i);
		double f3 = counts_pheno->get_freq(2, i);
		double m = (f1+f2+f3)/3.0;
		if (m < 1e-16) m = 1e-5;
		if (m > 1-1e-16)m = 1-1e-5;
		xa_pheno.push_back( m );
	}
}

double PhenoBF::normal_ldens(double x, double mu, double sigma2){
	double toreturn;
	toreturn = -log(sqrt(sigma2));
	toreturn += -log(sqrt( 2* M_PI));
	toreturn += - ( (x-mu)*(x-mu) )/(2*sigma2);

	return toreturn;
}

void PhenoBF::single_iteration(){
	for (int i = 0; i < ncontrol; i++)	update_xa_control(i);
	for (int i = 0; i < npheno; i++)	update_xa_pheno(i);
	update_branches();
	if (models) update_s();
	nit++;
}

void PhenoBF::update_xa_control(int whichsnp){
	// new x_a = old_xa + N(0, sd)

	double previous = xa_control[whichsnp];
	double lk_previous = single_llk0_control(whichsnp);
	double toadd = gsl_ran_gaussian(r, xa_update_sd);
	double n = previous+toadd;

	// reflect on boundaries
	if (n > 1) n =2-n;
	if (n < 0) n = -n;
	if (n > 1 || n < 0){
		xa_control[whichsnp] = previous;
		return;
	}

	// update
	xa_control[whichsnp] = n;
	double lk_new = single_llk0_control(whichsnp);
	bool a = accept(lk_new, log(gsl_ran_beta_pdf(n, beta_prior, beta_prior)), lk_previous,  log(gsl_ran_beta_pdf(previous, beta_prior, beta_prior)));
	if (!a) xa_control[whichsnp] = previous;
}

void PhenoBF::update_branches(){
	double previous = log(c1);
	double lk_previous = llk1();
	double toadd = gsl_ran_gaussian(r, branch_update_sd);
	double n = previous+toadd;
	c1 = exp(n);
	double lk_new = llk1();
	//cout << exp(previous)<< " "<< c1 << " c1\n";
	bool a = accept(lk_new, log(gsl_ran_gaussian_pdf(n, branch_prior)), lk_previous,  log(gsl_ran_gaussian_pdf(previous, branch_prior)));
	if (!a) {
		c1 = exp(previous);
		lk_new = lk_previous;
	}

	//c2
	previous = log(c2);
	lk_previous = lk_new;
	toadd = gsl_ran_gaussian(r, branch_update_sd);
	n = previous+toadd;
	c2 = exp(n);
	lk_new = llk1();
	//cout << exp(previous)<< " "<< c2 << " c2\n";
	a = accept(lk_new, log(gsl_ran_gaussian_pdf(n, branch_prior)), lk_previous,  log(gsl_ran_gaussian_pdf(previous, branch_prior)));
	if (!a) {
		c2 = exp(previous);
		lk_new = lk_previous;
	}

	//c3
	previous = log(c3);
	lk_previous = lk_new;
	toadd = gsl_ran_gaussian(r, branch_update_sd);
	n = previous+toadd;
	c3 = exp(n);
	lk_new = llk1();
	//cout << exp(previous)<< " "<< c3 << " c3\n";
	a = accept(lk_new, log(gsl_ran_gaussian_pdf(n, branch_prior)), lk_previous,  log(gsl_ran_gaussian_pdf(previous, branch_prior)));
	if (!a) {
		c3 = exp(previous);
		lk_new = lk_previous;
	}

}

bool PhenoBF::accept(double lk1, double p1, double lk0, double p0){
	double logratio = p1+lk1 - p0 -lk0;
	if (logratio > 0) return true;
	double prob = exp(logratio);
	//cout << "prob : "<< prob << "\n";
	double runif = gsl_rng_uniform(r);
	if (runif < prob) return true;
	else return false;
}

void PhenoBF::update_xa_pheno(int whichsnp){
	double previous = xa_pheno[whichsnp];
	double lk_previous = single_llk1(whichsnp);
	double toadd = gsl_ran_gaussian(r, xa_update_sd);
	double n = previous+toadd;
	if (n > 1 || n <0){
		xa_pheno[whichsnp] = previous;
		return;
	}
	xa_pheno[whichsnp] = n;
	double lk_new = single_llk1(whichsnp);
	bool a = accept(lk_new, log(gsl_ran_beta_pdf(n, beta_prior, beta_prior)), lk_previous,  log(gsl_ran_beta_pdf(previous, beta_prior, beta_prior)));
	//cout << a << "\n";
	if (!a) xa_pheno[whichsnp] = previous;
}

void PhenoBF::update_s(){
	double previous = s;
	double lk_previous = llk1();
	double toadd = gsl_ran_gaussian(r, s_update_sd);
	double n = previous+toadd;
	s = n;
	double lk_new = llk1();
	//cout << previous<< " "<< s << " s\n";
	bool a = accept(lk_new, log(gsl_ran_gaussian_pdf(n, s_prior)), lk_previous,  log(gsl_ran_gaussian_pdf(previous, s_prior)));
	//cout << a << "\n";
	if (!a) s = previous;

}

void PhenoBF::run_MCMC(int totalit){
	stored_c1.clear();
	stored_c2.clear();
	stored_c3.clear();
	stored_s.clear();
	stored_llk.clear();
	for (int i = 0; i < totalit; i++){
		single_iteration();
		if (i % params->sfreq == 0){
			double lk1 = llk1();
			cout << i << " "<< c1 << " "<< c2 << " "<< c3 << " "<< s << " "<< lk1 << "\n";
			stored_c1.push_back(c1);
			stored_c2.push_back(c2);
			stored_c3.push_back(c3);
			stored_s.push_back(s);
			stored_llk.push_back( lk1 );
		}
	}
}

void PhenoBF::print_stored(string outfile){
	ogzstream out(outfile.c_str());
	vector<double> sorted_c1, sorted_c2, sorted_c3, sorted_s;
	for (int i = 0; i < stored_c1.size(); i++){
		if (i * params->sfreq < params->nburn){
			out << "# "<< i* params->sfreq;
		}
		else{
			out << i*params->sfreq;
			sorted_c1.push_back(stored_c1[i]);
			sorted_c2.push_back(stored_c2[i]);
			sorted_c3.push_back(stored_c3[i]);
			sorted_s.push_back(stored_s[i]);
		}
		out << " "<< stored_c1[i] << " "<< stored_c2[i] << " "<< stored_c3[i] << " "<< stored_s[i] << " "<< stored_llk[i] << "\n";
	}
	sort(sorted_c1.begin(), sorted_c1.end());
	sort(sorted_c2.begin(), sorted_c2.end());
	sort(sorted_c3.begin(), sorted_c3.end());
	sort(sorted_s.begin(), sorted_s.end());
	double c1_05 = sorted_c1[sorted_c1.size()*0.025];
	double c1_5 = sorted_c1[sorted_c1.size()*0.5];
	double c1_95 = sorted_c1[sorted_c1.size()*0.975];

	double c2_05 = sorted_c2[sorted_c2.size()*0.025];
	double c2_5 = sorted_c2[sorted_c2.size()*0.5];
	double c2_95 = sorted_c2[sorted_c2.size()*0.975];

	double c3_05 = sorted_c3[sorted_c3.size()*0.025];
	double c3_5 = sorted_c3[sorted_c3.size()*0.5];
	double c3_95 = sorted_c3[sorted_c3.size()*0.975];

	double s_05 = sorted_s[sorted_s.size()*0.025];
	double s_5 = sorted_s[sorted_s.size()*0.5];
	double s_95 = sorted_s[sorted_s.size()*0.975];

	double sum = accumulate(sorted_s.begin(), sorted_s.end(), 0.0);
	double mean = sum / sorted_s.size();

	vector<double> diff(sorted_s.size());
	transform(sorted_s.begin(), sorted_s.end(), diff.begin(),
	               bind2nd(minus<double>(), mean));

	double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = sqrt(sq_sum / (sorted_s.size()-1));
	out << "# c1 "<< c1_05 << " "<< c1_5 << " "<< c1_95 << "\n";
	out << "# c2 "<< c2_05 << " "<< c2_5 << " "<< c2_95 << "\n";
	out << "# c3 "<< c3_05 << " "<< c3_5 << " "<< c3_95 << "\n";
	out << "# s "<< s_05 << " "<< s_5 << " "<< s_95 << "\n";
	out << "# s-stats (mean/sd/Z): "<< mean << " "<< stdev << " "<< mean/stdev << "\n";
}
