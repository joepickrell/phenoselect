/*
 * test.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: pickrell
 */

#include "PhenoBF.h"
using namespace std;

int main(){
	PhyloPop_params p;
	p.pops[0] = "TSI";
	p.pops[1] = "JPT";
	p.pops[2] = "CEU";
	ofstream out("outtest");
	PhenoBF t("../testdata/random_ctj_forall.gz", "../testdata/height2_ctj.gz", &p);
	for (double f = -0.1; f< 1.1; f+= 0.01){
		out << f<< " "<< t.normal_ldens(f, 0.9, 0.5) << "\n";
	}
	//cout << t.normal_ldens(0, 10, 0.01) << "\n";
	//cout << t.llk0() << "\n";
	//t.c1 = 0.01;
	//t.c2 = 0.02;
//	t.c3 = 0.03;
	//cout << t.llk0() << "\n";
	//t.run_MCMC(1000000);
	//for (int i = 0; i < t.stored_c1.size(); i++){
	//	cout << t.stored_c1[i]<< " "<< t.stored_c2[i] << " "<< t.stored_c3[i]<< " " << t.stored_s[i]<< " "<<  t.stored_llk[i]<< "\n";
	//}
	//cout << t.llk0() << "\n";
	return 0;

}
