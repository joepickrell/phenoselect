/*
 * test.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: pickrell
 */

#include "PhenoBF.h"
using namespace std;

int main(){

	PhenoBF t("test_countin.gz", "test_pheno.gz");
	//cout << t.normal_ldens(0, 10, 0.01) << "\n";
	cout << t.llk0() << "\n";
	t.c1 = 0.01;
	t.c2 = 0.02;
	t.c3 = 0.03;
	cout << t.llk0() << "\n";

	return 0;

}
