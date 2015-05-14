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
	CountData(string, PhyloPop_params*);
	PhyloPop_params* params;
	void read_counts(string);
	map<string, int> head2id;
	map<int, string> id2head;
	map<int, double> mean_hzy;
	map<int, double> mean_ninds;
	map<int, int> id2nsnp;
	double get_freq(int, int);
	pair<int, int> get_counts(int, int);
	vector<vector<pair<int, int> > > allele_counts;
	int npop, nsnp;
	string get_pop_in_index(int); //return the name of the population at a given index

};

#endif /* COUNTDATA_H_ */
